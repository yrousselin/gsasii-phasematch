# -*- coding: utf-8 -*-
"""Recherche de phases à l'aveugle (search-match) dans GSAS-II — fenêtre interactive.

Menu : Export -> Entire project as -> (user) Phase search (search-match)

Le programme classe les références de l'index (CSD locale, dossiers de CIF) d'après le
diagramme ouvert ; l'utilisateur regarde la superposition, retient une phase, et la recherche
suivante porte sur ce qui reste. Fichiers nécessaires dans le même dossier : pm_session.py,
pm_search.py, pm_index.py, pm_core.py, pm_f0.py, pm_gpx.py, pm_plot.py, csd_pm_worker.py.
"""
import os
import sys
import threading
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

try:
    from .. import GSASIIfiles as G2fil
except Exception:
    try:
        from GSASII import GSASIIfiles as G2fil
    except Exception:
        import GSASIIfiles as G2fil
import wx

DEFAULT_DB = os.path.expanduser("~/.GSASII/pm_index.sqlite")


def tree_patterns(G2frame):
    """Histogrammes de poudre de l'arbre du projet ouvert (sans sauvegarde). -> {nom: (data, sous-éléments)}"""
    tree = G2frame.GPXtree
    out = {}
    item, cookie = tree.GetFirstChild(G2frame.root)
    while item:
        name = tree.GetItemText(item)
        if name.startswith("PWDR"):
            sub = {}
            item2, cookie2 = tree.GetFirstChild(item)
            while item2:
                sub[tree.GetItemText(item2)] = tree.GetItemPyData(item2)
                item2, cookie2 = tree.GetNextChild(item, cookie2)
            out[name] = (tree.GetItemPyData(item)[1], sub)
        item, cookie = tree.GetNextChild(G2frame.root, cookie)
    return out


class PhaseMatchFrame(wx.Frame):
    def __init__(self, parent, pat, db, outdir, project=""):
        wx.Frame.__init__(self, parent, title="Recherche de phases — " + pat["name"][:80], size=(1500, 880))
        import pm_session
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, NavigationToolbar2WxAgg
        self.session = pm_session.Session(pat, db)
        self.outdir = outdir
        self.project = project
        self._g2frame = parent          # pour retrouver le chemin du projet s'il est enregistré après coup
        # résultats relus par le rapport d'affinement : <projet>_phases.json, tenu à jour à chaque phase retenue ou retirée
        self.results_path = (os.path.splitext(project)[0] + "_phases.json") if project else None
        self._wrote = False
        self.cands = []
        self.busy = False
        self._thread = None
        self.CreateStatusBar()

        split = wx.SplitterWindow(self, style=wx.SP_LIVE_UPDATE)
        left = wx.Panel(split)
        right = wx.Panel(split)

        # ---- colonne de gauche ----
        v = wx.BoxSizer(wx.VERTICAL)
        v.Add(wx.StaticText(left, label="Phases retenues"), 0, wx.ALL, 4)
        self.acc_list = wx.ListCtrl(left, style=wx.LC_REPORT | wx.LC_SINGLE_SEL, size=(-1, 120))
        for i, (lab, w) in enumerate((("n°", 30), ("Identifiant", 110), ("Nom", 200), ("Part %", 60), ("≈ % mass.", 75))):
            self.acc_list.InsertColumn(i, lab, width=w)
        v.Add(self.acc_list, 0, wx.EXPAND | wx.LEFT | wx.RIGHT, 4)
        h = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_remove = wx.Button(left, label="Retirer la phase sélectionnée")
        self.btn_export = wx.Button(left, label="Exporter les phases retenues…")
        h.Add(self.btn_remove, 0, wx.RIGHT, 6)
        h.Add(self.btn_export, 0)
        v.Add(h, 0, wx.ALL, 4)

        self.cand_label = wx.StaticText(left, label="Candidats")
        v.Add(self.cand_label, 0, wx.ALL, 4)
        self.cand_list = wx.ListCtrl(left, style=wx.LC_REPORT | wx.LC_SINGLE_SEL)
        for i, (lab, w) in enumerate((("n°", 30), ("Identifiant", 100), ("Nom", 170), ("Validé %", 62), ("Couv. %", 58), ("Contraste", 64), ("T (K)", 48), ("Équiv.", 48))):
            self.cand_list.InsertColumn(i, lab, width=w)
        v.Add(self.cand_list, 1, wx.EXPAND | wx.LEFT | wx.RIGHT, 4)
        h = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_accept = wx.Button(left, label="Retenir cette phase")
        self.btn_mercury = wx.Button(left, label="Voir dans Mercury")
        h.Add(self.btn_accept, 0, wx.RIGHT, 6)
        h.Add(self.btn_mercury, 0)
        v.Add(h, 0, wx.ALL, 4)

        self.info = wx.TextCtrl(left, style=wx.TE_MULTILINE | wx.TE_READONLY | wx.TE_DONTWRAP, size=(-1, 230))
        self.info.SetFont(wx.Font(11, wx.FONTFAMILY_TELETYPE, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL))
        v.Add(self.info, 0, wx.EXPAND | wx.LEFT | wx.RIGHT, 4)

        h = wx.BoxSizer(wx.HORIZONTAL)
        h.Add(wx.StaticText(left, label="Éléments exigés :"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        self.txt_elements = wx.TextCtrl(left, value="", size=(110, -1))
        self.chk_strict = wx.CheckBox(left, label="aucun autre")
        h.Add(self.txt_elements, 0, wx.RIGHT, 6)
        h.Add(self.chk_strict, 0, wx.ALIGN_CENTER_VERTICAL)
        v.Add(h, 0, wx.ALL, 4)
        h = wx.BoxSizer(wx.HORIZONTAL)
        h.Add(wx.StaticText(left, label="Candidats par tour :"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        self.spin_n = wx.SpinCtrl(left, min=4, max=30, initial=10, size=(60, -1))
        self.btn_search = wx.Button(left, label="Relancer la recherche")
        h.Add(self.spin_n, 0, wx.RIGHT, 8)
        h.Add(self.btn_search, 0)
        v.Add(h, 0, wx.ALL, 4)
        left.SetSizer(v)

        # ---- tracé ----
        self.fig = Figure(figsize=(9, 7), dpi=100)
        self.canvas = FigureCanvasWxAgg(right, -1, self.fig)
        self.toolbar = NavigationToolbar2WxAgg(self.canvas)
        self.toolbar.Realize()
        rv = wx.BoxSizer(wx.VERTICAL)
        rv.Add(self.canvas, 1, wx.EXPAND)
        rv.Add(self.toolbar, 0, wx.EXPAND)
        right.SetSizer(rv)
        split.SplitVertically(left, right, 600)
        split.SetMinimumPaneSize(300)

        self.btn_accept.Bind(wx.EVT_BUTTON, self.on_accept)
        self.btn_remove.Bind(wx.EVT_BUTTON, self.on_remove)
        self.btn_export.Bind(wx.EVT_BUTTON, self.on_export)
        self.btn_mercury.Bind(wx.EVT_BUTTON, self.on_mercury)
        self.btn_search.Bind(wx.EVT_BUTTON, self.on_search)
        self.cand_list.Bind(wx.EVT_LIST_ITEM_SELECTED, self.on_select)
        self.cand_list.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.on_accept)
        self.Bind(wx.EVT_CLOSE, self.on_close)
        self._alive = True
        self.redraw()
        self.start(self._load_and_propose, self.show_candidates)

    # ---- travail en arrière-plan ----------------------------------------------------------
    def status(self, msg):
        if self._alive:
            wx.CallAfter(self._set_status, msg)

    def _set_status(self, msg):
        if self._alive:
            self.SetStatusText(msg)

    def start(self, work, done):
        if self.busy:
            return
        self.busy = True
        self.enable(False)

        def target():
            res, err = None, None
            try:
                res = work()
            except Exception:
                err = traceback.format_exc()
            wx.CallAfter(self._finish, done, res, err)

        self._thread = threading.Thread(target=target, daemon=True)
        self._thread.start()

    def _finish(self, done, res, err):
        if not self._alive:
            return
        self.busy = False
        self.enable(True)
        if err:
            self.SetStatusText("Erreur pendant le calcul.")
            wx.MessageBox(err[-1500:], "Recherche de phases : erreur", wx.OK | wx.ICON_ERROR, self)
            return
        done(res)

    def enable(self, on):
        for b in (self.btn_accept, self.btn_remove, self.btn_export, self.btn_mercury, self.btn_search):
            b.Enable(on)

    def _load_and_propose(self):
        if self.session.ix is None:
            n = self.session.load(self.status)
            self.status("Index : %d références." % n)
        return self.session.propose(self.status)

    def _propose(self):
        return self.session.propose(self.status)

    # ---- affichage ------------------------------------------------------------------------
    def current(self):
        i = self.cand_list.GetFirstSelected()
        return self.cands[i] if 0 <= i < len(self.cands) else None

    def redraw(self, cand=None, keep_zoom=False):
        import pm_plot
        xlim = None
        if keep_zoom and self.fig.axes:
            xlim = self.fig.axes[0].get_xlim()
        pm_plot.draw(self.fig, self.session, cand, keep_xlim=xlim)
        self.canvas.draw()

    def show_accepted(self):
        self.acc_list.DeleteAllItems()
        for i, a in enumerate(self.session.accepted):
            self.acc_list.InsertItem(i, str(i + 1))
            self.acc_list.SetItem(i, 1, a.get("label") or a["ident"])
            self.acc_list.SetItem(i, 2, a["name"])
            self.acc_list.SetItem(i, 3, "%.0f" % (100 * a.get("share", 0)))
            self.acc_list.SetItem(i, 4, "%.0f" % (100 * a["wt"]) if a.get("wt") is not None else "-")
            self.acc_list.SetItemTextColour(i, wx.Colour(a["color"]))

    def show_candidates(self, cands):
        self.cands = cands or []
        self.cand_list.DeleteAllItems()
        for i, c in enumerate(self.cands):
            self.cand_list.InsertItem(i, str(i + 1))
            self.cand_list.SetItem(i, 1, c.get("label") or c["ident"])
            self.cand_list.SetItem(i, 2, c["name"])
            self.cand_list.SetItem(i, 3, "%.0f" % (100 * c["s3"]))
            self.cand_list.SetItem(i, 4, "%.0f" % (100 * c["cover"]))
            self.cand_list.SetItem(i, 5, "%.1f" % c["contrast"])
            self.cand_list.SetItem(i, 6, "%.0f" % c["temp_K"] if c["temp_K"] else "-")
            self.cand_list.SetItem(i, 7, str(len(c["equiv"])))
        self.cand_label.SetLabel("Candidats — tour %d (sur ce qui reste à expliquer)" % self.session.round)
        self.show_accepted()
        if self.cands:
            self.cand_list.Select(0)
            self.on_select(None)              # au cas où la sélection par programme ne déclenche pas l'événement
            self.SetStatusText("Tour %d : %d phases candidates. Sélectionner pour superposer, « Retenir » pour passer à la suite." % (self.session.round, len(self.cands)))
        else:
            self.info.SetValue("")
            self.redraw()
            self.SetStatusText("Aucun candidat : plus assez de pics dans ce qui reste, ou filtre trop strict.")

    def on_select(self, event):
        c = self.current()
        if c is None:
            return
        txt = ["Validé    %5.1f %%  part de ce qui reste expliquée par la structure, maille ajustée" % (100 * c["s3"]),
               "Couverture %4.0f %%  part des raies calculées du candidat présentes dans le diagramme" % (100 * c["cover"]),
               "Contraste  %4.1f    rapport au meilleur des autres candidats (au-dessus de 2 : net)" % c["contrast"],
               "Maille ajustée  %.4f %.4f %.4f  %.3f %.3f %.3f" % tuple(c["cell"]), "", c["info"]]
        if c["equiv"]:
            txt += ["", "Autres entrées de la même phase : " + ", ".join(c["equiv"][:12])]
        self.info.SetValue("\n".join(txt))
        self.redraw(c, keep_zoom=True)

    # ---- résultats pour le rapport ----------------------------------------------------------
    def save_results(self, quiet=True):
        """Écrit l'état courant dans <projet>_phases.json. Si le projet n'avait pas de fichier à l'ouverture de la
        fenêtre mais a été enregistré depuis, le chemin est repris maintenant."""
        if not self.results_path:
            proj = getattr(self._g2frame, "GSASprojectfile", "") or ""
            if proj:
                self.project = proj
                self.results_path = os.path.splitext(proj)[0] + "_phases.json"
        if not self.results_path:
            self.SetStatusText("Projet non enregistré : résultats non disponibles pour le rapport d'affinement "
                               "(File > Save project, puis retenir de nouveau la phase).")
            return False
        try:
            self.session.save_results(self.results_path, self.project)
            self._wrote = True
            if not quiet:
                self.SetStatusText("Résultats enregistrés pour le rapport : " + self.results_path)
            return True
        except Exception as ex:
            self.SetStatusText("Résultats non enregistrés (%s)" % ex)
            return False

    def sync_results(self):
        """Après « Retenir » ou « Retirer » : fichier mis à jour ; s'il ne reste aucune phase, le fichier écrit
        pendant cette séance est retiré (la question « aucune phase trouvée » est posée à la fermeture)."""
        if self.session.accepted:
            self.save_results(quiet=False)
        elif self._wrote and self.results_path and os.path.exists(self.results_path):
            try:
                os.remove(self.results_path)
            except OSError:
                pass
            self._wrote = False

    # ---- actions --------------------------------------------------------------------------
    def on_accept(self, event):
        c = self.current()
        if c is None or self.busy:
            return
        self.session.accept(c)
        self.sync_results()
        self.show_accepted()
        self.cand_list.DeleteAllItems()
        self.cands = []
        self.redraw(keep_zoom=True)
        self.start(self._propose, self.show_candidates)

    def on_remove(self, event):
        i = self.acc_list.GetFirstSelected()
        if i < 0:
            i = len(self.session.accepted) - 1
        if i < 0 or self.busy:
            return
        self.session.remove(i)
        self.sync_results()
        self.show_accepted()
        self.cand_list.DeleteAllItems()
        self.cands = []
        self.redraw(keep_zoom=True)
        self.start(self._propose, self.show_candidates)

    def on_search(self, event):
        if self.busy:
            return
        els = self.txt_elements.GetValue().replace(",", " ").split()
        els = [e[0].upper() + e[1:].lower() for e in els]
        strict = self.chk_strict.GetValue()
        self.session.n3 = int(self.spin_n.GetValue())

        def work():
            if self.session.ix is None:
                self.session.load(self.status)
            n = self.session.set_filter(None, els or None, strict)
            self.status("%d références après filtrage." % n)
            return self.session.propose(self.status)

        self.cand_list.DeleteAllItems()
        self.cands = []
        self.start(work, self.show_candidates)

    def on_mercury(self, event):
        c = self.current()
        if c is None:
            i = self.acc_list.GetFirstSelected()
            c = self.session.accepted[i] if 0 <= i < len(self.session.accepted) else None
        if c is None:
            return
        import pm_index
        try:
            self.SetStatusText(pm_index.open_in_mercury(c["cif"], c["ident"]))
        except Exception as ex:
            wx.MessageBox(str(ex), "Mercury", wx.OK | wx.ICON_WARNING, self)

    def on_export(self, event):
        if not self.session.accepted:
            wx.MessageBox("Aucune phase retenue.", "Export", wx.OK | wx.ICON_INFORMATION, self)
            return
        dlg = wx.DirDialog(self, "Dossier où écrire les CIF des phases retenues", defaultPath=self.outdir if os.path.isdir(os.path.dirname(self.outdir)) else "")
        try:
            if dlg.ShowModal() != wx.ID_OK:
                return
            out = dlg.GetPath()
        finally:
            dlg.Destroy()
        try:
            files = self.session.export(out)
            self.save_results()
        except Exception as ex:
            wx.MessageBox(str(ex), "Export", wx.OK | wx.ICON_ERROR, self)
            return
        wx.MessageBox("%d fichiers CIF (maille ajustée) et phases_retenues.csv écrits dans :\n%s\n\nPour l'affinement : Import > Phase > from CIF file." % (len(files), out),
                      "Export", wx.OK | wx.ICON_INFORMATION, self)

    def on_close(self, event):
        try:
            s = self.session
            if self.results_path and s.searched and not s.accepted and not self.busy:
                old = os.path.exists(self.results_path)
                msg = ("Aucune phase n'a été retenue.\n\nEnregistrer « aucune phase concordante trouvée » pour le rapport d'affinement ?"
                       + ("\n(le résultat d'une recherche précédente sera remplacé)" if old else ""))
                if wx.MessageBox(msg, "Recherche de phases", wx.YES_NO | wx.NO_DEFAULT | wx.ICON_QUESTION, self) == wx.YES:
                    self.save_results()
        except Exception:
            traceback.print_exc()
        self._alive = False
        event.Skip()


class ExportPhaseMatch(G2fil.ExportBaseclass):
    """Ouvre la fenêtre de recherche de phases pour un histogramme de poudre du projet."""

    def __init__(self, G2frame):
        super(self.__class__, self).__init__(G2frame=G2frame, formatName="Phase search (search-match)", extension=".csv",
                                             longFormatName="Recherche de phases à l'aveugle dans l'index CSD (fenêtre interactive)")
        self.exporttype = ["project"]

    def Exporter(self, event=None):
        G2frame = self.G2frame
        try:
            import pm_gpx
            if not os.path.exists(DEFAULT_DB):
                wx.MessageBox("Index introuvable : %s\n\nLe construire une fois dans le Terminal :\n  cd ~/.GSASII/exports\n  caffeinate -i ~/g2main/bin/python pm_index.py build --csd" % DEFAULT_DB,
                              "Recherche de phases", wx.OK | wx.ICON_WARNING, G2frame)
                return
            pats = tree_patterns(G2frame)
            names = [n for n in pats if not n.startswith("PWDR Autobkg")] or list(pats)
            if not names:
                wx.MessageBox("Aucun histogramme de poudre dans le projet.", "Recherche de phases", wx.OK | wx.ICON_WARNING, G2frame)
                return
            name = names[0]
            if len(names) > 1:
                dlg = wx.SingleChoiceDialog(G2frame, "Histogramme à analyser", "Recherche de phases", names)
                try:
                    if dlg.ShowModal() != wx.ID_OK:
                        return
                    name = dlg.GetStringSelection()
                finally:
                    dlg.Destroy()
            data, sub = pats[name]
            pat = pm_gpx.pattern_from_parts(name, data, sub["Instrument Parameters"][0], sub["Limits"][1])
            proj = getattr(G2frame, "GSASprojectfile", "") or ""
            outdir = (os.path.splitext(proj)[0] + "_pm") if proj else os.path.expanduser("~/pm_resultats")
            frame = PhaseMatchFrame(G2frame, pat, DEFAULT_DB, outdir, proj)
            frame.Show()
            if not proj:
                wx.MessageBox("Le projet n'est pas encore enregistré.\n\nVous pouvez chercher les phases, mais pour que le "
                              "résultat soit repris par le rapport d'affinement (<projet>_phases.json), enregistrez le projet "
                              "(File > Save project) avant de retenir les phases - ou retenez-les de nouveau après "
                              "l'enregistrement.", "Recherche de phases", wx.OK | wx.ICON_INFORMATION, frame)
        except Exception:
            wx.MessageBox(traceback.format_exc()[-1500:], "Recherche de phases : erreur", wx.OK | wx.ICON_ERROR, G2frame)
