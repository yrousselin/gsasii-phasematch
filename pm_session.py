"""pm_session : session interactive de recherche de phases (moteur de la fenêtre GSAS-II).

  s = Session(pat, db); s.load(); cands = s.propose(); s.accept(cands[0]); cands = s.propose(); ...
Le programme classe les candidats ; c'est l'utilisateur qui retient une phase. Les phases retenues
sont ajustées ensemble (échelles non négatives) et la recherche suivante porte sur ce qui reste.
"""
import csv, datetime, json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
from scipy.optimize import nnls

import pm_core as pc
import pm_index as pi
import pm_search as ps

COLORS = ["#1f77b4", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b", "#17becf", "#bcbd22", "#e377c2"]
RESULTS_FORMAT = 1

# masses atomiques (g/mol) pour l'estimation des proportions massiques
_MASS = dict(zip(
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd "
    "Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra "
    "Ac Th Pa U Np Pu".split(),
    [1.008, 4.0026, 6.94, 9.0122, 10.81, 12.011, 14.007, 15.999, 18.998, 20.180, 22.990, 24.305, 26.982, 28.085, 30.974, 32.06,
     35.45, 39.948, 39.098, 40.078, 44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.38, 69.723, 72.630,
     74.922, 78.971, 79.904, 83.798, 85.468, 87.62, 88.906, 91.224, 92.906, 95.95, 98.0, 101.07, 102.91, 106.42, 107.87, 112.41,
     114.82, 118.71, 121.76, 127.60, 126.90, 131.29, 132.91, 137.33, 138.91, 140.12, 140.91, 144.24, 145.0, 150.36, 151.96, 157.25,
     158.93, 162.50, 164.93, 167.26, 168.93, 173.05, 174.97, 178.49, 180.95, 183.84, 186.21, 190.23, 192.22, 195.08, 196.97, 200.59,
     204.38, 207.2, 208.98, 209.0, 210.0, 222.0, 223.0, 226.0, 227.0, 232.04, 231.04, 238.03, 237.0, 244.0]))
_MASS["D"] = 2.014


def formula_mass(formula):
    """Masse molaire d'une formule brute simple (« C12 H22 O11 ») ; None si elle n'est pas lisible sans ambiguïté."""
    f = (formula or "").strip()
    if not f or "(" in f or "," in f or "[" in f or f[0].isdigit():
        return None
    tot = 0.0
    for el, num in re.findall(r"([A-Z][a-z]?)(\d*\.?\d*)", re.sub(r"\d*[+-]", "", f)):
        if el not in _MASS:
            return None
        tot += _MASS[el] * (float(num) if num else 1.0)
    return tot or None


def cell_mass(st, formula=None, z=None):
    """Masse du contenu de la maille (g/mol). Somme sur les atomes de la structure ; si la formule × Z est cohérente
    avec cette somme (à 25 % près), c'est elle qui est prise : elle compte les hydrogènes absents du modèle."""
    els, _, occ, _ = pc.expand_atoms(st["ops"], st["atoms"])
    m_atoms = float(sum(_MASS.get(e, 0.0) * o for e, o in zip(els, occ)))
    mf = formula_mass(formula)
    try:
        m_form = mf * float(z) if mf and z else None
    except (TypeError, ValueError):
        m_form = None
    if m_form and m_atoms > 0 and 0.8 < m_atoms / m_form < 1.25:
        return m_form
    return m_atoms


def database_stats(db):
    """Nombre d'entrées de chaque base dans l'index : {source: {"entries": n, "usable": n_ok, "version": ..}}."""
    con = pi.open_db(db)
    out = {}
    for src, n, ok in con.execute("SELECT source, COUNT(*), SUM(status='ok') FROM entries GROUP BY source"):
        out[src] = {"entries": int(n), "usable": int(ok or 0)}
    meta = dict(con.execute("SELECT key, value FROM meta").fetchall())
    if "CSD" in out and meta.get("csd_version"):
        out["CSD"]["version"] = str(meta["csd_version"])
    for src in out:
        last = con.execute("SELECT MAX(added) FROM entries WHERE source=?", (src,)).fetchone()[0]
        if last:
            out[src]["updated"] = str(last)[:10]
    return out


class Session:
    def __init__(self, pat, db, sources=None, elements=None, strict=False, n3=10, n2=400, tol=0.4, max_shift=(-0.6, 1.8)):
        self.pat, self.db = pat, db
        self.sources, self.elements, self.strict = sources, elements, strict
        self.n3, self.n2, self.tol, self.max_shift = n3, n2, tol, max_shift
        self.obs = ps.Obs(pat)
        self.net = np.clip(self.obs.net, 0, None)
        self.accepted = []
        self.round = 1
        self.ix = None
        self.db_stats = {}
        self.searched = False          # une recherche a-t-elle abouti (liste de candidats proposée) ?

    # -- index ------------------------------------------------------------------------------
    def load(self, progress=lambda m: None):
        progress("Lecture de l'index (plus long la première fois)…")
        con = pi.open_db(self.db)
        self.ix = pi.load_index(self.db)
        try:
            self.db_stats = database_stats(self.db)
        except Exception:
            self.db_stats = {}
        self.set_filter(self.sources, self.elements, self.strict, con)
        return len(self.ix["ids"])

    def set_filter(self, sources=None, elements=None, strict=False, con=None):
        self.sources, self.elements, self.strict = sources, elements, strict
        con = con or pi.open_db(self.db)
        self.rows = ps.select_rows(con, self.ix, sources, elements, strict)
        return len(self.rows)

    # -- recherche sur le résidu courant --------------------------------------------------
    def propose(self, progress=lambda m: None):
        obs, ix, rows = self.obs, self.ix, self.rows
        con = pi.open_db(self.db)
        o, thr, pk_bins, pk_h = obs.grid()
        if len(pk_bins) < 4 or len(rows) == 0:
            self.searched = True
            return []
        tol_bins = int(round(self.tol / 100 / pc.LOG_STEP))
        shifts = np.arange(int(self.max_shift[0] / 100 / pc.LOG_STEP), int(self.max_shift[1] / 100 / pc.LOG_STEP) + 1, max(1, tol_bins // 2))
        progress("Tour %d : tri de %d références…" % (self.round, len(rows)))
        s1, sh = ps.stage1(ix, rows, obs, o, thr, tol_bins, shifts)
        order = np.argsort(-s1)[:self.n2]
        order = order[s1[order] > 0]
        progress("Tour %d : tri fin de %d candidats…" % (self.round, len(order)))
        res2 = ps.stage2(ix, rows[order], sh[order], obs, o, pk_bins)
        sc2 = np.array([r[0] for r in res2]) if res2 else np.zeros(0)
        o2 = np.argsort(-sc2)
        taken = {a["id"] for a in self.accepted}
        leaders, fam, chosen = [], {}, []
        for k in o2[:80]:
            g = next((i for i, c in enumerate(leaders) if ps.same_lines(c, res2[k][6])), None)
            if g is None:
                leaders.append(res2[k][6])
                g = len(leaders) - 1
            fam.setdefault(g, []).append(k)
            if len(fam[g]) <= (2 if g == 0 else 1) and len(chosen) < self.n3 and int(ix["ids"][rows[order[k]]]) not in taken:
                chosen.append((k, g))
        dbrows = {k: dict(pi.get_row(con, int(ix["ids"][rows[order[k]]]))) for k, g in chosen}
        progress("Tour %d : lecture de %d structures…" % (self.round, len(chosen)))
        cifs = pi.get_cifs(con, list(dbrows.values()))
        cands = []
        ss = float(obs.res @ obs.res) or 1.0
        for n, (k, g) in enumerate(chosen):
            row = dbrows[k]
            progress("Tour %d : validation %d/%d (%s)…" % (self.round, n + 1, len(chosen), row["ident"]))
            r = rows[order[k]]
            cell2 = ps.A_to_cell(res2[k][3] @ res2[k][2]) or list(ix["cell"][r])
            try:
                txt = cifs.get(row["id"], "")
                if not txt:
                    raise ValueError("structure illisible")
                st = pc.structure_from_cif(txt)
                fit = ps.Fit3(obs, st, cell2)
                expl, ycalc, cell3, wfit = fit.run()
                stt, sti = fit.sticks
                sabs = float(fit.s) / fit.imax
                try:
                    mass = cell_mass(st, row["formula"], json.loads(row["info"] or "{}").get("z"))
                except Exception:
                    mass = 0.0
            except Exception:
                continue
            tot = float(ycalc.sum())
            others = [int(ix["ids"][rows[order[j]]]) for j in fam[g] if j != k][:12]
            eq = [x[0] for x in con.execute("SELECT ident FROM entries WHERE id IN (%s)" % ",".join("?" * len(others)), others)] if others else []
            cands.append({"id": row["id"], "row": row, "ident": row["ident"], "name": row["name"] or "", "formula": row["formula"] or "",
                          "label": row["ident"] if row["source"] == "CSD" else "%s %s" % (row["source"], os.path.basename(row["ident"])),
                          "sg": row["sg"] or "", "temp_K": row["temp_K"], "s1": float(s1[order[k]]), "s2": float(sc2[k]), "s3": float(expl),
                          "cover": float(np.minimum(ycalc, obs.res).sum() / tot) if tot > 0 else 0.0, "cell": cell3, "ycalc": ycalc,
                          "sticks": (stt, sti), "sabs": sabs, "mass": mass, "cif": txt, "fam": g, "equiv": eq, "info": pi.format_info(row),
                          "shift": 100 * (np.exp(sh[order[k]] * pc.LOG_STEP) - 1)})
        cands.sort(key=lambda c: -c["s3"])
        # regroupement : même famille de raies, ou profils ajustés quasi identiques
        out = []
        for c in cands:
            n = float(np.sqrt(c["ycalc"] @ c["ycalc"])) or 1.0
            lead = next((L for L in out if L["fam"] == c["fam"] or float(L["ycalc"] @ c["ycalc"]) / (n * L["_n"]) >= 0.8), None)
            if lead is None:
                c["_n"] = n
                out.append(c)
            else:
                lead["equiv"] = [c["ident"]] + [e for e in lead["equiv"] if e != c["ident"]]
        best_other = lambda i: max([c["s3"] for j, c in enumerate(out) if j != i] or [0.0])
        for i, c in enumerate(out):
            c["contrast"] = c["s3"] / max(best_other(i), 0.02)
        self.searched = True
        return out

    # -- phases retenues ------------------------------------------------------------------
    def _refit(self):
        obs = self.obs
        if not self.accepted:
            obs.res = self.net.copy()
            return
        P = np.stack([a["ycalc0"] for a in self.accepted], 1)
        w = 1.0 / np.maximum(obs.yo, 1.0) ** ps.Fit3.WPOW
        s, _ = nnls(P * w[:, None], self.net * w)
        for a, k in zip(self.accepted, s):
            a["ycalc"] = k * a["ycalc0"]
            a["sticks"] = (a["sticks0"][0], k * a["sticks0"][1])
        obs.res = np.clip(self.net - P @ s, 0, None)
        tot = sum(float(a["ycalc"].sum()) for a in self.accepted) or 1.0
        for a in self.accepted:
            a["share"] = float(a["ycalc"].sum()) / tot
        # proportions massiques approchées : w ∝ échelle absolue × masse de la maille × volume (Hill & Howard)
        q = [float(k) * a.get("sabs", 0.0) * a.get("mass", 0.0) * pc.cell_volume(a["cell"]) for a, k in zip(self.accepted, s)]
        ok = all(a.get("sabs", 0) > 0 and a.get("mass", 0) > 0 for a in self.accepted) and sum(q) > 0
        for a, v in zip(self.accepted, q):
            a["wt"] = v / sum(q) if ok else None

    def accept(self, cand):
        a = dict(cand)
        a["ycalc0"] = cand["ycalc"].copy()
        a["sticks0"] = (cand["sticks"][0].copy(), cand["sticks"][1].copy())
        a["round"] = self.round
        a["color"] = COLORS[len(self.accepted) % len(COLORS)]
        self.accepted.append(a)
        self.round += 1
        self._refit()

    def remove(self, i):
        self.accepted.pop(i)
        self.round = len(self.accepted) + 1
        for k, a in enumerate(self.accepted):
            a["color"] = COLORS[k % len(COLORS)]
        self._refit()

    def unexplained(self):
        """Part du signal net non expliquée par les phases retenues (en intensité)."""
        return float(self.obs.res.sum() / max(self.net.sum(), 1e-30))

    # -- résultats pour le rapport --------------------------------------------------------
    def results(self, project=""):
        """État de la recherche sous forme de dict (écrit en JSON pour le rapport d'affinement)."""
        obs = self.obs
        r3 = lambda v: [round(float(t), 3) for t in v]
        phases = []
        for i, a in enumerate(self.accepted):
            row = a["row"]
            info = json.loads(row["info"] or "{}")
            keys = ("a", "b", "c", "al", "be", "ga")
            phases.append({
                "order": i + 1, "source": row["source"], "ident": os.path.basename(a["ident"]) if row["source"] != "CSD" else a["ident"],
                "name": a["name"], "formula": a["formula"], "sg": a["sg"], "temp_K": a["temp_K"],
                "share": round(a.get("share", 0.0), 4), "wt": round(a["wt"], 4) if a.get("wt") is not None else None,
                "score": round(a["s3"], 4), "cover": round(a["cover"], 4), "contrast": round(a.get("contrast", 0.0), 3),
                "cell_fit": [round(float(v), 5) for v in a["cell"]],
                "cell_db": [row[k] for k in keys] if row["a"] else None,
                "authors": info.get("authors") or "", "journal": info.get("journal") or "", "volume": info.get("volume") or "",
                "page": info.get("page") or "", "year": int(info["year"]) if info.get("year") else None, "doi": info.get("doi") or "",
                "ccdc_number": info.get("ccdc_number"), "z": info.get("z"), "equiv": list(a.get("equiv") or [])[:12],
                "color": a["color"], "ycalc": [round(float(v), 1) for v in a["ycalc"]],
                "sticks": [r3(a["sticks"][0]), [round(float(v), 1) for v in a["sticks"][1]]]})
        return {"format": RESULTS_FORMAT, "kind": "pattern", "date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
                "project": os.path.basename(project or ""), "histogram": self.pat.get("name", ""),
                "databases": self.db_stats, "sources_searched": list(self.sources) if self.sources else sorted(self.db_stats),
                "elements": list(self.elements or []), "strict": bool(self.strict),
                "unexplained": round(self.unexplained(), 4), "phases": phases,
                "x": [round(float(v), 4) for v in obs.x], "net": [round(float(v), 1) for v in self.net]}

    def save_results(self, path, project=""):
        """Écrit <projet>_phases.json (écriture atomique). Renvoie le chemin."""
        tmp = path + ".tmp"
        with open(tmp, "w", encoding="utf8") as f:
            json.dump(self.results(project), f, ensure_ascii=False)
        os.replace(tmp, path)
        return path

    # -- export ---------------------------------------------------------------------------
    def export(self, outdir):
        os.makedirs(outdir, exist_ok=True)
        files = []
        with open(os.path.join(outdir, "phases_retenues.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["ordre", "source", "identifiant", "nom", "formule", "groupe", "T_K", "part_intensite_calculee", "score_valide", "couverture",
                        "a", "b", "c", "alpha", "beta", "gamma", "fichier_cif"])
            for i, a in enumerate(self.accepted):
                name = "phase%d_%s.cif" % (i + 1, "".join(ch if ch.isalnum() else "_" for ch in os.path.basename(a["ident"])))
                txt, ndup = pc.clean_cif(a["cif"], a["cell"], title=os.path.basename(a["ident"]),
                                         comment="%s [%s] %s\nMaille ajustée sur le diagramme par la recherche de phases ; atomes en double retirés : %d\n"
                                                 "Uiso fixés à 0,038 (0,05 pour H) quand la source n'en donne pas" % (a["ident"], a["row"]["source"], a["name"], 0))
                txt = txt.replace("atomes en double retirés : 0", "atomes en double retirés : %d" % ndup)
                with open(os.path.join(outdir, name), "w") as fh:
                    fh.write(txt)
                files.append(os.path.join(outdir, name))
                w.writerow([i + 1, a["row"]["source"], a["ident"], a["name"], a["formula"], a["sg"], a["temp_K"], round(a.get("share", 0), 4),
                            round(a["s3"], 4), round(a["cover"], 4)] + [round(v, 4) for v in a["cell"]] + [name])
        return files
