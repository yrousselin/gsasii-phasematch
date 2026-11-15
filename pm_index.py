#!/usr/bin/env python3
"""pm_index : index de raies pour la recherche de phases à l'aveugle (SQLite).

  python pm_index.py build  --db index.sqlite --csd [--csd-python CHEMIN] [--workers N] [--limit N]
  python pm_index.py build  --db index.sqlite --cif-dir DOSSIER [--source NOM]
  python pm_index.py cod-sync [--cod-dir DOSSIER] [--archive]  (copie locale de la COD, puis mises à jour)
  python pm_index.py build  --db index.sqlite --cod [--cod-dir DOSSIER]
  python pm_index.py stats  --db index.sqlite
  python pm_index.py info   --db index.sqlite IDENTIFIANT
  python pm_index.py view   --db index.sqlite IDENTIFIANT      (ouvre la structure dans Mercury)
  python pm_index.py bench  [--limit 300]                      (chronomètre la lecture de la CSD, sans rien écrire)

« build » est incrémental : seules les entrées absentes de l'index sont calculées. On peut
l'interrompre (Ctrl-C) et le relancer ; relancé après une mise à jour de la CSD, il n'ajoute
que les nouveaux refcodes. L'index ne contient pas les structures, seulement les raies les
plus fortes et la fiche descriptive ; la structure est relue dans sa source à la demande.
"""
import argparse, base64, glob, json, os, queue, re, shutil, sqlite3, subprocess, sys, tempfile, threading, time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")      # calculs en parallèle par processus, pas par fils

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import pm_core as pc

SCHEMA = """
CREATE TABLE IF NOT EXISTS entries (
  id INTEGER PRIMARY KEY, source TEXT NOT NULL, ident TEXT NOT NULL, status TEXT NOT NULL,
  stamp TEXT, name TEXT, formula TEXT, elements TEXT, sg TEXT,
  a REAL, b REAL, c REAL, al REAL, be REAL, ga REAL, vol REAL, temp_K REAL,
  nlines INTEGER, lines BLOB, info TEXT, added TEXT,
  UNIQUE(source, ident));
CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT);
"""
COLS = ("source", "ident", "status", "stamp", "name", "formula", "elements", "sg", "a", "b", "c", "al", "be", "ga",
        "vol", "temp_K", "nlines", "lines", "info", "added")


def open_db(path):
    con = sqlite3.connect(path)
    con.executescript(SCHEMA)
    con.execute("PRAGMA journal_mode=WAL")
    return con


def record_from_structure(source, ident, st, info=None, stamp=""):
    """Calcule l'enregistrement d'index d'une structure (dict prêt à insérer)."""
    meta = dict(st.get("meta") or {})
    meta.update(info or {})
    rec = {"source": source, "ident": ident, "stamp": stamp, "status": "ok", "name": meta.get("name", ""),
           "formula": meta.get("formula", ""), "sg": meta.get("sg", ""), "temp_K": meta.get("temp_K"),
           "elements": "", "vol": None, "nlines": 0, "lines": None}
    cell = st.get("cell")
    if cell and all(v is not None for v in cell):
        rec.update(dict(zip(("a", "b", "c", "al", "be", "ga"), cell)))
        try:
            rec["vol"] = pc.cell_volume(cell)
        except Exception:
            pass
    try:
        if not st.get("atoms"):
            raise ValueError("nocoords")
        rec["elements"] = " " + " ".join(pc.elements_of(st)) + " "
        if not rec["formula"]:
            rec["formula"] = pc.formula_from_atoms(st)
        b, inten, hkl, dmin_eff = pc.stored_lines(st)
        rec["lines"] = pc.pack_lines(b, inten, hkl)
        rec["nlines"] = int(len(b))
        if dmin_eff > pc.D_MIN_STORE + 1e-6:
            meta["dmin_calc"] = round(dmin_eff, 3)
        if meta.get("nosym"):
            raise ValueError("symétrie indéterminée (ni opérateurs ni symbole exploitable)")
        if meta.get("sym_from"):
            meta["warning"] = "; ".join(x for x in (meta.get("warning"), "symétrie déduite du %s" % meta.pop("sym_from")) if x)
    except ValueError as e:
        rec["status"] = "nocoords" if str(e) in ("nocoords", "pas de coordonnées atomiques") else "error: %s" % e
        if "symétrie" in str(e) and not pc.have_gemmi():
            rec["stamp"] = ""           # sera retentée une fois le module gemmi installé
    except Exception as e:  # entrée inexploitable : notée pour ne pas être retentée à chaque mise à jour
        rec["status"] = "error: %s: %s" % (type(e).__name__, e)
    meta.pop("nosym", None)
    rec["info"] = json.dumps(meta, ensure_ascii=False)
    return rec


def insert(con, recs):
    now = time.strftime("%Y-%m-%d")
    rows = []
    for r in recs:
        r = dict(r)
        r["added"] = now
        if isinstance(r.get("lines"), str):
            r["lines"] = base64.b64decode(r["lines"])
        rows.append(tuple(r.get(c) for c in COLS))
    con.executemany("INSERT OR REPLACE INTO entries (%s) VALUES (%s)" % (",".join(COLS), ",".join("?" * len(COLS))), rows)
    con.commit()


def known(con, source):
    return {i: s for i, s in con.execute("SELECT ident, stamp FROM entries WHERE source=?", (source,))}


# ---- source : dossier de CIF -------------------------------------------------------------
def _cif_job(args):
    source, root, rel, stamp = args
    try:
        with open(os.path.join(root, rel), "r", errors="replace") as f:
            st = pc.structure_from_cif(f.read())
    except Exception as e:
        return {"source": source, "ident": rel, "stamp": stamp, "status": "error: %s" % e, "info": "{}"}
    rec = record_from_structure(source, rel, st, stamp=stamp)
    if not rec["name"]:
        rec["name"] = os.path.splitext(os.path.basename(rel))[0]
    return rec


def build_cif_dir(con, root, source, workers, limit=None):
    root = os.path.abspath(root)
    files = sorted(glob.glob(os.path.join(root, "**", "*.cif"), recursive=True))
    have = known(con, source)
    todo = []
    for p in files:
        rel = os.path.relpath(p, root)
        s = os.stat(p)
        stamp = "%d:%d" % (s.st_size, int(s.st_mtime))
        if have.get(rel) != stamp:
            todo.append((source, root, rel, stamp))
    if limit:
        todo = todo[:limit]
    print("[%s] %d fichiers CIF, %d déjà dans l'index, %d à traiter" % (source, len(files), len(files) - len(todo), len(todo)))
    con.execute("INSERT OR REPLACE INTO meta VALUES (?,?)", ("root:" + source, root))
    con.commit()
    if not todo:
        return 0
    t0 = time.time()
    done = 0
    batch = []
    if workers > 1 and len(todo) > 50:
        import multiprocessing as mp
        pool = mp.Pool(workers)
        it = pool.imap_unordered(_cif_job, todo, chunksize=20)
    else:
        pool, it = None, map(_cif_job, todo)
    try:
        for rec in it:
            batch.append(rec)
            done += 1
            if len(batch) >= 500:
                insert(con, batch)
                batch = []
                print("   %d/%d  (%.0f/s)" % (done, len(todo), done / (time.time() - t0)), flush=True)
    finally:
        insert(con, batch)
        if pool:
            pool.terminate()
    print("[%s] %d entrées ajoutées en %.1f s" % (source, done, time.time() - t0))
    return done


# ---- source : COD (copie locale de la collection de CIF) ---------------------------------
COD_RSYNC = "rsync://www.crystallography.net/cif/"
COD_ARCHIVE = "http://www.crystallography.net/archives/cod-cifs-mysql.tgz"
COD_DEFAULT_DIR = os.path.join(os.path.expanduser("~"), "COD", "cif")
COD_MARK = ".pm_cod_mirror"
_COD_ID = re.compile(r"^(\d{7})\.cif$")


def _clean_txt(v):
    if v is None:
        return ""
    v = " ".join(str(v).split())
    return "" if v in ("?", ".") else v


def cod_dir(con, hint=None):
    """Dossier de la copie locale : --cod-dir, sinon celui déjà enregistré dans l'index, sinon ~/COD/cif."""
    if hint:
        return os.path.abspath(os.path.expanduser(hint))
    r = con.execute("SELECT value FROM meta WHERE key='root:COD'").fetchone() if con is not None else None
    return r[0] if r else COD_DEFAULT_DIR


def cod_info(st):
    """Fiche descriptive et motif d'exclusion éventuel d'une entrée COD. -> (info, nom, exclusion|None)"""
    tg = st.get("tags") or {}
    g = lambda k: _clean_txt(tg.get(k))
    names = [x for x in (g("_chemical_name_mineral"), g("_chemical_name_common"), g("_chemical_name_systematic")) if x]
    info = {"synonyms": names[1:], "authors": "; ".join(_clean_txt(a) for a in (st.get("authors") or [])[:12]),
            "journal": g("_journal_name_full"), "volume": g("_journal_volume"), "page": g("_journal_page_first"),
            "year": pc._num(tg.get("_journal_year")), "doi": g("_journal_paper_doi"),
            "z": (lambda z: int(z) if z is not None and float(z).is_integer() else z)(pc._num(tg.get("_cell_formula_units_z"))), "density": pc._num(tg.get("_exptl_crystal_density_diffrn")),
            "color": g("_exptl_crystal_colour"), "remarks": g("_publ_section_title")[:300]}
    r = pc._num(tg.get("_refine_ls_r_factor_gt") or tg.get("_refine_ls_r_factor_all"))
    if r is not None:
        info["rfac"] = round(100 * r, 2) if r < 1 else r
    p = pc._num(tg.get("_diffrn_ambient_pressure"))
    if p is not None and abs(p - 101.3) > 5:
        info["pressure"] = "%g kPa" % p
    method = g("_cod_struct_determination_method").lower()
    flag = g("_cod_error_flag").lower()
    dup = g("_cod_duplicate_entry")
    warn = []
    if flag == "warnings":
        warn.append("entrée signalée « warnings » par la COD")
    if g("_cod_related_optimal_struct"):
        warn.append("la COD indique une meilleure détermination : %s" % g("_cod_related_optimal_struct"))
    if "powder" in method:
        warn.append("structure déterminée sur poudre")
    if warn:
        info["warning"] = " ; ".join(warn)
    info = {k: v for k, v in info.items() if v not in (None, "", [])}
    excl = None
    if flag in ("errors", "retracted"):
        excl = flag
    elif "theor" in method or "predict" in method:
        excl = "theoretical"
    elif dup:
        excl = "duplicate of %s" % dup
    return info, (names[0] if names else ""), excl


def _cod_job(args):
    root, rel, stamp = args
    ident = os.path.splitext(os.path.basename(rel))[0]
    base = {"source": "COD", "ident": ident, "stamp": stamp}
    try:
        with open(os.path.join(root, rel), "r", errors="replace") as f:
            st = pc.structure_from_cif(f.read())
    except Exception as e:
        return dict(base, status="error: %s" % e, info=json.dumps({"relpath": rel}))
    info, name, excl = cod_info(st)
    info["relpath"] = rel
    if excl:        # écartée : la fiche est gardée, les raies ne sont pas calculées
        meta = dict(st.get("meta") or {})
        rec = dict(base, status="excluded: " + excl, name=name or meta.get("name", ""), formula=meta.get("formula", ""),
                   sg=meta.get("sg", ""), temp_K=meta.get("temp_K"), elements="", nlines=0, lines=None)
        cell = st.get("cell")
        if cell and all(v is not None for v in cell):
            rec.update(dict(zip(("a", "b", "c", "al", "be", "ga"), cell)))
        rec["info"] = json.dumps(info, ensure_ascii=False)
        return rec
    rec = record_from_structure("COD", ident, st, info=info, stamp=stamp)
    if name:
        rec["name"] = name
    rec["name"] = _clean_txt(rec.get("name"))
    return rec


def scan_cod(root):
    """Parcourt la copie locale. -> {identifiant: (chemin relatif, taille:date)}"""
    out = {}
    stack = [root]
    while stack:
        d = stack.pop()
        try:
            it = list(os.scandir(d))
        except OSError:
            continue
        for e in it:
            if e.is_dir(follow_symlinks=False):
                if not e.name.startswith("."):
                    stack.append(e.path)
            else:
                m = _COD_ID.match(e.name)
                if m:
                    s = e.stat()
                    out[m.group(1)] = (os.path.relpath(e.path, root), "%d:%d" % (s.st_size, int(s.st_mtime)))
    return out


def build_cod(con, root, workers, limit=None, rebuild=False):
    root = cod_dir(con, root)
    if not os.path.isdir(root):
        sys.exit("copie locale de la COD introuvable (%s) : lancer d'abord « pm_index.py cod-sync »" % root)
    if not pc.have_gemmi():
        print("[COD] attention : module gemmi absent. Les fichiers sans opérateurs de symétrie (quelques pour cent)\n"
              "      seront laissés de côté ; installer avec :  %s -m pip install gemmi  puis relancer." % sys.executable, flush=True)
    print("[COD] parcours de %s…" % root, flush=True)
    files = scan_cod(root)
    if not files:
        sys.exit("aucun fichier CIF de la COD (7 chiffres .cif) dans %s" % root)
    have = {} if rebuild else known(con, "COD")
    todo = [(root, rel, stamp) for ident, (rel, stamp) in sorted(files.items()) if have.get(ident) != stamp]
    nmiss = len(todo)
    if limit:
        todo = todo[::max(1, len(todo) // limit)][:limit]
    gone = [i for i in known(con, "COD") if i not in files]
    print("[COD] %d fichiers, %d déjà dans l'index, %d à traiter%s" % (
        len(files), len(files) - nmiss, len(todo), ", %d retirés de la COD" % len(gone) if gone else ""))
    con.execute("INSERT OR REPLACE INTO meta VALUES (?,?)", ("root:COD", root))
    if gone and len(gone) < 0.2 * len(files):     # entrées supprimées en amont (pas un dossier incomplet)
        con.executemany("DELETE FROM entries WHERE source='COD' AND ident=?", [(i,) for i in gone])
    con.commit()
    if not todo:
        return 0
    t0 = time.time()
    done, batch = 0, []
    if workers > 1 and len(todo) > 50:
        import multiprocessing as mp
        pool = mp.Pool(workers)
        it = pool.imap_unordered(_cod_job, todo, chunksize=40)
    else:
        pool, it = None, map(_cod_job, todo)
    try:
        for rec in it:
            batch.append(rec)
            done += 1
            if len(batch) >= 1000:
                insert(con, batch)
                batch = []
                rate = done / (time.time() - t0)
                print("   %d/%d  (%.0f/s, reste ~%.0f min)" % (done, len(todo), rate, (len(todo) - done) / rate / 60), flush=True)
    except KeyboardInterrupt:
        print("\ninterrompu : ce qui est fait est conservé, relancer la même commande pour continuer")
    finally:
        insert(con, batch)
        if pool:
            pool.terminate()
    dt = time.time() - t0
    print("[COD] %d entrées traitées en %.1f min (%.0f par seconde)" % (done, dt / 60, done / max(dt, 1e-9)))
    return done


def _mirror_ok(dest):
    """Le dossier peut-il recevoir la copie ? (rsync --delete efface ce qui n'est pas dans la COD)"""
    if not os.path.isdir(dest):
        return True
    names = [n for n in os.listdir(dest) if not n.startswith(".")]
    if os.path.exists(os.path.join(dest, COD_MARK)) or not names:
        return True
    return all(re.match(r"^\d$", n) for n in names)        # déjà une arborescence COD (1 … 9)


def cod_sync(dest, archive=False, url=None):
    """Crée ou met à jour la copie locale des CIF de la COD (rsync ; --archive : archive HTTP)."""
    dest = os.path.abspath(os.path.expanduser(dest))
    if not _mirror_ok(dest):
        sys.exit("%s contient autre chose que la COD : choisir un dossier vide ou dédié (--cod-dir)" % dest)
    os.makedirs(dest, exist_ok=True)
    with open(os.path.join(dest, COD_MARK), "w") as f:
        f.write("copie locale de la COD, entretenue par pm_index.py cod-sync\n")
    free = shutil.disk_usage(dest).free / 1e9
    print("[COD] dossier : %s   (espace libre : %.0f Go)" % (dest, free), flush=True)
    if archive:
        return _cod_sync_archive(dest, url or COD_ARCHIVE)
    rs = shutil.which("rsync")
    if not rs:
        sys.exit("rsync introuvable : relancer avec --archive")
    cmd = [rs, "-a", "-v", "--delete", "--exclude", COD_MARK, COD_RSYNC, dest + os.sep]
    print("[COD] " + " ".join(cmd), flush=True)
    print("[COD] première copie : longue (toute la collection) ; les suivantes ne prennent que les nouveautés.", flush=True)
    t0, n, last = time.time(), 0, time.time()
    try:
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, errors="replace")
        tail = []
        for ln in p.stdout:
            ln = ln.rstrip()
            if ln.endswith(".cif"):
                n += 1
                if time.time() - last > 10:
                    last = time.time()
                    print("   %d fichiers reçus (%.0f min) … %s" % (n, (time.time() - t0) / 60, ln[-24:]), flush=True)
            elif ln:
                tail = (tail + [ln])[-8:]
        rc = p.wait()
    except KeyboardInterrupt:
        p.terminate()
        print("\ninterrompu : relancer la même commande pour reprendre là où la copie s'est arrêtée")
        return 1
    if rc != 0:
        print("\n".join(tail))
        print("[COD] rsync a échoué (code %d). Si le port 873 est bloqué par le réseau, essayer :\n"
              "      pm_index.py cod-sync --archive   (archive HTTP, mise à jour par la COD tous les 3-4 mois)" % rc)
        return rc
    print("[COD] copie à jour : %d fichiers reçus ou modifiés en %.1f min" % (n, (time.time() - t0) / 60))
    return 0


def _download_resume(url, path):
    """Télécharge url dans path en reprenant là où une coupure s'est produite (requêtes HTTP Range)."""
    import urllib.request, urllib.error
    total, fails, last, t0 = None, 0, time.time(), time.time()
    while True:
        have = os.path.getsize(path) if os.path.exists(path) else 0
        if total is not None and have >= total:
            return True
        req = urllib.request.Request(url, headers={"User-Agent": "pm_index (GSAS-II user script)"})
        if have:
            req.add_header("Range", "bytes=%d-" % have)
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                if have and resp.status != 206:           # le serveur ignore la reprise : on repart de zéro
                    have = 0
                rng = resp.headers.get("Content-Range")   # bytes a-b/total
                if rng and "/" in rng and rng.split("/")[-1].isdigit():
                    total = int(rng.split("/")[-1])
                elif resp.headers.get("Content-Length"):
                    total = have + int(resp.headers["Content-Length"])
                with open(path, "ab" if have else "wb") as f:
                    while True:
                        buf = resp.read(1 << 20)
                        if not buf:
                            break
                        f.write(buf)
                        fails = 0
                        if time.time() - last > 10:
                            last = time.time()
                            done = f.tell()
                            print("   %.2f Go%s  (%.0f min)" % (done / 1e9, " sur %.2f (%.0f %%)" % (total / 1e9, 100.0 * done / total)
                                                              if total else "", (time.time() - t0) / 60), flush=True)
        except urllib.error.HTTPError as e:
            if e.code == 416:                             # plage hors fichier : tout est déjà là
                return True
            fails += 1
        except KeyboardInterrupt:
            print("\ninterrompu : relancer la même commande pour reprendre le téléchargement")
            return False
        except Exception as e:
            fails += 1
            print("   coupure (%s) : reprise…" % str(e)[:60], flush=True)
        if total is None and fails == 0:                  # taille inconnue et fin de flux normale
            return True
        if fails >= 30:
            print("[COD] trop de coupures successives : relancer la même commande plus tard (le téléchargement reprendra)")
            return False
        if fails:
            time.sleep(min(60, 3 * fails))


def _cod_sync_archive(dest, url, keep=False):
    """Télécharge l'archive de la COD (avec reprise) à côté du dossier, puis en extrait les CIF."""
    import tarfile
    if url == COD_ARCHIVE:
        try:
            import lzma  # noqa: F401
            url = COD_ARCHIVE[:-4] + ".txz"               # même contenu, nettement plus compact
        except ImportError:
            pass
    arch = os.path.join(os.path.dirname(dest.rstrip(os.sep)) or dest, os.path.basename(url) + ".part")
    print("[COD] téléchargement de %s\n      vers %s (reprise automatique après coupure)" % (url, arch), flush=True)
    if not _download_resume(url, arch):
        return 1
    print("[COD] téléchargement terminé (%.2f Go) ; extraction des CIF…" % (os.path.getsize(arch) / 1e9), flush=True)
    t0, n, nskip, last = time.time(), 0, 0, time.time()
    part = re.compile(r"^[0-9A-Za-z_.-]+$")
    try:
        with tarfile.open(arch, mode="r|*") as tar:
            for m in tar:
                if not m.isfile() or not _COD_ID.match(os.path.basename(m.name)):
                    continue
                parts = m.name.replace("\\", "/").split("/")
                if "cif" in parts[:-1]:
                    parts = parts[len(parts) - parts[::-1].index("cif"):]
                else:
                    parts = parts[-4:]
                if any(not part.match(x) or x in (".", "..") for x in parts):
                    continue
                path = os.path.join(dest, *parts)
                try:
                    st = os.stat(path)
                    if st.st_size == m.size and int(st.st_mtime) == int(m.mtime):
                        nskip += 1
                        continue
                except OSError:
                    pass
                os.makedirs(os.path.dirname(path), exist_ok=True)
                src = tar.extractfile(m)
                with open(path, "wb") as f:
                    shutil.copyfileobj(src, f)
                if m.mtime:
                    os.utime(path, (m.mtime, m.mtime))
                n += 1
                if time.time() - last > 10:
                    last = time.time()
                    print("   %d fichiers extraits (%.0f min)" % (n + nskip, (time.time() - t0) / 60), flush=True)
    except (tarfile.ReadError, EOFError, OSError) as e:
        print("[COD] archive illisible ou incomplète (%s) : supprimer %s et relancer" % (e, arch))
        return 1
    print("[COD] %d fichiers CIF extraits (%d déjà présents) en %.1f min" % (n, nskip, (time.time() - t0) / 60))
    if (n + nskip) and not keep:
        os.remove(arch)
    return 0 if (n + nskip) else 1


# ---- source : CSD locale (API Python du CCDC, en sous-processus) -------------------------
def find_csd_python(hint=None):
    cands = [hint] if hint else []
    cands += [os.environ.get("CSD_PYTHON", "")]
    home = os.path.expanduser("~")
    for pat in ("CCDC/ccdc-software/csd-python-api/miniconda/bin/python", "CCDC/*/csd-python-api/miniconda/bin/python",
                "CCDC/Python_API_*/miniconda/bin/python", "CCDC/**/miniconda/bin/python"):
        cands += sorted(glob.glob(os.path.join(home, pat), recursive=True))
    cands += sorted(glob.glob("/Applications/CCDC/**/miniconda/bin/python", recursive=True))
    for c in cands:
        if c and os.path.isfile(c):
            return c
    return None


def build_csd(con, csd_python, workers, limit=None, rebuild=False):
    py = find_csd_python(csd_python)
    if not py:
        sys.exit("Python du CCDC introuvable : donner --csd-python CHEMIN")
    worker = os.path.join(HERE, "csd_pm_worker.py")
    print("[CSD] lecture de la liste des refcodes (%s)…" % py, flush=True)
    out = subprocess.run([py, worker, "list"], capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit("échec du worker CSD :\n" + out.stderr[-2000:])
    ids, version = [], ""
    for ln in out.stdout.splitlines():
        if ln.startswith("#version"):
            version = ln.split(None, 1)[1] if " " in ln else ""
        elif ln.strip():
            ids.append(ln.strip())
    have = {} if rebuild else known(con, "CSD")
    todo = [(k, i) for k, i in enumerate(ids) if i not in have]
    nmiss = len(todo)
    if limit:
        # échantillon réparti dans toute la base (les premiers refcodes, anciens, ne sont pas représentatifs)
        todo = todo[::max(1, len(todo) // limit)][:limit]
    print("[CSD %s] %d refcodes, %d déjà dans l'index, %d à traiter" % (version, len(ids), len(ids) - nmiss, len(todo)))
    con.execute("INSERT OR REPLACE INTO meta VALUES (?,?)", ("csd_version", version))
    con.execute("INSERT OR REPLACE INTO meta VALUES (?,?)", ("csd_python", py))
    con.commit()
    if not todo:
        return 0
    workers = max(1, min(workers, len(todo) // 200 + 1))
    tmp = tempfile.mkdtemp(prefix="pm_csd_")
    q = queue.Queue(maxsize=5000)
    procs = []

    def pump(p):
        for ln in p.stdout:
            q.put(ln)
        q.put(None)

    for w in range(workers):
        tf = os.path.join(tmp, "todo_%d.txt" % w)
        with open(tf, "w") as f:
            f.write("\n".join("%d %s" % t for t in todo[w::workers]) + "\n")
        p = subprocess.Popen([py, worker, "build", "--todo", tf], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
        threading.Thread(target=pump, args=(p,), daemon=True).start()
        procs.append(p)
    t0 = time.time()
    done, alive, batch = 0, workers, []
    try:
        while alive:
            ln = q.get()
            if ln is None:
                alive -= 1
                continue
            try:
                rec = json.loads(ln)
            except ValueError:
                continue
            rec["source"] = "CSD"
            rec["stamp"] = version
            rec["info"] = json.dumps(rec.get("info") or {}, ensure_ascii=False)
            batch.append(rec)
            done += 1
            if len(batch) >= 1000:
                insert(con, batch)
                batch = []
                rate = done / (time.time() - t0)
                print("   %d/%d  (%.0f/s, reste ~%.0f min)" % (done, len(todo), rate, (len(todo) - done) / rate / 60), flush=True)
    except KeyboardInterrupt:
        print("\ninterrompu : ce qui est fait est conservé, relancer la même commande pour continuer")
    finally:
        insert(con, batch)
        for p in procs:
            p.terminate()
    dt = time.time() - t0
    print("[CSD] %d entrées ajoutées en %.1f min (%.0f par seconde avec %d workers)" % (done, dt / 60, done / max(dt, 1e-9), workers))
    return done


# ---- lecture de l'index pour la recherche ------------------------------------------------
def load_index(db):
    """-> dict : ids, bins (n, N_STORE) int32 (-1 = vide), inten float32, hkl int8 (n, N_STORE, 3), cell (n, 6).
    Cache .npz à côté de la base."""
    con = open_db(db)
    n, mx = con.execute("SELECT COUNT(*), COALESCE(MAX(id),0) FROM entries WHERE status='ok'").fetchone()
    tot = con.execute("SELECT COUNT(*) FROM entries").fetchone()[0]
    key = "%d:%d:%d" % (n, mx, tot)
    cache = db + ".cache.npz"
    if os.path.exists(cache):
        try:
            z = np.load(cache)
            if str(z["key"]) == key:
                return {k: z[k] for k in ("ids", "bins", "inten", "hkl", "cell")}
        except Exception:
            pass
    ids = np.zeros(n, np.int64)
    bins = np.full((n, pc.N_STORE), -1, np.int32)
    inten = np.zeros((n, pc.N_STORE), np.float32)
    hkl = np.zeros((n, pc.N_STORE, 3), np.int8)
    cell = np.zeros((n, 6), np.float64)
    for i, row in enumerate(con.execute("SELECT id, lines, a, b, c, al, be, ga FROM entries WHERE status='ok' ORDER BY id")):
        b, x, h = pc.unpack_lines(row[1])
        m = min(len(b), pc.N_STORE)
        ids[i] = row[0]
        bins[i, :m] = b[:m]
        inten[i, :m] = x[:m]
        hkl[i, :m] = h[:m]
        cell[i] = row[2:8]
    out = {"ids": ids, "bins": bins, "inten": inten, "hkl": hkl, "cell": cell}
    try:
        np.savez(cache, key=key, **out)
    except OSError:
        pass
    return out


def get_row(con, ident_or_id):
    cur = con.cursor()
    cur.row_factory = sqlite3.Row
    if isinstance(ident_or_id, (int, np.integer)):
        return cur.execute("SELECT * FROM entries WHERE id=?", (int(ident_or_id),)).fetchone()
    r = cur.execute("SELECT * FROM entries WHERE ident=? ORDER BY source", (ident_or_id,)).fetchone()
    return r or cur.execute("SELECT * FROM entries WHERE ident LIKE ? ORDER BY source", ("%" + ident_or_id + "%",)).fetchone()


def get_cifs(con, rows):
    """Texte CIF de plusieurs entrées, relu dans leur source (un seul lancement du Python du CCDC). -> {id: texte}"""
    out = {}
    csd = [r for r in rows if r["source"] == "CSD"]
    if csd:
        py = (con.execute("SELECT value FROM meta WHERE key='csd_python'").fetchone() or [find_csd_python()])[0]
        res = subprocess.run([py, os.path.join(HERE, "csd_pm_worker.py"), "cif"] + [r["ident"] for r in csd], capture_output=True, text=True)
        by_ident = {}
        for ln in res.stdout.splitlines():
            if ln.startswith("{"):
                try:
                    d = json.loads(ln)
                    by_ident[d["ident"]] = d.get("cif", "")
                except ValueError:
                    pass
        for r in csd:
            out[r["id"]] = by_ident.get(r["ident"], "")
    for r in rows:
        if r["source"] != "CSD":
            root = con.execute("SELECT value FROM meta WHERE key=?", ("root:" + r["source"],)).fetchone()[0]
            rel = r["ident"]
            if r["source"] == "COD":
                i = r["ident"]
                rel = json.loads(r["info"] or "{}").get("relpath") or os.path.join(i[0], i[1:3], i[3:5], i + ".cif")
            try:
                with open(os.path.join(root, rel), "r", errors="replace") as f:
                    out[r["id"]] = f.read()
            except OSError:
                out[r["id"]] = ""
    return out


def get_structure(con, row):
    """Relit la structure complète dans sa source (dossier de CIF ou CSD). -> (structure, texte CIF)."""
    txt = get_cifs(con, [row])[row["id"]]
    if not txt:
        raise ValueError("structure illisible dans sa source")
    return pc.structure_from_cif(txt), txt


def open_in_mercury(cif_text, ident):
    """Écrit le CIF dans un fichier temporaire et l'ouvre dans Mercury. -> message"""
    path = os.path.join(tempfile.gettempdir(), "pm_%s.cif" % "".join(c if c.isalnum() else "_" for c in ident))
    with open(path, "w") as f:
        f.write(cif_text)
    apps = glob.glob(os.path.expanduser("~/CCDC/**/[Mm]ercury.app"), recursive=True) + glob.glob("/Applications/**/[Mm]ercury.app", recursive=True)
    if sys.platform == "darwin" and apps:
        subprocess.Popen(["open", "-a", apps[0], path])
        return "Ouvert dans Mercury : " + path
    if sys.platform == "darwin":
        subprocess.Popen(["open", path])
        return "Mercury introuvable dans ~/CCDC : CIF ouvert avec l'application par défaut (%s)" % path
    return "Mercury introuvable ; CIF écrit ici : " + path


def format_info(row):
    info = json.loads(row["info"] or "{}")
    L = ["%s  [%s]" % (row["ident"], row["source"])]

    def add(lab, v):
        if v not in (None, "", [], "?"):
            L.append("  %-22s %s" % (lab, v))

    add("Nom", row["name"])
    add("Synonymes", "; ".join(info.get("synonyms") or []))
    add("Formule", row["formula"])
    add("Groupe d'espace", row["sg"])
    if row["a"]:
        add("Maille (Å, °)", "%.4f %.4f %.4f  %.3f %.3f %.3f" % tuple(row[k] for k in ("a", "b", "c", "al", "be", "ga")))
        add("Volume (Å³)", "%.2f" % row["vol"] if row["vol"] else None)
    add("Z, Z'", " / ".join(str(info[k]) for k in ("z", "zprime") if info.get(k) is not None))
    add("Densité calc.", info.get("density"))
    add("Température (K)", row["temp_K"])
    add("Pression", info.get("pressure"))
    add("Facteur R (%)", info.get("rfac"))
    add("Polymorphe", info.get("polymorph"))
    add("Couleur / faciès", " / ".join(x for x in (info.get("color"), info.get("habit")) if x))
    add("Point de fusion", info.get("melting_point"))
    add("Désordre", info.get("disorder"))
    add("Remarques", info.get("remarks"))
    add("Auteurs", info.get("authors"))
    add("Référence", " ".join(str(x) for x in (info.get("journal"), info.get("volume"), info.get("page"), "(%s)" % int(info["year"]) if info.get("year") else None) if x))
    add("DOI", info.get("doi"))
    add("N° CCDC", info.get("ccdc_number"))
    if row["source"] == "COD":
        add("Fiche COD", "https://www.crystallography.net/cod/%s.html" % row["ident"])
    add("Avertissement", info.get("warning"))
    add("État dans l'index", "%s, %s raies, ajoutée le %s" % (row["status"], row["nlines"], row["added"]))
    return "\n".join(L)


def view_in_mercury(con, row):
    st, txt = get_structure(con, row)
    return open_in_mercury(txt, row["ident"])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=("build", "stats", "info", "view", "bench", "cod-sync"))
    ap.add_argument("ident", nargs="?")
    ap.add_argument("--db", default=os.path.expanduser("~/.GSASII/pm_index.sqlite"))
    ap.add_argument("--csd", action="store_true")
    ap.add_argument("--csd-python")
    ap.add_argument("--cif-dir")
    ap.add_argument("--source", help="nom donné au dossier de CIF (défaut : nom du dossier)")
    ap.add_argument("--cod", action="store_true", help="indexer la copie locale de la COD (voir cod-sync)")
    ap.add_argument("--cod-dir", help="dossier de la copie locale de la COD (défaut : ~/COD/cif)")
    ap.add_argument("--archive", action="store_true", help="cod-sync : archive HTTP au lieu de rsync")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--limit", type=int, help="ne traiter que N nouvelles entrées (essai)")
    ap.add_argument("--rebuild", action="store_true", help="avec --csd : recalculer toutes les entrées, y compris celles déjà dans l'index")
    a = ap.parse_args()
    if a.cmd == "bench":
        py = find_csd_python(a.csd_python) or sys.exit("Python du CCDC introuvable : donner --csd-python CHEMIN")
        subprocess.run([py, os.path.join(HERE, "csd_pm_worker.py"), "bench", str(a.limit or 300)])
        return
    con = open_db(a.db)
    if a.cmd == "cod-sync":
        dest = cod_dir(con, a.cod_dir)
        rc = cod_sync(dest, archive=a.archive)
        if rc == 0:
            con.execute("INSERT OR REPLACE INTO meta VALUES (?,?)", ("root:COD", dest))
            con.commit()
            print("Étape suivante :  pm_index.py build --cod")
        sys.exit(rc)
    if a.cmd == "build":
        if a.cod:
            build_cod(con, a.cod_dir, a.workers, a.limit, a.rebuild)
        if a.cif_dir:
            build_cif_dir(con, a.cif_dir, a.source or os.path.basename(os.path.abspath(a.cif_dir)), a.workers, a.limit)
        if a.csd:
            build_csd(con, a.csd_python, a.workers, a.limit, a.rebuild)
        if not (a.cif_dir or a.csd or a.cod):
            sys.exit("donner --csd, --cod et/ou --cif-dir DOSSIER")
    if a.cmd in ("build", "stats"):
        con.execute("PRAGMA wal_checkpoint(TRUNCATE)")
        print("Index %s (%.1f Mo)" % (a.db, os.path.getsize(a.db) / 1e6))
        for src, n, ok, noc, exc in con.execute("SELECT source, COUNT(*), SUM(status='ok'), SUM(status='nocoords'), "
                                                "SUM(status LIKE 'excluded%') FROM entries GROUP BY source"):
            print("  %-12s %8d entrées : %d utilisables, %d sans coordonnées, %d écartées, %d en erreur" % (src, n, ok, noc, exc, n - ok - noc - exc))
        for st, n in con.execute("SELECT status, COUNT(*) FROM entries WHERE source='COD' AND status LIKE 'excluded%' "
                                 "AND status NOT LIKE 'excluded: duplicate%' GROUP BY status"):
            print("       COD %-28s %d" % (st, n))
        ndup = con.execute("SELECT COUNT(*) FROM entries WHERE source='COD' AND status LIKE 'excluded: duplicate%'").fetchone()[0]
        if ndup:
            print("       COD excluded: duplicate          %d" % ndup)
    if a.cmd in ("info", "view"):
        row = get_row(con, a.ident)
        if row is None:
            sys.exit("entrée inconnue : %s" % a.ident)
        print(format_info(row))
        if a.cmd == "view":
            print(view_in_mercury(con, row))


if __name__ == "__main__":
    main()
