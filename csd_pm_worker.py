#!/usr/bin/env python
"""Worker lancé dans le Python du CCDC (API Python de la CSD). Ne pas lancer à la main.

  list                 -> '#version X' puis un refcode par ligne
  build --todo FICHIER -> une ligne JSON par refcode (raies + fiche)
  cif REFCODE [...]    -> une ligne JSON {"ident":..., "cif":...} par refcode

NON TESTÉ hors du Mac : écrit d'après la documentation de l'API ; chaque attribut est lu
avec précaution pour qu'une différence de version ne fasse pas échouer l'ensemble.
"""
import base64, json, os, re, sys, time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"          # un seul fil par worker : le parallélisme vient du nombre de workers

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def _get(obj, name, default=None):
    try:
        v = getattr(obj, name)
        return default if v is None else v
    except Exception:
        return default


def _reader():
    from ccdc import io
    return io.EntryReader("CSD")


def cmd_list():
    from ccdc import io
    try:
        ver = io.csd_version()
    except Exception:
        ver = ""
    sys.stdout.write("#version %s\n" % ver)
    rd = _reader()
    try:
        n = len(rd)
        for i in range(n):
            sys.stdout.write(rd.identifier(i) + "\n")
    except Exception:
        for e in rd:
            sys.stdout.write(e.identifier + "\n")


def _entry(rd, pos, ident):
    """Accès par position (rapide) avec vérification du refcode ; à défaut, accès par refcode."""
    if pos is not None:
        try:
            e = rd[pos]
            if e.identifier == ident:
                return e
        except Exception:
            pass
    return rd.entry(ident)


def _structure(entry):
    import pm_core as pc
    cr = entry.crystal
    cell = list(cr.cell_lengths) + list(cr.cell_angles)
    ops = [pc.parse_symop(s) for s in cr.symmetry_operators]
    mol = None
    for attr in ("disordered_molecule", "asymmetric_unit_molecule", "molecule"):
        mol = _get(cr, attr)
        if mol is not None:
            break
    atoms = []
    for a in (mol.atoms if mol is not None else []):
        fc = _get(a, "fractional_coordinates")
        if fc is None:
            continue
        el = pc.element_of(_get(a, "atomic_symbol", "") or _get(a, "label", ""))
        if el is None:
            continue
        occ = _get(a, "occupancy", 1.0)
        atoms.append((el, float(fc[0]), float(fc[1]), float(fc[2]), float(occ if occ is not None else 1.0), None))
    return {"cell": cell, "ops": ops, "atoms": atoms, "meta": {}}


def _info(entry):
    cr = _get(entry, "crystal")
    pub = _get(entry, "publication")
    t = _get(entry, "temperature", "")
    m = re.search(r"(-?\d+(?:\.\d+)?)\s*(K|deg\.?\s*C)?", str(t) if t else "")
    temp = None
    if m:
        temp = float(m.group(1)) + (273.15 if (m.group(2) or "K").lower().startswith("deg") else 0.0)
    elif _get(entry, "has_3d_structure", False):
        temp = None      # non précisée dans la CSD = ambiante par convention ; laissé vide
    authors = _get(pub, "authors", "") if pub is not None else ""
    if isinstance(authors, (list, tuple)):
        authors = ", ".join(authors)
    info = {
        "name": _get(entry, "chemical_name", "") or "",
        "synonyms": list(_get(entry, "synonyms", []) or [])[:6],
        "formula": _get(entry, "formula", "") or "",
        "sg": _get(cr, "spacegroup_symbol", "") if cr is not None else "",
        "temp_K": temp,
        "temp_txt": str(t) if t else "",
        "rfac": _get(entry, "r_factor"),
        "z": _get(cr, "z_value") if cr is not None else None,
        "zprime": _get(cr, "z_prime") if cr is not None else None,
        "density": _get(cr, "calculated_density") if cr is not None else None,
        "polymorph": _get(entry, "polymorph", ""),
        "color": _get(entry, "color", ""),
        "habit": _get(entry, "habit", ""),
        "melting_point": _get(entry, "melting_point", ""),
        "pressure": _get(entry, "pressure", ""),
        "disorder": _get(entry, "disorder_details", "") or ("oui" if _get(entry, "is_disordered", False) else ""),
        "remarks": _get(entry, "remarks", ""),
        "authors": authors,
        "journal": str(_get(_get(pub, "journal"), "name", _get(pub, "journal", "")) or "") if pub is not None else "",
        "volume": _get(pub, "volume", "") if pub is not None else "",
        "page": _get(pub, "first_page", "") if pub is not None else "",
        "year": _get(pub, "year") if pub is not None else None,
        "doi": _get(pub, "doi", "") if pub is not None else "",
        "ccdc_number": _get(entry, "ccdc_number"),
    }
    for k, v in list(info.items()):
        if isinstance(v, float):
            info[k] = round(v, 4)
        elif v is not None and not isinstance(v, (str, int, float, list)):
            info[k] = str(v)
    return info


def cmd_build(todo):
    import pm_index as pi
    rd = _reader()
    ids = []
    with open(todo) as f:
        for x in f:
            x = x.split()
            if len(x) == 2:
                ids.append((int(x[0]), x[1]))
            elif len(x) == 1:
                ids.append((None, x[0]))
    for pos, ident in ids:
        try:
            e = _entry(rd, pos, ident)
            info = _info(e)
            if not _get(e, "has_3d_structure", False):
                rec = {"ident": ident, "status": "nocoords", "name": info["name"], "formula": info["formula"], "sg": info["sg"],
                       "temp_K": info["temp_K"], "info": info, "nlines": 0}
                try:
                    cr = e.crystal
                    cell = list(cr.cell_lengths) + list(cr.cell_angles)
                    rec.update(dict(zip(("a", "b", "c", "al", "be", "ga"), cell)))
                    rec["vol"] = cr.cell_volume
                except Exception:
                    pass
            else:
                st = _structure(e)
                rec = pi.record_from_structure("CSD", ident, st, info=info)
                rec["info"] = json.loads(rec["info"])
                if rec.get("lines") is not None:
                    rec["lines"] = base64.b64encode(rec["lines"]).decode()
                rec.pop("source", None)
        except Exception as ex:
            rec = {"ident": ident, "status": "error: %s: %s" % (type(ex).__name__, ex), "info": {}}
        sys.stdout.write(json.dumps(rec) + "\n")
        sys.stdout.flush()


def cmd_cif(idents):
    rd = _reader()
    for ident in idents:
        try:
            out = {"ident": ident, "cif": rd.entry(ident).crystal.to_string("cif")}
        except Exception as ex:
            out = {"ident": ident, "cif": "", "error": "%s: %s" % (type(ex).__name__, ex)}
        sys.stdout.write(json.dumps(out) + "\n")


def cmd_bench(n):
    """Chronomètre chaque étape sur n entrées réparties dans toute la CSD."""
    import platform
    import numpy as np
    import pm_core as pc
    import pm_index as pi
    T = {}

    def tick(k, t):
        T[k] = T.get(k, 0.0) + time.time() - t

    t = time.time(); rd = _reader(); tick("ouverture de la CSD (une fois)", t)
    tot = len(rd)
    t = time.time(); ids = [rd.identifier(i) for i in range(0, tot, max(1, tot // 2000))][:2000]; tick("liste de 2000 refcodes (une fois)", t)
    print("machine %s, Python %s, numpy %s, %d entrées dans la CSD" % (platform.machine(), platform.python_version(), np.__version__, tot))
    pos = list(range(0, tot, max(1, tot // n)))[:n]
    slow, nat, ok = [], [], 0
    for i in pos:
        t = time.time(); e = rd[i]; ident = e.identifier; tick("accès par position", t)
        t = time.time(); rd.entry(ident); tick("accès par refcode", t)
        t = time.time(); info = _info(e); tick("fiche", t)
        if not _get(e, "has_3d_structure", False):
            continue
        t0 = time.time()
        try:
            t = time.time(); cr = e.crystal; cell = list(cr.cell_lengths) + list(cr.cell_angles); ops = cr.symmetry_operators; tick("maille et symétrie", t)
            t = time.time(); mol = _get(cr, "disordered_molecule"); tick("molécule (disordered_molecule)", t)
            t = time.time(); st = _structure(e); tick("structure complète (dont boucle sur les atomes)", t)
            t = time.time(); pi.record_from_structure("CSD", ident, st, info=info); tick("calcul des raies", t)
            ok += 1
            nat.append(len(st["atoms"]) * len(st["ops"]))
        except Exception as ex:
            print("  erreur %s : %s" % (ident, ex))
        slow.append((time.time() - t0, ident))
    m = max(1, len(pos))
    print("%d entrées lues, %d avec structure ; atomes par maille : médiane %d, max %d" % (len(pos), ok, np.median(nat) if nat else 0, max(nat) if nat else 0))
    for k, v in T.items():
        print("  %-48s %8.1f ms%s" % (k, 1000 * v / (1 if "une fois" in k else m), "" if "une fois" in k else " par entrée"))
    slow.sort(reverse=True)
    print("  les plus lentes : " + ", ".join("%s %.2f s" % (i, s) for s, i in slow[:5]))


if __name__ == "__main__":
    if sys.argv[1] == "bench":
        cmd_bench(int(sys.argv[2]) if len(sys.argv) > 2 else 300)
    elif sys.argv[1] == "list":
        cmd_list()
    elif sys.argv[1] == "build":
        cmd_build(sys.argv[sys.argv.index("--todo") + 1])
    elif sys.argv[1] == "cif":
        cmd_cif(sys.argv[2:])
