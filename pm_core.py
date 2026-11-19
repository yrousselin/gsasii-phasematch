"""pm_core : lecture de structures et calcul des raies de diffraction (numpy seul).

Utilisé par l'index (pm_index.py), le worker CSD (csd_pm_worker.py) et la recherche
(pm_search.py). Aucune dépendance à GSAS-II ni à l'API du CCDC.
"""
import os
import re
import sys
import numpy as np

from pm_f0 import F0

# ---- grille commune en ln(d) : un pas = 0,025 % en d ------------------------------------
D_MIN_STORE = 1.5        # Å : raies conservées dans l'index (2θ ≤ 61,8° en Cu Kα)
D_MAX_STORE = 60.0
LOG_STEP = 2.5e-4
N_STORE = 48             # raies les plus fortes conservées par entrée
MAX_HKL_STORE = 3000     # réflexions calculées au plus par entrée de l'index (dmin relevé pour les grandes mailles)
B_OVERALL = 3.0          # Å², agitation thermique globale (Uiso ≈ 0,038)
CU_KA1 = 1.540598


def shared_db_path():
    """Emplacement commun à tous les utilisateurs de l'ordinateur."""
    if os.name == "nt":
        return os.path.join(os.environ.get("ProgramData", r"C:\ProgramData"), "GSASII", "pm_index.sqlite")
    return "/Users/Shared/GSASII/pm_index.sqlite" if sys.platform == "darwin" else "/usr/local/share/GSASII/pm_index.sqlite"


def default_db():
    """Index utilisé par défaut : variable PM_INDEX_DB, sinon celui de l'utilisateur (~/.GSASII), sinon l'index commun."""
    env = os.environ.get("PM_INDEX_DB")
    if env:
        return os.path.expanduser(env)
    own = os.path.expanduser("~/.GSASII/pm_index.sqlite")
    if os.path.exists(own):
        return own
    shared = shared_db_path()
    return shared if os.path.exists(shared) else own


def d_to_bin(d):
    return np.rint(np.log(np.asarray(d, float) / D_MIN_STORE) / LOG_STEP).astype(np.int64)


def bin_to_d(b):
    return D_MIN_STORE * np.exp(np.asarray(b, float) * LOG_STEP)


N_BINS = int(d_to_bin(D_MAX_STORE)) + 1


# ---- CIF --------------------------------------------------------------------------------
_TOK = re.compile(r"""'(.*?)'(?=\s|$)|"(.*?)"(?=\s|$)|(\#.*)|(\S+)""")


def _tokens(text):
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith(";"):
            buf = [ln[1:]]
            i += 1
            while i < len(lines) and not lines[i].startswith(";"):
                buf.append(lines[i])
                i += 1
            i += 1
            yield ("v", "\n".join(buf).strip())
            continue
        for m in _TOK.finditer(ln):
            if m.group(3) is not None:
                break
            if m.group(1) is not None:
                yield ("v", m.group(1))
            elif m.group(2) is not None:
                yield ("v", m.group(2))
            else:
                t = m.group(4)
                yield ("t", t)
        i += 1


def parse_cif_blocks(text):
    """Liste de blocs : chacun = (nom, {tag: valeur}, [ {tag: [valeurs]} ... ])."""
    blocks = []
    tags, loops, name = None, None, None
    toks = list(_tokens(text))
    i, n = 0, len(toks)

    def is_tag(k):
        return toks[k][0] == "t" and toks[k][1].startswith("_")

    def is_kw(k):
        if toks[k][0] != "t":
            return False
        low = toks[k][1].lower()
        return low == "loop_" or low.startswith("data_") or low.startswith("save_") or low == "global_"

    while i < n:
        kind, t = toks[i]
        low = t.lower() if kind == "t" else ""
        if kind == "t" and low.startswith("data_"):
            tags, loops, name = {}, [], t[5:]
            blocks.append((name, tags, loops))
            i += 1
        elif tags is None:
            i += 1
        elif kind == "t" and low == "loop_":
            i += 1
            cols = []
            while i < n and is_tag(i):
                cols.append(toks[i][1].lower())
                i += 1
            vals = []
            while i < n and not (is_tag(i) or is_kw(i)):
                vals.append(toks[i][1])
                i += 1
            if cols:
                nrow = len(vals) // len(cols)
                loops.append({c: [vals[r * len(cols) + k] for r in range(nrow)] for k, c in enumerate(cols)})
        elif kind == "t" and t.startswith("_"):
            if i + 1 < n and not (is_tag(i + 1) or is_kw(i + 1)):
                tags[t.lower()] = toks[i + 1][1]
                i += 2
            else:
                i += 1
        else:
            i += 1
    return blocks


def _num(s, default=None):
    if s is None:
        return default
    s = s.strip()
    if s in ("?", ".", ""):
        return default
    m = re.match(r"^([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)", s)
    return float(m.group(1)) if m else default


def element_of(sym, allowed=None):
    """Symbole d'élément à partir d'un type ou d'une étiquette d'atome (C1, O2-, Fe3+, D, Wat1...).
    allowed : éléments de la formule brute, pour trancher les étiquettes ambiguës (HO1 = H ou Ho ?)."""
    sym = sym.strip()
    m = re.match(r"^([A-Za-z]{1,2})", sym)
    if not m:
        return None
    if re.match(r"^(wat|h2o|ow)", sym, re.I):      # molécule d'eau des fichiers de minéraux (AMCSD)
        return "O"
    s = m.group(1)
    two = s[0].upper() + s[1:].lower()
    one = s[0].upper()
    if one in ("D", "T") and not (len(two) == 2 and two in F0 and (allowed is None or two in allowed)):
        return "H"
    cands = []
    if len(two) == 2 and two in F0 and two not in ("Va",):
        cands.append(two)
    if one in F0 and one != "Q":
        cands.append(one)
    if allowed:
        for c in cands:
            if c in allowed:
                return c
    return cands[0] if cands else None


def formula_elements(formula):
    """Éléments cités dans une formule brute (D compté comme H). -> set, vide si illisible."""
    out = set()
    for e in re.findall(r"[A-Z][a-z]?", formula or ""):
        if e in ("D", "T"):
            e = "H"
        if e in F0:
            out.add(e)
    return out


def have_gemmi():
    try:
        import gemmi  # noqa: F401
        return True
    except ImportError:
        return False


_HALL_CACHE = {}
_KEY_CACHE = {}


def _op_key(op):
    R, t = op
    return (tuple(np.rint(R).astype(int).ravel()), tuple(np.rint((np.asarray(t) % 1.0) * 24).astype(int) % 24))


def _is_group(ops):
    """La liste d'opérateurs est-elle fermée pour la composition (modulo les translations de réseau) ?"""
    keys = {_op_key(o) for o in ops}
    if len(keys) != len(ops):
        return False
    for R1, t1 in ops:
        for R2, t2 in ops:
            if _op_key((R1 @ R2, R1 @ np.asarray(t2) + np.asarray(t1))) not in keys:
                return False
    return True


def _formula_count(formula):
    """Nombre d'atomes autres que H dans une formule brute (0 si illisible)."""
    n = 0.0
    for el, num in re.findall(r"([A-Z][a-z]?)(\d*\.?\d*)", formula or ""):
        if el in ("H", "D", "T") or el not in F0:
            continue
        try:
            n += float(num) if num else 1.0
        except ValueError:
            n += 1.0
    return n


def _content_count(ops, atoms):
    """Nombre d'atomes autres que H dans la maille (somme des occupations après symétrie)."""
    els, xyz, occ, _ = expand_atoms(ops, [a for a in atoms if a[0] != "H"])
    return float(np.sum(occ)) if len(els) else 0.0


def _op_keys(trip):
    k = tuple(trip)
    if k not in _KEY_CACHE:
        _KEY_CACHE[k] = {_op_key(parse_symop(t)) for t in trip}
    return _KEY_CACHE[k]



def ops_from_symbol(hall, hm, cell):
    """Opérateurs de symétrie déduits du symbole de Hall ou d'Hermann-Mauguin quand le CIF ne les
    donne pas (module gemmi nécessaire). -> liste de triplets 'x,y,z' ou None."""
    key = (hall, hm, round(cell[3]), round(cell[5]))
    if key in _HALL_CACHE:
        return _HALL_CACHE[key]
    out = None
    try:
        import gemmi
    except ImportError:
        return None
    for h in ([hall] if hall else []):
        h = h.replace("_", " ").strip().strip('"').strip()
        if h and h != "?":
            try:
                g = gemmi.symops_from_hall(h)
                out = [op.triplet() for op in g]
                break
            except Exception:
                pass
    if out is None and hm and hm.strip() not in ("?", ""):
        name = hm.strip()
        tries = [name, name.replace("_", " ")]
        if " " not in name:
            tries.append(name)
        # ancienne notation des groupes cubiques : m 3 m -> m -3 m, m 3 -> m -3
        tries.append(re.sub(r"(?<=[A-Za-z] )3(?= |$)", "-3", name) if re.search(r"[mnda] 3( |$)", name) else name)
        tries.append(re.sub(r"\s+", " ", name.replace(":", " :")))
        for t in tries:
            try:
                sg = gemmi.find_spacegroup_by_name(t, alpha=cell[3], gamma=cell[5])
            except Exception:
                sg = None
            if sg is not None:
                out = [op.triplet() for op in sg.operations()]
                break
    _HALL_CACHE[key] = out
    return out


def parse_symop(s):
    """'-x,1/2+y,-z' -> (R 3x3, t 3)."""
    R = np.zeros((3, 3))
    t = np.zeros(3)
    parts = s.lower().replace(" ", "").split(",")
    if len(parts) != 3:
        raise ValueError("opérateur de symétrie illisible : %r" % s)
    for i, p in enumerate(parts):
        for m in re.finditer(r"([+-]?)(?:(\d+\.?\d*|\.\d+)(?:/(\d+))?\*?)?([xyz])|([+-]?)(\d+\.?\d*|\.\d+)(?:/(\d+))?", p):
            if m.group(4):
                c = float(m.group(2)) / float(m.group(3) or 1) if m.group(2) else 1.0
                R[i, "xyz".index(m.group(4))] += -c if m.group(1) == "-" else c
            elif m.group(6):
                c = float(m.group(6)) / float(m.group(7) or 1)
                t[i] += -c if m.group(5) == "-" else c
    snap = np.rint(t * 24) / 24         # translations écrites en décimales tronquées (0.3333, 0.6667)
    t = np.where(np.abs(t - snap) < 2e-3, snap, t)
    return R, t


SYMOP_TAGS = ("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz", "_space_group_symop.operation_xyz")


def structure_from_cif(text):
    """Renvoie un dict structure : cell(6), ops[(R,t)], atoms[(el, x, y, z, occ, uiso|None)], meta."""
    best = None
    for name, tags, loops in parse_cif_blocks(text):
        cell = [_num(tags.get("_cell_length_" + k)) for k in "abc"] + [_num(tags.get("_cell_angle_" + k)) for k in ("alpha", "beta", "gamma")]
        if any(v is None for v in cell):
            continue
        ops = []
        ops_xyz = []
        atoms = []
        labels = []
        allowed = formula_elements(tags.get("_chemical_formula_sum")) or None
        for lp in loops:
            for tg in SYMOP_TAGS:
                if tg in lp:
                    ops = [parse_symop(s) for s in lp[tg]]
                    ops_xyz = list(lp[tg])
            if "_atom_site_fract_x" in lp:
                nrow = len(lp["_atom_site_fract_x"])
                for r in range(nrow):
                    xyz = [_num(lp["_atom_site_fract_" + k][r]) for k in "xyz"]
                    if any(v is None for v in xyz):
                        continue
                    if lp.get("_atom_site_calc_flag", [""] * nrow)[r].lower() == "dum":
                        continue
                    sym = lp.get("_atom_site_type_symbol", lp.get("_atom_site_label", ["?"] * nrow))[r]
                    el = element_of(sym, allowed) or element_of(lp.get("_atom_site_label", ["?"] * nrow)[r], allowed)
                    if el is None:
                        continue
                    occ = _num(lp.get("_atom_site_occupancy", ["1"] * nrow)[r], 1.0)
                    u = _num(lp.get("_atom_site_u_iso_or_equiv", [None] * nrow)[r])
                    if u is None:
                        b = _num(lp.get("_atom_site_b_iso_or_equiv", [None] * nrow)[r])
                        u = b / (8 * np.pi ** 2) if b is not None else None
                    atoms.append((el, xyz[0], xyz[1], xyz[2], occ, u))
                    labels.append(lp.get("_atom_site_label", [el] * nrow)[r])
        if not ops:
            for tg in SYMOP_TAGS:
                if tg in tags:
                    ops = [parse_symop(tags[tg])]
                    ops_xyz = [tags[tg]]
        sg = tags.get("_space_group_name_h-m_alt") or tags.get("_symmetry_space_group_name_h-m") or tags.get("_space_group.name_h-m_alt") or ""
        sym_from = ""
        hall = tags.get("_space_group_name_hall") or tags.get("_symmetry_space_group_name_hall") or ""
        if not ops and sg.replace(" ", "").upper() not in ("P1",):
            trip = ops_from_symbol(hall, sg, cell)
            if trip:
                ops = [parse_symop(t) for t in trip]
                ops_xyz = list(trip)
                sym_from = "symbole"
        elif ops and (hall or sg):
            # liste d'opérateurs incomplète (centrage ou inversion oubliés) : le symbole fait foi
            # quand ses opérateurs contiennent strictement ceux du CIF
            trip = ops_from_symbol(hall, sg, cell)
            cifkeys = {_op_key(o) for o in ops} if trip else set()
            if trip and cifkeys != _op_keys(trip):
                use = False
                if len(ops) <= 48 and not _is_group(ops):
                    use = True          # liste du CIF incohérente (ce n'est pas un groupe)
                elif len(trip) > len(ops) and cifkeys < _op_keys(trip):
                    # sous-groupe du groupe annoncé : liste incomplète, ou symbole erroné ?
                    # le contenu de la maille comparé à Z x formule tranche ; sans cela le CIF fait foi
                    want = _formula_count(tags.get("_chemical_formula_sum")) * (_num(tags.get("_cell_formula_units_z")) or 0)
                    if want > 0 and atoms:
                        n_cif = _content_count(ops, atoms)
                        n_sym = _content_count([parse_symop(t) for t in trip], atoms)
                        use = abs(n_sym - want) < abs(n_cif - want)
                if use:
                    ops = [parse_symop(t) for t in trip]
                    ops_xyz = list(trip)
                    sym_from = "symbole (liste d'opérateurs du CIF incomplète ou incohérente)"
        nosym = not ops
        if nosym:
            ops = [(np.eye(3), np.zeros(3))]
        meta = {
            "name": tags.get("_chemical_name_systematic") or tags.get("_chemical_name_common") or tags.get("_chemical_name_mineral") or "",
            "formula": tags.get("_chemical_formula_sum") or "",
            "sg": sg.strip(),
            "temp_K": _num(tags.get("_diffrn_ambient_temperature") or tags.get("_cell_measurement_temperature")),
            "rfac": _num(tags.get("_refine_ls_r_factor_gt") or tags.get("_refine_ls_r_factor_all")),
            "journal": tags.get("_journal_name_full") or "",
            "year": _num(tags.get("_journal_year")),
            "doi": tags.get("_journal_paper_doi") or tags.get("_citation_doi") or "",
            "block": name,
            "nosym": nosym and sg.replace(" ", "").upper() not in ("P1", ""),
        }
        if sym_from:
            meta["sym_from"] = sym_from
        for k in ("name", "formula", "journal", "doi"):
            if meta[k] in ("?", "."):
                meta[k] = ""
        authors = next((lp["_publ_author_name"] for lp in loops if "_publ_author_name" in lp), None)
        if authors is None:
            authors = [tags["_publ_author_name"]] if tags.get("_publ_author_name") else []
        st = {"cell": cell, "ops": ops, "atoms": atoms, "meta": meta, "labels": labels, "ops_xyz": ops_xyz or ["x,y,z"],
              "tags": tags, "authors": authors}
        if atoms:
            return st
        best = best or st
    if best is None:
        raise ValueError("pas de maille dans le CIF")
    return best


# ---- géométrie --------------------------------------------------------------------------
def metric(cell):
    a, b, c, al, be, ga = cell
    ca, cb, cg = np.cos(np.radians([al, be, ga]))
    G = np.array([[a * a, a * b * cg, a * c * cb], [a * b * cg, b * b, b * c * ca], [a * c * cb, b * c * ca, c * c]])
    return G


def cell_volume(cell):
    return float(np.sqrt(abs(np.linalg.det(metric(cell)))))


def gen_hkl(cell, dmin):
    """Demi-sphère d'indices (une réflexion par paire de Friedel) avec d >= dmin."""
    G = metric(cell)
    Gs = np.linalg.inv(G)
    # bornes : |h| <= a_i / dmin projeté (longueur réelle / dmin)
    hmax = [int(np.floor(np.sqrt(G[i, i]) / dmin)) for i in range(3)]
    h = np.arange(0, hmax[0] + 1)
    k = np.arange(-hmax[1], hmax[1] + 1)
    l = np.arange(-hmax[2], hmax[2] + 1)
    H, K, L = np.meshgrid(h, k, l, indexing="ij")
    hkl = np.stack([H.ravel(), K.ravel(), L.ravel()], 1)
    keep = (hkl[:, 0] > 0) | ((hkl[:, 0] == 0) & (hkl[:, 1] > 0)) | ((hkl[:, 0] == 0) & (hkl[:, 1] == 0) & (hkl[:, 2] > 0))
    hkl = hkl[keep]
    q = np.einsum("ij,jk,ik->i", hkl, Gs, hkl)
    ok = q <= 1.0 / dmin ** 2 + 1e-12
    return hkl[ok], 1.0 / np.sqrt(q[ok])


def expand_atoms(ops, atoms):
    """Applique la symétrie et retire les doublons. -> els, xyz, occ, uiso.
    Doublons retirés : positions spéciales, mais aussi atomes répétés dans la liste d'origine
    (les structures polymériques de la CSD sont livrées avec des copies par symétrie ou par
    translation de certains atomes). Deux atomes du même élément au même endroit ne comptent
    qu'une fois ; des éléments différents sur un même site (occupation mixte) sont conservés."""
    els, pos, occ, uis = [], [], [], []
    Rs = np.array([o[0] for o in ops])
    ts = np.array([o[1] for o in ops])
    seen = {}
    for el, x, y, z, oc, u in atoms:
        p = (Rs @ np.array([x, y, z]) + ts) % 1.0
        key = np.rint(p * 1000).astype(int) % 1000
        _, idx = np.unique(key, axis=0, return_index=True)
        s = seen.setdefault(el, set())
        for j in idx:
            k = tuple(key[j])
            if k in s:
                continue
            s.add(k)
            els.append(el)
            pos.append(p[j])
            occ.append(oc)
            uis.append(u)
    return els, np.array(pos).reshape(-1, 3), np.array(occ), uis


def unique_atoms(ops, atoms):
    """Indices des atomes à garder : ceux dont l'orbite n'est pas déjà couverte par un atome précédent du même élément."""
    Rs = np.array([o[0] for o in ops])
    ts = np.array([o[1] for o in ops])
    seen, keep = {}, []
    for i, (el, x, y, z, oc, u) in enumerate(atoms):
        key = np.rint(((Rs @ np.array([x, y, z]) + ts) % 1.0) * 1000).astype(int) % 1000
        s = seen.setdefault(el, set())
        ks = {tuple(k) for k in key}
        if ks & s:
            continue
        s |= ks
        keep.append(i)
    return keep


def f0_table(elements, s2):
    """elements : liste de symboles ; s2 = (sinθ/λ)² ; -> dict el -> f(s2)."""
    out = {}
    for el in set(elements):
        c = F0.get(el) or F0["C"]
        a, b, cc = np.array(c[1:5]), np.array(c[5:9]), c[9]
        out[el] = (a[None, :] * np.exp(-b[None, :] * s2[:, None])).sum(1) + cc
    return out


def structure_factors(st, dmin, b_overall=B_OVERALL, use_adp=False, budget=6.0e7, max_hkl=None):
    """|F|² pour chaque hkl de la demi-sphère. -> hkl, d, F2, dmin_effectif.
    max_hkl : plafond du nombre de réflexions (dmin relevé pour les grandes mailles)."""
    els, xyz, occ, uis = expand_atoms(st["ops"], st["atoms"])
    nat = len(els)
    if nat == 0:
        raise ValueError("pas de coordonnées atomiques")
    vol = cell_volume(st["cell"])
    if not (5.0 < vol < 5.0e6):
        raise ValueError("volume de maille invraisemblable")
    # nb de réflexions ~ (2π/3) V / d³ : on relève dmin si la maille est grande
    if max_hkl and 2.094 * vol / dmin ** 3 > max_hkl:
        dmin = float((2.094 * vol / max_hkl) ** (1.0 / 3.0))
    if 2.094 * vol / dmin ** 3 * nat > budget:
        dmin = float((2.094 * vol * nat / budget) ** (1.0 / 3.0))
    hkl, d = gen_hkl(st["cell"], dmin)
    if len(hkl) == 0:
        raise ValueError("aucune réflexion")
    s2 = 1.0 / (4.0 * d * d)
    ftab = f0_table(els, s2)
    FA = np.zeros(len(hkl))
    FB = np.zeros(len(hkl))
    el_arr = np.array(els)
    h32 = (2 * np.pi * hkl).astype(np.float32)
    for el in sorted(set(els)):
        m = el_arr == el
        p = xyz[m].astype(np.float32)
        w = occ[m].astype(np.float32)
        if use_adp:
            u = np.array([(x if x is not None else b_overall / (8 * np.pi ** 2)) for x, mm in zip(uis, m) if mm], np.float32)
        ca = np.zeros(len(hkl), np.float32)
        sa = np.zeros(len(hkl), np.float32)
        step = max(1, int(2.0e6 // max(1, len(p))))
        for i0 in range(0, len(hkl), step):
            ph = h32[i0:i0 + step] @ p.T
            if use_adp:
                ww = w[None, :] * np.exp(-8 * np.pi ** 2 * u[None, :] * s2[i0:i0 + step, None].astype(np.float32))
                ca[i0:i0 + step] = (np.cos(ph) * ww).sum(1)
                sa[i0:i0 + step] = (np.sin(ph) * ww).sum(1)
            else:
                ca[i0:i0 + step] = np.cos(ph) @ w
                sa[i0:i0 + step] = np.sin(ph) @ w
        FA += ftab[el] * ca
        FB += ftab[el] * sa
    f2 = FA * FA + FB * FB
    if not use_adp:
        f2 = f2 * np.exp(-2 * b_overall * s2)
    return hkl, d, f2, dmin


def lp_factor(d, lam):
    """Lorentz-polarisation (faisceau non polarisé), poudre, longueur d'onde constante."""
    st = np.clip(lam / (2.0 * np.asarray(d, float)), 1e-6, 0.999999)
    th = np.arcsin(st)
    return (1 + np.cos(2 * th) ** 2) / (np.sin(th) ** 2 * np.cos(th))


def merged_lines(hkl, d, f2):
    """Regroupe les réflexions de même d (à un pas de grille près).
    -> bins, I (somme des F², sphère entière), hkl de la réflexion la plus forte de chaque groupe."""
    b = d_to_bin(d)
    ok = (b >= 0) & (b < N_BINS)
    b, f2, hkl = b[ok], f2[ok], hkl[ok]
    if len(b) == 0:
        return np.zeros(0, np.int64), np.zeros(0), np.zeros((0, 3), np.int64)
    acc = np.bincount(b, weights=2.0 * f2, minlength=N_BINS)
    nz = np.nonzero(acc)[0]
    # étiquette de groupe : cases adjacentes fusionnées par paquets de 3 au plus
    lab = np.zeros(len(nz), np.int64)
    start = 0
    for i in range(1, len(nz)):
        if nz[i] - nz[i - 1] > 1 or nz[i] - nz[start] >= 3:
            start = i
        lab[i] = start
    grp_of_bin = np.full(N_BINS, -1, np.int64)
    grp_of_bin[nz] = lab
    uniq, inv = np.unique(lab, return_inverse=True)
    I = np.bincount(inv, weights=acc[nz])
    cen = np.bincount(inv, weights=acc[nz] * nz) / np.maximum(I, 1e-300)
    g = np.searchsorted(uniq, grp_of_bin[b])
    order = np.lexsort((-f2, g))
    first = np.r_[True, g[order][1:] != g[order][:-1]]
    rep = np.zeros((len(uniq), 3), np.int64)
    rep[g[order][first]] = hkl[order][first]
    return np.rint(cen).astype(np.int64), I, rep


def stored_lines(st, n_store=N_STORE):
    """Raies conservées dans l'index : les n plus fortes (classement avec LP Cu Kα).
    -> bins uint16, intens float16 (somme m·F², 1 pour la plus forte), hkl int8 (n,3), dmin_effectif."""
    hkl, d, f2, dmin_eff = structure_factors(st, D_MIN_STORE, budget=1.5e7, max_hkl=MAX_HKL_STORE)
    b, inten, rep = merged_lines(hkl, d, f2)
    if len(b) == 0:
        raise ValueError("aucune raie")
    rank = inten * lp_factor(bin_to_d(b), CU_KA1)
    o = np.argsort(-rank)[:n_store]
    b, inten, rep = b[o], inten[o], rep[o]
    inten = inten / inten.max()
    return b.astype(np.uint16), inten.astype(np.float16), np.clip(rep, -127, 127).astype(np.int8), dmin_eff


def pack_lines(b, inten, hkl):
    return b.astype("<u2").tobytes() + inten.astype("<f2").tobytes() + hkl.astype("i1").tobytes()


def unpack_lines(blob):
    n = len(blob) // 7
    return (np.frombuffer(blob, "<u2", n), np.frombuffer(blob, "<f2", n, offset=2 * n),
            np.frombuffer(blob, "i1", 3 * n, offset=4 * n).reshape(n, 3))


def elements_of(st):
    seen = []
    for a in st["atoms"]:
        if a[0] not in seen:
            seen.append(a[0])
    return sorted(seen, key=lambda e: F0[e][0])


def formula_from_atoms(st):
    """Contenu de la maille (somme des occupations après symétrie)."""
    els, xyz, occ, _ = expand_atoms(st["ops"], st["atoms"])
    tot = {}
    for e, o in zip(els, occ):
        tot[e] = tot.get(e, 0.0) + o
    order = [e for e in ("C", "H") if e in tot] + sorted(e for e in tot if e not in ("C", "H"))
    return " ".join("%s%g" % (e, round(tot[e], 2)) for e in order)


def write_cif(path, name, cell, ops_xyz, atoms, sg="", extra=None):
    """Écrit un CIF minimal. ops_xyz : liste de chaînes 'x,y,z' ; atoms : (label, el, x, y, z, occ, uiso)."""
    L = ["data_%s" % re.sub(r"\W", "_", name)]
    for k, v in (extra or {}).items():
        L.append("%s '%s'" % (k, v))
    for k, v in zip(("_cell_length_a", "_cell_length_b", "_cell_length_c", "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"), cell):
        L.append("%s %.5f" % (k, v))
    if sg:
        L.append("_symmetry_space_group_name_H-M '%s'" % sg)
    L += ["loop_", "_symmetry_equiv_pos_as_xyz"] + ["'%s'" % o for o in ops_xyz]
    L += ["loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z", "_atom_site_occupancy", "_atom_site_U_iso_or_equiv"]
    for lab, el, x, y, z, occ, u in atoms:
        L.append("%s %s %.5f %.5f %.5f %.3f %.4f" % (lab, el, x, y, z, occ, u))
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")


def clean_cif(txt, cell=None, title="", comment=""):
    """Réécrit un CIF minimal et propre pour l'affinement : maille donnée (sinon celle du CIF), opérateurs de
    symétrie, atomes sans doublons, Uiso raisonnable quand le CIF n'en donne pas (0,038 ; 0,05 pour H)."""
    st = structure_from_cif(txt)
    keep = unique_atoms(st["ops"], st["atoms"])
    m = st["meta"]
    L = ["# " + c for c in comment.splitlines() if c]
    L.append("data_%s" % (re.sub(r"\W", "_", title or m.get("block") or "phase")))
    if m.get("name"):
        L.append("_chemical_name_common '%s'" % m["name"].replace("'", " "))
    for k, v in zip(("_cell_length_a", "_cell_length_b", "_cell_length_c", "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"), cell or st["cell"]):
        L.append("%s %.4f" % (k, v))
    if m.get("sg"):
        L.append("_symmetry_space_group_name_H-M '%s'" % m["sg"])
    L += ["loop_", "_symmetry_equiv_pos_as_xyz"] + ["'%s'" % o.strip() for o in st["ops_xyz"]]
    L += ["loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z",
          "_atom_site_occupancy", "_atom_site_U_iso_or_equiv"]
    for i in keep:
        el, x, y, z, occ, u = st["atoms"][i]
        if u is None or u <= 0:
            u = 0.05 if el == "H" else 0.038
        L.append("%s %s %.5f %.5f %.5f %.3f %.4f" % (re.sub(r"\s", "", st["labels"][i]) or el, el, x, y, z, occ, u))
    return "\n".join(L) + "\n", len(st["atoms"]) - len(keep)
