"""Lecture d'un histogramme de poudre dans un projet GSAS-II (.gpx) sans importer GSAS-II.
Le .gpx est une suite de pickles : seules les classes numpy et les types de base sont acceptés."""
import pickle
import numpy as np

_OK = {("builtins", n) for n in ("set", "frozenset", "list", "dict", "tuple", "complex", "slice", "range", "bytearray", "object")}
_OK |= {("collections", "OrderedDict"), ("collections", "defaultdict"), ("copyreg", "_reconstructor"), ("_codecs", "encode")}
_NO = {"load", "loads", "save", "savez", "fromfile", "memmap", "DataSource"}


class _Stub:
    def __init__(self, *a, **k):
        pass

    def __setstate__(self, s):
        pass


class _Unpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if (module, name) in _OK or (module.split(".")[0] == "numpy" and name not in _NO):
            return super().find_class(module, name)
        return _Stub


def load_gpx(path):
    proj, order = {}, []
    with open(path, "rb") as f:
        while True:
            try:
                data = _Unpickler(f, encoding="latin-1").load()
            except EOFError:
                break
            proj[data[0][0]] = {"data": data[0][1], **{s[0]: s[1] for s in data[1:]}}
            order.append(data[0][0])
    return proj, order


def powder_histograms(path):
    return [n for n in load_gpx(path)[1] if n.startswith("PWDR")]


def pattern_from_parts(name, data, inst, limits):
    """Construit le dict du diagramme à partir des éléments de l'arbre GSAS-II :
    data = [x, yobs, ...], inst = dict des paramètres instrumentaux, limits = (bas, haut)."""
    x = np.ma.getdata(np.array(data[0])).astype(float)
    y = np.ma.getdata(np.array(data[1])).astype(float)
    typ = str(inst["Type"][0])
    if not (typ.startswith("PX") or typ.startswith("PN")) or "T" in typ[2:]:
        raise ValueError("histogramme de type %s non pris en charge (longueur d'onde constante uniquement)" % typ)
    if "Lam1" in inst:
        lam1, lam2, ratio = inst["Lam1"][1], inst["Lam2"][1], inst["I(L2)/I(L1)"][1]
    else:
        lam1, lam2, ratio = inst["Lam"][1], None, 0.0
    lo, hi = limits
    return {"name": name, "x": x, "y": y, "lam1": float(lam1), "lam2": lam2, "ratio": float(ratio), "limits": (float(lo), float(hi))}


def read_pattern(path, hist=None):
    """-> dict : name, x (2θ °), y, lam1, lam2, ratio (Kα2/Kα1), limits."""
    proj, order = load_gpx(path)
    names = [n for n in order if n.startswith("PWDR")]
    if not names:
        raise ValueError("aucun histogramme de poudre dans %s" % path)
    name = names[0] if hist is None else (names[hist] if isinstance(hist, int) else hist)
    h = proj[name]
    return pattern_from_parts(name, h["data"][1], h["Instrument Parameters"][0], h["Limits"][1])
