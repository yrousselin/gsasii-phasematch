#!/usr/bin/env python3
"""pm_search : recherche de phases à l'aveugle dans un mélange, à partir d'un projet GSAS-II.

  python pm_search.py projet.gpx --db index.sqlite [--hist N] [--max-phases 4]
         [--sources CSD,local] [--elements "C H N O"] [--strict] [--out DOSSIER]

Principe (répété sur le résidu tant qu'une phase est trouvée) :
  1. tri grossier de tout l'index : les raies fortes de chaque référence sont-elles présentes
     dans le diagramme (tolérance large, dilatation isotrope de maille autorisée) ?
  2. tri fin des meilleures : maille ajustée sur les pics observés (à partir des hkl stockés),
     puis comparaison au profil avec une tolérance serrée ;
  3. validation des premières avec la structure complète relue dans sa source
     (maille, échelle et largeur ajustées sur le profil) ; la meilleure est soustraite.
"""
import argparse, csv, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
from scipy.optimize import least_squares
from scipy.ndimage import maximum_filter1d
from scipy.signal import find_peaks, peak_widths

import pm_core as pc
import pm_index as pi
import pm_gpx

LN2 = np.log(2.0)


# ---- préparation du diagramme ------------------------------------------------------------
def snip(y, m):
    v = np.log(np.log(np.sqrt(np.maximum(y, 0) + 1) + 1) + 1)
    for p in range(m, 0, -1):
        a = np.r_[np.full(p, v[0]), v[:-p]]
        b = np.r_[v[p:], np.full(p, v[-1])]
        v = np.minimum(v, 0.5 * (a + b))
    return (np.exp(np.exp(v) - 1) - 1) ** 2 - 1


class Obs:
    def __init__(self, pat, tmin=None, tmax=None):
        lo, hi = pat["limits"]
        lo, hi = max(lo, tmin or lo), min(hi, tmax or hi)
        m = (pat["x"] >= lo) & (pat["x"] <= hi)
        self.x, self.yo = pat["x"][m], pat["y"][m]
        self.lam, self.lam2, self.ratio = pat["lam1"], pat["lam2"], pat["ratio"]
        step = float(np.median(np.diff(self.x)))
        self.bkg = snip(self.yo, max(5, int(round(1.2 / step))))
        self.net = self.yo - self.bkg
        self.sig = np.sqrt(np.maximum(self.yo, 1.0))
        self.res = np.clip(self.net, 0, None)          # résidu courant (2θ)
        self.ss0 = float((self.res ** 2).sum())
        d_lo = self.lam / (2 * np.sin(np.radians(self.x[-1] / 2)))
        d_hi = self.lam / (2 * np.sin(np.radians(self.x[0] / 2)))
        self.bmin = int(max(0, pc.d_to_bin(d_lo)))
        self.bmax = int(min(pc.N_BINS - 1, pc.d_to_bin(d_hi)))
        self.gb = np.arange(self.bmin, self.bmax + 1)
        self.gt = 2 * np.degrees(np.arcsin(np.clip(self.lam / (2 * pc.bin_to_d(self.gb)), 0, 1)))   # 2θ de chaque case
        self.w0, self.w1 = self._widths()

    def _widths(self):
        pk, pr = find_peaks(self.res, height=8 * np.median(self.sig), prominence=8 * np.median(self.sig))
        if len(pk) < 3:
            return 0.08, 0.0
        step = float(np.median(np.diff(self.x)))
        w = peak_widths(self.res, pk, rel_height=0.5)[0] * step
        t = np.tan(np.radians(self.x[pk] / 2))
        o = np.argsort(t)
        tt, ww = [], []
        for part in np.array_split(o, min(4, max(1, len(o) // 4))):
            tt.append(np.median(t[part]))
            ww.append(np.percentile(w[part], 25))
        if len(tt) >= 2:
            w1, w0 = np.polyfit(tt, ww, 1)
            if w1 < 0 or w0 < 0.02:
                w0, w1 = float(np.median(ww)), 0.0
        else:
            w0, w1 = ww[0], 0.0
        return float(np.clip(w0, 0.03, 0.4)), float(np.clip(w1, 0.0, 0.5))

    def fwhm(self, tth):
        return self.w0 + self.w1 * np.tan(np.radians(np.asarray(tth) / 2))

    def grid(self):
        """Résidu sur la grille ln(d), seuil de bruit, liste de pics (cases, hauteurs)."""
        o = np.interp(self.gt[::-1], self.x, self.res)[::-1]
        thr = 4.0 * np.interp(self.gt[::-1], self.x, self.sig)[::-1]
        pk, pr = find_peaks(o, height=thr, prominence=0.5 * thr)
        pos = pk.astype(float)
        for i, k in enumerate(pk):
            if 0 < k < len(o) - 1:
                den = o[k - 1] - 2 * o[k] + o[k + 1]
                if den < 0:
                    pos[i] = k + 0.5 * (o[k - 1] - o[k + 1]) / den
        return o, thr, pos + self.bmin, o[pk]


# ---- étape 1 : tri grossier de tout l'index ----------------------------------------------
def stage1(ix, rows, obs, o, thr, tol_bins, shifts, K=12, q=0.35, c_over=2.0, chunk=50000):
    A = np.zeros(pc.N_BINS, np.float32)
    A[obs.bmin:obs.bmax + 1] = maximum_filter1d(np.where(o > thr, o, 0.0), size=2 * tol_bins + 1, mode="constant")
    lp = pc.lp_factor(pc.bin_to_d(np.arange(pc.N_BINS)), obs.lam).astype(np.float32)
    n = len(rows)
    best = np.zeros(n, np.float32)
    bshift = np.zeros(n, np.int32)
    for i0 in range(0, n, chunk):
        r = rows[i0:i0 + chunk]
        B = ix["bins"][r]
        valid = B >= 0
        Bc = np.clip(B, 0, pc.N_BINS - 1)
        R0 = ix["inten"][r] * lp[Bc]
        ar = np.arange(len(r))
        for k in shifts:
            Bk = B + k
            inr = valid & (Bk >= obs.bmin) & (Bk <= obs.bmax)
            M = A[np.clip(Bk, 0, pc.N_BINS - 1)] * inr
            R = R0 * inr
            top = inr[:, :K] & (R[:, :K] > 0)
            ratio = np.where(top, M[:, :K] / np.maximum(R[:, :K], 1e-30), np.inf)
            nin = top.sum(1)
            ratio.sort(axis=1)
            s = ratio[ar, np.minimum((q * nin).astype(int), K - 1)]
            s[(nin < 4) | ~np.isfinite(s)] = 0.0
            pred = s[:, None] * R
            E = np.minimum(pred, M).sum(1)
            over = np.maximum(pred - c_over * M, 0).sum(1)
            sc = E * (1.0 - over / np.maximum(pred.sum(1), 1e-30)) ** 2
            better = sc > best[i0:i0 + chunk]
            best[i0:i0 + chunk][better] = sc[better]
            bshift[i0:i0 + chunk][better] = k
    return best, bshift


# ---- maille : paramètres du tenseur réciproque -------------------------------------------
def design(hkl):
    h, k, l = [hkl[:, i].astype(float) for i in range(3)]
    return np.stack([h * h, k * k, l * l, h * k, h * l, k * l], 1)


def cell_to_A(cell):
    Gs = np.linalg.inv(pc.metric(cell))
    return np.array([Gs[0, 0], Gs[1, 1], Gs[2, 2], 2 * Gs[0, 1], 2 * Gs[0, 2], 2 * Gs[1, 2]])


def A_to_cell(A):
    Gs = np.array([[A[0], A[3] / 2, A[4] / 2], [A[3] / 2, A[1], A[5] / 2], [A[4] / 2, A[5] / 2, A[2]]])
    if np.any(np.linalg.eigvalsh(Gs) <= 0):
        return None
    G = np.linalg.inv(Gs)
    a, b, c = np.sqrt(np.diag(G))
    al = np.degrees(np.arccos(np.clip(G[1, 2] / (b * c), -1, 1)))
    be = np.degrees(np.arccos(np.clip(G[0, 2] / (a * c), -1, 1)))
    ga = np.degrees(np.arccos(np.clip(G[0, 1] / (a * b), -1, 1)))
    return [float(a), float(b), float(c), float(al), float(be), float(ga)]


def sym_basis(cell):
    """Matrice C (6, p) : A = C·p, selon le système cristallin déduit de la maille."""
    a, b, c, al, be, ga = cell
    eq = lambda u, v: abs(u - v) < 2e-4 * max(u, v)
    r = lambda u, v=90.0: abs(u - v) < 0.01
    E = np.eye(6)
    if r(al) and r(be) and r(ga):
        if eq(a, b) and eq(b, c):
            return (E[0] + E[1] + E[2])[:, None]
        for i, j, k in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):
            if eq(cell[i], cell[j]):
                return np.stack([E[i] + E[j], E[k]], 1)
        return E[:, :3]
    if r(al) and r(be) and r(ga, 120.0) and eq(a, b):
        return np.stack([E[0] + E[1] + E[3], E[2]], 1)
    if eq(a, b) and eq(b, c) and abs(al - be) < 0.01 and abs(be - ga) < 0.01:
        return np.stack([E[0] + E[1] + E[2], E[3] + E[4] + E[5]], 1)
    n90 = [r(al), r(be), r(ga)]
    if sum(n90) == 2:
        cross = {0: 5, 1: 4, 2: 3}[n90.index(False)]
        return np.stack([E[0], E[1], E[2], E[cross]], 1)
    return E


def refine_cell_on_peaks(hkl, bline, wline, cell0, shift_bins, pk_bins, tols=(16, 8, 5)):
    """Ajuste la maille d'une référence sur les pics observés. -> (p, C, corr, nouvelles cases des raies, n appariés)."""
    D = design(hkl)
    A0 = cell_to_A(cell0)
    Qh0 = D @ A0
    ok = (Qh0 > 0) & (bline >= 0)
    Qline = 1.0 / pc.bin_to_d(bline) ** 2
    corr = np.where(ok, Qline / np.where(Qh0 > 0, Qh0, 1.0), 1.0)
    C = sym_basis(cell0)
    p0 = np.linalg.lstsq(C, A0, rcond=None)[0]
    eps = np.exp(shift_bins * pc.LOG_STEP) - 1.0
    p = p0 / (1 + eps) ** 2
    pstart = p.copy()
    diag = np.mean(np.abs(C.T @ np.array([1, 1, 1, 0, 0, 0.0]) * pstart)) or 1.0
    scale = np.maximum(np.abs(pstart), 0.05 * abs(diag))
    strong = ok & (wline >= 0.03 * wline.max())
    DC = (D @ C) * corr[:, None]
    nm = 0
    if len(pk_bins) == 0:
        return p, C, corr, bline + shift_bins, 0
    for tol in tols:
        Q = DC @ p
        bpred = np.where(Q > 0, np.log(np.maximum(Q, 1e-12) ** -0.5 / pc.D_MIN_STORE) / pc.LOG_STEP, -1e9)
        j = np.clip(np.searchsorted(pk_bins, bpred), 1, len(pk_bins) - 1)
        near = np.where(np.abs(pk_bins[j - 1] - bpred) < np.abs(pk_bins[j] - bpred), j - 1, j)
        m = strong & (np.abs(pk_bins[near] - bpred) <= tol)
        nm = int(m.sum())
        if nm < C.shape[1] + 2:
            break
        Qobs = 1.0 / pc.bin_to_d(pk_bins[near[m]]) ** 2
        rows = DC[m] / Qobs[:, None] / 5e-4
        rhs = np.full(nm, 1.0 / 5e-4)
        prior = np.diag(1.0 / (0.02 * scale))
        sol = np.linalg.lstsq(np.vstack([rows, prior]), np.r_[rhs, pstart / (0.02 * scale)], rcond=None)[0]
        if A_to_cell(C @ sol) is None:
            break
        p = sol
    Q = DC @ p
    bnew = np.where(ok & (Q > 0), np.log(np.maximum(Q, 1e-12) ** -0.5 / pc.D_MIN_STORE) / pc.LOG_STEP, -1.0)
    return p, C, corr, bnew, nm


def asym_fit(o, c, kappa=4.0, it=4):
    cc = float(c @ c)
    if cc <= 0:
        return 0.0, 0.0
    s = float(o @ c) / cc
    w = np.ones_like(o)
    for _ in range(it):
        w = np.where(s * c > o, kappa, 1.0)
        s = float((w * o * c).sum() / max((w * c * c).sum(), 1e-30))
    return s, 1.0 - float((w * (o - s * c) ** 2).sum()) / max(float(o @ o), 1e-30)


def stage2(ix, cand_rows, cand_shift, obs, o, pk_bins):
    """Pour chaque candidat : maille ajustée sur les pics, score sur le profil (grille ln d)."""
    lp_of = lambda b: pc.lp_factor(pc.bin_to_d(b), obs.lam)
    ng = len(o)
    out = []
    for r, sh in zip(cand_rows, cand_shift):
        B = ix["bins"][r].astype(float)
        valid = B >= 0
        I = ix["inten"][r].astype(float) * np.where(valid, lp_of(np.where(valid, B, 0)), 0.0)
        p, C, corr, bnew, nm = refine_cell_on_peaks(ix["hkl"][r].astype(int), np.where(valid, B, -1), I, list(ix["cell"][r]), int(sh), pk_bins)
        c = np.zeros(ng)
        inr = valid & (bnew >= obs.bmin) & (bnew <= obs.bmax)
        tth = 2 * np.degrees(np.arcsin(np.clip(obs.lam / (2 * pc.bin_to_d(np.where(inr, bnew, obs.bmin))), 0, 1)))
        w = obs.fwhm(tth)
        wb = np.radians(w) / (2 * np.tan(np.radians(tth / 2))) / pc.LOG_STEP       # largeur en cases
        for j in np.nonzero(inr)[0]:
            cen = bnew[j] - obs.bmin
            hw = int(3 * wb[j]) + 2
            i0, i1 = max(0, int(cen) - hw), min(ng, int(cen) + hw + 1)
            u = (np.arange(i0, i1) - cen) / wb[j]
            c[i0:i1] += I[j] / w[j] * np.exp(-4 * LN2 * u * u)
        s, sc = asym_fit(o, c)
        out.append((sc, s, p, C, nm, int(inr.sum()), (bnew, I, inr)))
    return out


# ---- étape 3 : structure complète, ajustement sur le profil ------------------------------
def profile(x, tth, area, fw, eta=0.5):
    y = np.zeros_like(x)
    step = (x[-1] - x[0]) / (len(x) - 1)
    for t, a, g in zip(tth, area, fw):
        i0 = max(0, int((t - 12 * g - x[0]) / step))
        i1 = min(len(x), int((t + 12 * g - x[0]) / step) + 2)
        if i1 <= i0:
            continue
        u = (x[i0:i1] - t) / g
        y[i0:i1] += a / g * (eta * (2 / np.pi) / (1 + 4 * u * u) + (1 - eta) * 2 * np.sqrt(LN2 / np.pi) * np.exp(-4 * LN2 * u * u))
    return y


class Fit3:
    WPOW = 0.25     # poids 1/yo**WPOW dans l'ajustement (0 : aucun ; 0,5 : statistique)

    def __init__(self, obs, st, cell_start):
        self.obs = obs
        self.w = 1.0 / np.maximum(obs.yo, 1.0) ** self.WPOW
        dmin = obs.lam / (2 * np.sin(np.radians(obs.x[-1] / 2))) * 0.995
        hkl, d, f2, _ = pc.structure_factors(st, max(dmin, 0.7))
        I = 2 * f2 * pc.lp_factor(d, obs.lam)
        keep = I > 2e-4 * I.max()
        self.D = design(hkl[keep])
        self.I = I[keep] / I.max()
        self.imax = float(I.max())      # échelle absolue des intensités (pour l'estimation des proportions)
        self.C = sym_basis(cell_start)
        self.p0 = np.linalg.lstsq(self.C, cell_to_A(cell_start), rcond=None)[0]
        diag = np.abs(self.p0).max()
        self.sc = np.maximum(np.abs(self.p0), 0.05 * diag)

    def calc(self, v):
        obs = self.obs
        p = self.p0 + v[:len(self.p0)] * self.sc
        Q = self.D @ (self.C @ p)
        if np.any(Q <= 0):
            return None
        d = Q ** -0.5
        w0, w1 = v[-2], v[-1]
        y = np.zeros_like(obs.x)
        st1 = np.clip(obs.lam / (2 * d), 0, 0.999)
        t1 = 2 * np.degrees(np.arcsin(st1))
        self._st = (t1, self.I * 0.788 / (w0 + w1 * np.tan(np.radians(t1 / 2))))       # hauteur de pic (Kα1)
        for lam, k in ((obs.lam, 1.0),) + (((obs.lam2, obs.ratio),) if obs.lam2 and obs.ratio > 0 else ()):
            st = lam / (2 * d)
            m = st < 0.999
            tth = 2 * np.degrees(np.arcsin(st[m]))
            mm = (tth > obs.x[0] - 1) & (tth < obs.x[-1] + 1)
            fw = w0 + w1 * np.tan(np.radians(tth[mm] / 2))
            y += profile(obs.x, tth[mm], k * self.I[m][mm], fw)
        return y

    def resid(self, v):
        c = self.calc(v)
        w = self.w
        if c is None:
            return self.obs.res * w * 10
        cc = float((w * c) @ (w * c))
        self.s = max(0.0, float((w * self.obs.res) @ (w * c)) / cc) if cc > 0 else 0.0
        return (self.obs.res - self.s * c) * w

    def run(self):
        n = len(self.p0)
        v0 = np.r_[np.zeros(n), max(self.obs.w0, 0.031), max(self.obs.w1, 0.001)]
        lo = np.r_[np.full(n, -0.03), 0.03, 0.0]
        hi = np.r_[np.full(n, 0.03), 0.6, 0.8]
        sol = least_squares(self.resid, v0, bounds=(lo, hi), x_scale=np.r_[np.full(n, 3e-4), 0.02, 0.02], diff_step=np.r_[np.full(n, 2e-5), 1e-3, 1e-3], max_nfev=40)
        r = self.resid(sol.x)
        cell = A_to_cell(self.C @ (self.p0 + sol.x[:n] * self.sc))
        expl = 1.0 - float(r @ r) / max(float((self.w * self.obs.res) @ (self.w * self.obs.res)), 1e-30)
        yc = self.s * self.calc(sol.x)
        t1, h1 = self._st
        m = (t1 >= self.obs.x[0]) & (t1 <= self.obs.x[-1]) & (h1 > 0.01 * h1.max())
        self.sticks = (t1[m], self.s * h1[m])
        return expl, yc, cell, sol.x[-2:]


def same_lines(l1, l2, tol=16, ntop=15, thr=0.7):
    """Deux références donnent-elles les mêmes raies (après ajustement de maille sur le diagramme) ?
    Indépendant du choix des axes : c'est la comparaison qui définit une « famille » (redéterminations,
    autres descriptions de la même maille, isotypes, surstructures)."""
    def frac(a, b):
        ba, ia, ina = a
        bb, _, inb = b
        if ina.sum() < 5 or inb.sum() < 5:
            return 0.0
        idx = np.nonzero(ina)[0]
        idx = idx[np.argsort(-ia[idx])[:ntop]]
        ref = np.sort(bb[inb])
        j = np.clip(np.searchsorted(ref, ba[idx]), 1, len(ref) - 1)
        dist = np.minimum(np.abs(ref[j - 1] - ba[idx]), np.abs(ref[j] - ba[idx]))
        return float((dist <= tol).mean())
    return max(frac(l1, l2), frac(l2, l1)) >= thr


# ---- boucle principale -------------------------------------------------------------------
def select_rows(con, ix, sources=None, elements=None, strict=False):
    if not sources and not elements:
        return np.arange(len(ix["ids"]))
    q, args = "SELECT id, elements FROM entries WHERE status='ok'", []
    if sources:
        q += " AND source IN (%s)" % ",".join("?" * len(sources))
        args += sources
    keep = []
    want = set(elements or [])
    for rid, els in con.execute(q, args):
        have = set((els or "").split())
        if want and not want <= have:
            continue
        if strict and want and not have <= want | {"H"}:
            continue
        keep.append(rid)
    return np.nonzero(np.isin(ix["ids"], np.array(keep, np.int64)))[0]


def search(gpx, db, hist=None, max_phases=4, sources=None, elements=None, strict=False, n2=400, n3=6,
           tol=0.4, max_shift=(-0.6, 1.8), min_gain=0.12, min_contrast=1.5, out=None, verbose=True, log=print, explore=False):
    t00 = time.time()
    pat = pm_gpx.read_pattern(gpx, hist)
    obs = Obs(pat)
    con = pi.open_db(db)
    ix = pi.load_index(db)
    rows = select_rows(con, ix, sources, elements, strict)
    log("Diagramme : %s" % pat["name"])
    log("  2θ %.2f–%.2f°, λ = %.5f Å, largeur de raie estimée %.3f° + %.3f·tanθ" % (obs.x[0], obs.x[-1], obs.lam, obs.w0, obs.w1))
    log("Index : %d entrées utilisables, %d retenues pour la recherche (lecture %.1f s)" % (len(ix["ids"]), len(rows), time.time() - t00))
    tol_bins = int(round(tol / 100 / pc.LOG_STEP))
    step = max(1, tol_bins // 2)
    shifts = np.arange(int(max_shift[0] / 100 / pc.LOG_STEP), int(max_shift[1] / 100 / pc.LOG_STEP) + 1, step)
    found, table = [], []
    for rnd in range(1, max_phases + 1):
        t0 = time.time()
        o, thr, pk_bins, pk_h = obs.grid()
        if len(pk_bins) < 4:
            log("\nTour %d : plus assez de pics dans le résidu (%d), arrêt." % (rnd, len(pk_bins)))
            break
        s1, sh = stage1(ix, rows, obs, o, thr, tol_bins, shifts)
        t1 = time.time()
        order = np.argsort(-s1)[:n2]
        order = order[s1[order] > 0]
        res2 = stage2(ix, rows[order], sh[order], obs, o, pk_bins)
        t2 = time.time()
        sc2 = np.array([r[0] for r in res2]) if res2 else np.zeros(0)
        o2 = np.argsort(-sc2)
        # familles de même maille (redéterminations, isotypes) : 2 validations au plus par famille
        leaders, fam, chosen = [], {}, []
        for k in o2[:80]:
            g = next((i for i, c in enumerate(leaders) if same_lines(c, res2[k][6])), None)
            if g is None:
                leaders.append(res2[k][6])
                g = len(leaders) - 1
            fam.setdefault(g, []).append(k)
            if len(fam[g]) <= 2 and len(chosen) < n3:
                chosen.append((k, g))
        cands = []
        for k, g in chosen:
            r = rows[order[k]]
            row = pi.get_row(con, int(ix["ids"][r]))
            cell2 = A_to_cell(res2[k][3] @ res2[k][2]) or list(ix["cell"][r])
            try:
                st, txt = pi.get_structure(con, row)
                expl, ycalc, cell3, wfit = Fit3(obs, st, cell2).run()
            except Exception as e:
                expl, ycalc, cell3, wfit, txt = -1.0, None, cell2, (0, 0), ""
                log("   (validation impossible pour %s : %s)" % (row["ident"], e))
            cands.append({"row": row, "s1": float(s1[order[k]]), "rank1": int(k) + 1, "s2": float(sc2[k]), "fam": g, "nfam": len(fam[g]),
                          "s3": float(expl), "shift": 100 * (np.exp(sh[order[k]] * pc.LOG_STEP) - 1), "cell": cell3, "ycalc": ycalc, "cif": txt})
        t3 = time.time()
        cands.sort(key=lambda c: -c["s3"])
        ss = float(obs.res @ obs.res)
        log("\nTour %d — signal restant : %.1f %% du signal initial ; %d pics dans le résidu" % (rnd, 100 * ss / obs.ss0, len(pk_bins)))
        log("  durées : tri grossier %.1f s (%d réf.), tri fin %.1f s (%d), validation %.1f s (%d)" % (t1 - t0, len(rows), t2 - t1, len(order), t3 - t2, len(cands)))
        log("  %-3s %-26s %-34s %-11s %8s %5s | %8s %6s %7s %7s  %s" % ("", "identifiant", "nom", "groupe", "V (Å³)", "T(K)", "grossier", "(rang)", "fin", "validé", "mêmes raies"))
        for i, c in enumerate(cands):
            row = c["row"]
            log("  %-3d %-26s %-34s %-11s %8.1f %5s | %8.0f %6d %6.1f%% %6.1f%%  %s" % (
                i + 1, row["ident"][-26:], (row["name"] or "")[:34], (row["sg"] or "")[:11], pc.cell_volume(c["cell"]) if c["cell"] else 0,
                "%.0f" % row["temp_K"] if row["temp_K"] else "-", c["s1"], c["rank1"], 100 * c["s2"], 100 * c["s3"],
                "famille %d (%d entrées dans les 80 premières)" % (c["fam"] + 1, c["nfam"])))
            table.append([rnd, i + 1, row["source"], row["ident"], row["name"], row["formula"], row["sg"], row["temp_K"], c["s1"], c["rank1"],
                          round(c["s2"], 4), round(c["s3"], 4), c["fam"] + 1] + [round(v, 4) for v in (c["cell"] or [0] * 6)])
        if not cands or cands[0]["ycalc"] is None:
            log("  -> aucun candidat : arrêt.")
            break
        # même phase que la première : mêmes raies au tri fin, ou profils ajustés quasi identiques
        y0 = cands[0]["ycalc"]
        n0 = float(np.sqrt(y0 @ y0)) or 1.0
        for c in cands:
            yc = c["ycalc"]
            c["cos"] = float(y0 @ yc) / (n0 * (float(np.sqrt(yc @ yc)) or 1.0)) if yc is not None else 0.0
            c["same"] = c["fam"] == cands[0]["fam"] or c["cos"] >= 0.8
        for c, line in zip(cands, table[-len(cands):]):
            line.append(int(c["same"]))
        same = [c["row"]["ident"] for c in cands[1:] if c["same"]]
        others = [c["s3"] for c in cands if not c["same"]]
        contrast = cands[0]["s3"] / max(max(others) if others else 0.0, 0.02)
        cands[0]["contrast"] = contrast
        if same:
            log("  même phase que le n° 1 (profils identiques) : %s" % ", ".join(same))
        log("  contraste du n° 1 par rapport aux autres phases : %.1f (seuil %.1f) ; part du résidu expliquée : %.0f %% (seuil %.0f %%)" % (
            contrast, min_contrast, 100 * cands[0]["s3"], 100 * min_gain))
        # retenue si elle se détache nettement (contraste >= 2), ou moins nettement (>= min_contrast)
        # mais en expliquant la majeure partie du résidu (>= 60 %)
        ok = cands[0]["s3"] >= min_gain and (contrast >= 2.0 or (contrast >= min_contrast and cands[0]["s3"] >= 0.60))
        if not ok and not explore:
            log("  -> pas de phase qui se détache nettement : arrêt. Meilleur candidat non retenu : %s (%s)" % (cands[0]["row"]["ident"], (cands[0]["row"]["name"] or "")[:50]))
            break
        cands[0]["conf"] = ("forte" if contrast >= 2.0 else "moyenne, à confirmer par affinement") if ok else "NON VALIDÉE (mode exploration)"
        best = cands[0]
        gain = float((obs.res ** 2).sum())
        obs.res = np.clip(obs.res - best["ycalc"], 0, None)
        best["share"] = (gain - float((obs.res ** 2).sum())) / obs.ss0
        best["area"] = float(best["ycalc"].sum())
        found.append(best)
        cl = best["cell"]
        log("  -> retenue : %s (%s) — confiance %s" % (best["row"]["ident"], best["row"]["name"], best["conf"]))
        log("     maille ajustée : %.4f %.4f %.4f  %.3f %.3f %.3f  (décalage isotrope initial %+.1f %%)" % (*cl, best["shift"]))
    log("\n=== Bilan (%.0f s) ===" % (time.time() - t00))
    tot = sum(f["area"] for f in found) or 1.0
    for i, f in enumerate(found):
        log("  phase %d : %-26s %-34s part de l'intensité calculée %4.0f %%   contraste %.1f (confiance %s)" % (i + 1, f["row"]["ident"][-26:], (f["row"]["name"] or "")[:34], 100 * f["area"] / tot, f["contrast"], f["conf"]))
    log("  signal non expliqué : %.1f %% (en intensité au carré)" % (100 * float((obs.res ** 2).sum()) / obs.ss0))
    if out:
        os.makedirs(out, exist_ok=True)
        with open(os.path.join(out, "pm_results.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["tour", "rang", "source", "identifiant", "nom", "formule", "groupe", "T_K", "score_grossier", "rang_grossier", "score_fin",
                        "score_valide", "famille", "a", "b", "c", "alpha", "beta", "gamma", "meme_phase_que_n1"])
            w.writerows(table)
        for i, f in enumerate(found):
            name = "phase%d_%s.cif" % (i + 1, "".join(ch if ch.isalnum() else "_" for ch in os.path.basename(f["row"]["ident"])))
            with open(os.path.join(out, name), "w") as fh:
                fh.write(f["cif"])
        log("  résultats : %s (tableau CSV + CIF des phases retenues)" % out)
    return found, obs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gpx")
    ap.add_argument("--db", default=pc.default_db())
    ap.add_argument("--hist", type=int)
    ap.add_argument("--max-phases", type=int, default=4)
    ap.add_argument("--sources", help="sources à interroger, séparées par des virgules (défaut : toutes)")
    ap.add_argument("--elements", help='éléments exigés, ex. "C H N O" (facultatif)')
    ap.add_argument("--strict", action="store_true", help="avec --elements : aucun autre élément (H toléré)")
    ap.add_argument("--tol", type=float, default=0.4, help="tolérance du tri grossier sur d, en %% (défaut 0,4)")
    ap.add_argument("--min-gain", type=float, default=0.12, help="part minimale du résidu qu'une phase doit expliquer (défaut 0,12)")
    ap.add_argument("--min-contrast", type=float, default=1.5, help="rapport minimal entre la meilleure référence et les autres phases (défaut 1,5 ; confiance forte à partir de 2)")
    ap.add_argument("--explore", action="store_true", help="exploration : à chaque tour, soustraire le meilleur candidat même s'il n'est pas validé")
    ap.add_argument("--n3", type=int, default=6, help="nombre de candidats validés avec la structure complète à chaque tour (défaut 6)")
    ap.add_argument("--out")
    a = ap.parse_args()
    out = a.out or os.path.splitext(a.gpx)[0] + "_pm"
    search(a.gpx, a.db, hist=a.hist, max_phases=a.max_phases, sources=a.sources.split(",") if a.sources else None,
           elements=a.elements.split() if a.elements else None, strict=a.strict, tol=a.tol, min_gain=a.min_gain, min_contrast=a.min_contrast, out=out, explore=a.explore, n3=a.n3)


if __name__ == "__main__":
    main()
