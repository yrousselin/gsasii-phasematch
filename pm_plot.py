"""pm_plot : tracé de la session de recherche de phases dans une figure matplotlib (sans dépendance à wx)."""
import numpy as np

CAND_COLOR = "#d62728"


def draw(fig, session, cand=None, keep_xlim=None):
    """Haut : diagramme mesuré (fond soustrait), phases retenues et candidat sélectionné.
    Milieu : positions des raies. Bas : ce qui reste à expliquer, avec le candidat."""
    fig.clear()
    obs = session.obs
    gs = fig.add_gridspec(3, 1, height_ratios=[5, 0.9, 2.2], hspace=0.05, left=0.07, right=0.985, top=0.95, bottom=0.09)
    ax = fig.add_subplot(gs[0])
    at = fig.add_subplot(gs[1], sharex=ax)
    ar = fig.add_subplot(gs[2], sharex=ax)
    x = obs.x
    ax.plot(x, session.net, color="black", lw=0.8, label="mesuré (fond soustrait)", zorder=5)
    tot = np.zeros_like(x)
    for a in session.accepted:
        ax.fill_between(x, a["ycalc"], color=a["color"], alpha=0.35, lw=0)
        ax.plot(x, a["ycalc"], color=a["color"], lw=0.9, label="%s  %s" % (a["ident"], a["name"][:28]))
        tot += a["ycalc"]
    rows = list(session.accepted)
    if cand is not None:
        ax.plot(x, cand["ycalc"], color=CAND_COLOR, lw=1.1, label="candidat : %s  %s" % (cand["ident"], cand["name"][:28]), zorder=6)
        rows = rows + [dict(cand, color=CAND_COLOR)]
    # positions des raies, une ligne par phase
    n = max(1, len(rows))
    for i, a in enumerate(rows):
        t, h = a["sticks"]
        if len(t):
            keep = h >= 0.02 * h.max()
            at.vlines(t[keep], n - i - 0.9, n - i - 0.1, color=a["color"], lw=0.8)
    at.set_ylim(0, n)
    at.set_yticks([])
    at.tick_params(labelbottom=False)
    ax.tick_params(labelbottom=False)
    # reste à expliquer
    ar.plot(x, obs.res, color="0.35", lw=0.8, label="reste à expliquer")
    if cand is not None:
        ar.plot(x, cand["ycalc"], color=CAND_COLOR, lw=1.0, alpha=0.9)
    ar.set_xlabel("2θ (°)")
    ar.legend(loc="upper right", fontsize=8, frameon=False)
    ax.set_ylabel("Intensité")
    ax.legend(loc="upper right", fontsize=8, frameon=False)
    ymax = float(session.net.max()) or 1.0
    ax.set_ylim(-0.02 * ymax, 1.05 * ymax)
    rmax = max(float(obs.res.max()), float(cand["ycalc"].max()) if cand is not None else 0.0) or 1.0
    ar.set_ylim(-0.03 * rmax, 1.08 * rmax)
    ax.set_xlim(*(keep_xlim or (x[0], x[-1])))
    ttl = "Tour %d — %d phase(s) retenue(s) ; reste à expliquer : %.0f %% de l'intensité" % (session.round, len(session.accepted), 100 * session.unexplained())
    ax.set_title(ttl, fontsize=10)
    return ax, at, ar
