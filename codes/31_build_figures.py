import csv
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.lines import Line2D
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figfont

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
DATA = os.path.join(ROOT, "data")
FIGS = os.path.join(ROOT, "figures")
os.makedirs(FIGS, exist_ok=True)

S1, S2, S3, S4 = "#146b4f", "#b0257a", "#8a6d0f", "#5c3a7d"
GREY = "#b9b8b3"
WHITE = "#ffffff"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
TITLE_SIZE = 9.5
TITLE_PAD = 7
DPI = 400

figfont.use_arial()
plt.rcParams.update({
    "figure.facecolor": WHITE, "axes.facecolor": WHITE, "savefig.facecolor": WHITE,
    "font.family": "Arial", "font.size": 8.5,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.7,
    "axes.labelcolor": INK, "axes.titlesize": TITLE_SIZE,
    "axes.titleweight": "bold",
    "xtick.color": INK2, "ytick.color": INK2,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.frameon": False, "legend.fontsize": 8,
    "grid.color": "#e9e8e4", "grid.linewidth": 0.5,
    "lines.linewidth": 1.8, "lines.markersize": 4.5,
})


def die(msg):
    sys.exit("STOP: %s" % msg)


def jload(*parts):
    p = os.path.join(DATA, *parts)
    if not os.path.isfile(p):
        die("missing %s" % p)
    with open(p) as fh:
        return json.load(fh)


def tload(*parts):
    p = os.path.join(DATA, *parts)
    if not os.path.isfile(p):
        die("missing %s" % p)
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def tidy(ax, ylab=None, xlab=None, title=None, grid="y"):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    if grid:
        ax.set_axisbelow(True)
        ax.grid(axis=grid, which="major")
    if ylab:
        ax.set_ylabel(ylab)
    if xlab:
        ax.set_xlabel(xlab)
    if title:
        ax.set_title(title, loc="left", color=INK, pad=TITLE_PAD)


def save(fig, name):
    p = os.path.join(FIGS, name + ".png")
    fig.savefig(p, dpi=DPI, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(os.path.join(FIGS, name + ".svg"), bbox_inches="tight",
                pad_inches=0.04)
    plt.close(fig)
    print("  %-14s %7.0f kB" % (name, os.path.getsize(p) / 1024.0))
    return p


N = jload("derived_values.json")
PSET = tload("pharmacophore_position_set.tsv")
THR = N["sweep"]["threshold"]


def figure1():
    fig = plt.figure(figsize=(7.0, 5.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.72, 1.0], hspace=0.40, wspace=0.34)

    ax = fig.add_subplot(gs[0, :])
    ax.set_xlim(0, 10.2); ax.set_ylim(0, 3.0); ax.axis("off")
    ax.set_title("a",
                 loc="left", color=INK, pad=2)
    boxes = [(0.10, "Human PCFT\nSLC46A1, Q96NT5",
              "sequence only,\nno experimental structure", GREY),
             (3.75, "Chicken PCFT\n7BC6 and 7BC7",
              "cryo-EM coordinates,\n%.1f%% identical" % N["map"]["identity_pct"], S1),
             (7.40, "Humanised surrogate",
              "%d positions changed\non the chicken scaffold"
              % N["construct"]["n_positions_changed"], S3)]
    for x, t, sub, c in boxes:
        ax.add_patch(FancyBboxPatch((x, 0.70), 2.70, 1.25, boxstyle="round,pad=0.05",
                                    linewidth=1.1, edgecolor=c, facecolor=WHITE))
        ax.text(x + 1.35, 1.63, t, ha="center", va="center", fontsize=8.4,
                weight="bold", color=INK)
        ax.text(x + 1.35, 1.08, sub, ha="center", va="center", fontsize=7.4,
                color=INK2, linespacing=1.4)
    for x0, x1, lab in ((2.80, 3.75, "read the set onto\nthe scaffold"),
                        (6.45, 7.40, "substitute, then\nrescue what breaks")):
        ax.add_patch(FancyArrowPatch((x0, 1.33), (x1, 1.33), arrowstyle="-|>",
                                     mutation_scale=10, linewidth=1.0, color=INK2))
        ax.text((x0 + x1) / 2.0, 2.38, lab, ha="center", va="center",
                fontsize=7.0, color=INK2, linespacing=1.4)
    ax.text(5.1, 0.14,
            "%d pharmacophore positions, %d of them different between the two species"
            % (N["set"]["n_positions"], N["set"]["n_primaries"]),
            ha="center", va="center", fontsize=8.0, color=INK)

    ax = fig.add_subplot(gs[1, 0])
    rows = sorted(PSET, key=lambda r: int(r["chicken_pos"]))
    cp = np.array([int(r["chicken_pos"]) for r in rows])
    off = np.array([int(r["chicken_pos"]) - int(r["human_pos"]) for r in rows])
    named = np.array([r["in_parker_named"] == "1" for r in rows])
    ax.scatter(cp[~named], off[~named], s=15, color=GREY, zorder=3,
               label="position set")
    ax.scatter(cp[named], off[named], s=26, color=S2, zorder=4,
               label="stated by Parker et al.")
    ax.set_yticks(sorted(set(off.tolist())))
    ax.set_ylim(min(off) - 2.2, max(off) + 2.6)
    tidy(ax, ylab="chicken minus human position", xlab="chicken position",
         title="b")
    ax.legend(loc="center left", handletextpad=0.3, borderpad=0.1)

    ax = fig.add_subplot(gs[1, 1])
    have = [(int(r["chicken_pos"]), float(r["min_pemetrexed_dist_A"]),
             r["differs"] == "1")
            for r in PSET if r["min_pemetrexed_dist_A"] not in ("", None)]
    have.sort(key=lambda t: t[1])
    d = np.array([h[1] for h in have])
    prim = np.array([h[2] for h in have])
    idx = np.arange(len(have))
    ax.barh(idx[~prim], d[~prim], color=GREY, height=0.74, label="conserved")
    ax.barh(idx[prim], d[prim], color=S2, height=0.74, label="humanising")
    ax.set_yticks(idx)
    ax.set_yticklabels([str(h[0]) for h in have], fontsize=5.8)
    ax.set_ylim(len(have) - 0.3, -2.6)
    tidy(ax, xlab="minimum distance to pemetrexed, " + u"\u00c5",
         ylab="chicken position", grid="x",
         title="c")
    ax.legend(loc="upper right", handletextpad=0.3, borderpad=0.1, ncol=2,
              columnspacing=0.8)
    return save(fig, "fig1_engraftment_and_map")


def figure2():
    fig = plt.figure(figsize=(7.0, 2.6))
    gs = fig.add_gridspec(1, 3, wspace=0.42, width_ratios=[1.0, 1.0, 1.15])
    ladder = N["coupled"]["ladder"]
    cuts = sorted(float(k) for k in ladder)

    ax = fig.add_subplot(gs[0, 0])
    cp = [ladder["%.1f" % c]["coupled"] for c in cuts]
    pr = [ladder["%.1f" % c]["pairs"] for c in cuts]
    ax.plot(cuts, cp, "o-", color=S1, label="positions, of %d"
            % N["coupled"]["n_primary"])
    ax.plot(cuts, pr, "s-", color=S2, label="pairs")
    for x, y in zip(cuts, cp):
        ax.annotate(str(y), (x, y), textcoords="offset points", xytext=(0, 6),
                    ha="center", fontsize=7, color=S1)
    tidy(ax, ylab="count", xlab="contact cutoff, " + u"Å",
         title="a")
    ax.set_ylim(0, max(pr) * 1.32)
    ax.legend(loc="upper left", handletextpad=0.3, borderpad=0.1)

    ax = fig.add_subplot(gs[0, 1])
    cut = "%.1f" % N["coupled"]["reported_cutoff_A"]
    edges = [(e[0], e[1], e[2]) for e in N["coupled"]["cutoffs"][cut]["edges"]]
    sep = sorted(e[2] for e in edges)
    ax.bar(range(len(sep)), sep, color=[S3 if s > N["coupled"]["window"] else GREY
                                        for s in sep], width=0.78)
    ax.axhline(N["coupled"]["window"], color=S2, linewidth=1.0, linestyle="--")
    ax.text(0.3, N["coupled"]["window"] * 1.5,
            "published search window, " + u"±%d" % N["coupled"]["window"],
            fontsize=7.2, color=S2)
    ax.set_yscale("log")
    tidy(ax, ylab="sequence separation, residues",
         xlab="coupled pair, ordered by separation",
         title="b")

    ax = fig.add_subplot(gs[0, 2])
    comps = N["coupled"]["cutoffs"][cut]["components"]
    iso = N["coupled"]["cutoffs"][cut].get("isolated_positions") or []
    pos = {}
    for ci, comp in enumerate(comps):
        n = len(comp)
        for k, p in enumerate(sorted(comp)):
            ang = 2 * np.pi * k / max(n, 1)
            pos[p] = (ci * 2.4 + 0.95 * np.cos(ang), 0.95 * np.sin(ang))
    for k, p in enumerate(sorted(iso)):
        pos[p] = (len(comps) * 2.4 + 0.45, 1.0 - 0.9 * k)
    for a, b, _ in edges:
        if a in pos and b in pos:
            ax.plot([pos[a][0], pos[b][0]], [pos[a][1], pos[b][1]],
                    color=GREY, linewidth=0.9, zorder=1)
    for p, (x, y) in pos.items():
        ax.scatter([x], [y], s=120, color=WHITE, edgecolor=S1 if p not in iso else GREY,
                   linewidth=1.1, zorder=3)
        ax.text(x, y, str(p), ha="center", va="center", fontsize=5.8, color=INK,
                zorder=4)
    ax.axis("off")
    ax.set_title("c", loc="left", color=INK, pad=TITLE_PAD)
    return save(fig, "fig2_coupled_set")


def figure3():
    fig = plt.figure(figsize=(7.0, 2.6))
    gs = fig.add_gridspec(1, 3, wspace=0.40)
    sweep = jload("rescue_sweep_report.json")
    pub = [r for r in sweep["summary"] if r["policy"] == "published"]

    ax = fig.add_subplot(gs[0, 0])
    lab = [str(r["radius"]) for r in pub]
    cross = [r["rescuers_crossing"] for r in pub]
    resc = [r["targets_rescued"] for r in pub]
    x = np.arange(len(pub))
    ax.bar(x, cross, color=S1, width=0.68)
    for xi, c in zip(x, cross):
        ax.annotate(str(c), (xi, c), textcoords="offset points", xytext=(0, 3),
                    ha="center", fontsize=7, color=INK2)
    ax.set_xticks(x); ax.set_xticklabels(lab)
    ax.set_ylim(0, max(cross) * 1.22)
    tidy(ax, ylab="substitutions crossing the threshold",
         xlab="search radius, residues",
         title="a")
    ax2 = ax.twinx()
    ax2.plot(x, resc, "o-", color=S2)
    ax2.set_ylim(0, max(resc) + 2)
    ax2.set_ylabel("targets rescued", color=S2)
    ax2.tick_params(axis="y", colors=S2)
    for side in ("top",):
        ax2.spines[side].set_visible(False)

    ax = fig.add_subplot(gs[0, 1])
    items = [("chicken %d" % N["pairs_407"]["target"], N["pairs_407"]["baseline"],
              N["pairs_407"]["best_single_p"],
              N["pairs_407"]["best_double"]["p_target"],
              N["pairs_407"]["best_double"]),
             ("chicken %d" % N["pairs_144"]["target"], N["pairs_144"]["baseline"],
              N["pairs_144"]["best_double"]["best_single_of_pair"],
              N["pairs_144"]["best_double"]["p_target"],
              N["pairs_144"]["best_double"])]
    w = 0.26
    for i, (name, base, single, double, bd) in enumerate(items):
        ax.bar(i - w, base, width=w, color=GREY)
        ax.bar(i, single, width=w, color=S1)
        ax.bar(i + w, double, width=w, color=S3)
        ax.annotate("%s%d%s with %s%d%s" % (bd["from_i"], bd["pos_i"], bd["to_i"],
                                            bd["from_j"], bd["pos_j"], bd["to_j"]),
                    (i + w, double), textcoords="offset points", xytext=(0, 4),
                    ha="center", fontsize=6.4, color=S3)
    ax.axhline(THR, color=S2, linewidth=1.0, linestyle="--")
    ax.text(-0.42, THR * 1.25, "threshold %g" % THR, fontsize=7.2, color=S2)
    ax.set_yscale("log")
    ax.set_xticks(range(len(items)))
    ax.set_xticklabels([i[0] for i in items])
    tidy(ax, ylab="sequence likelihood of the human residue",
         title="b")
    ax.legend(handles=[Line2D([], [], color=GREY, linewidth=5, label="unrescued"),
                       Line2D([], [], color=S1, linewidth=5, label="best single"),
                       Line2D([], [], color=S3, linewidth=5, label="best pair")],
              loc="upper left", handletextpad=0.4, borderpad=0.1)

    ax = fig.add_subplot(gs[0, 2])
    dbl = tload("esm2_step8", "c407_double_rescue.tsv")
    fold = np.array([float(r["fold_vs_best_single_of_pair"]) for r in dbl])
    ax.hist(fold, bins=np.logspace(np.log10(max(fold.min(), 1e-3)),
                                   np.log10(fold.max()), 44),
            color=S1, edgecolor=WHITE, linewidth=0.3)
    ax.axvline(1.0, color=INK2, linewidth=0.9, linestyle=":")
    ax.axvline(np.median(fold), color=S2, linewidth=1.1)
    ax.annotate("median %.2f" % np.median(fold), (np.median(fold), 0.86),
                xycoords=("data", "axes fraction"), xytext=(5, 0), fontsize=7.2,
                textcoords="offset points", color=S2)
    ax.set_xscale("log"); ax.set_yscale("log")
    tidy(ax, ylab="pairs", xlab="pair over the better of its two halves",
         title="c")
    return save(fig, "fig3_rescue_search")


def figure4():
    fig = plt.figure(figsize=(7.0, 2.9))
    gs = fig.add_gridspec(1, 3, wspace=0.56)

    ax = fig.add_subplot(gs[0, 0])
    ts = sorted(N["set_effect"]["targets"], key=lambda t: t["target"])
    x = np.arange(len(ts))
    alone = [t["alone"] for t in ts]
    inset = [t["primaries_rescuers_primary_wins"] for t in ts]
    for xi, a, b in zip(x, alone, inset):
        ax.plot([xi, xi], [a, b], color=GREY, linewidth=1.2, zorder=1)
    ax.scatter(x, alone, s=34, color=S1, zorder=3, label="measured alone")
    ax.scatter(x, inset, s=34, color=S2, zorder=3, label="with the set present")
    ax.axhline(THR, color=INK2, linewidth=0.9, linestyle="--")
    ax.text(-0.45, THR * 1.22, "threshold", fontsize=7.2, color=INK2)
    ax.set_yscale("log")
    ax.set_xticks(x); ax.set_xticklabels([str(t["target"]) for t in ts])
    tidy(ax, ylab="sequence likelihood", xlab="chicken position",
         title="a")
    ax.legend(loc="lower left", handletextpad=0.3, borderpad=0.1)

    ax = fig.add_subplot(gs[0, 1])
    scan = tload("esm2_step6", "construct_rescue_scan.tsv")
    f = np.array([float(r["fold_construct_over_wt"]) for r in scan
                  if r["fold_construct_over_wt"] not in ("", "nan")])
    ax.hist(f, bins=np.logspace(np.log10(max(f.min(), 1e-3)), np.log10(f.max()), 50),
            color=S1, edgecolor=WHITE, linewidth=0.3)
    ax.axvline(1.0, color=INK2, linewidth=0.9, linestyle=":")
    ax.axvline(np.median(f), color=S2, linewidth=1.1)
    ax.annotate("median %.2f\n%s substitutions" % (np.median(f), f"{len(f):,}"),
                (0.03, 0.83), xycoords="axes fraction", fontsize=7.2, color=S2)
    ax.set_xscale("log")
    tidy(ax, ylab="substitutions", xlab="construct over wild type",
         title="b")

    ax = fig.add_subplot(gs[0, 2])
    sing = tload("esm2_step8", "c407_singles.tsv")
    paired = [r for r in sing
              if r["p_wt_step2"] not in ("", "nan", None)
              and float(r["p_wt_step2"]) > 0]
    g = np.array([float(r["p_construct"]) / float(r["p_wt_step2"]) for r in paired])
    ax.hist(g, bins=36, color=S3, edgecolor=WHITE, linewidth=0.3)
    ax.axvline(1.0, color=INK2, linewidth=0.9, linestyle=":")
    ax.axvline(np.median(g), color=S2, linewidth=1.1)
    ax.annotate("median %.2f\n%d candidates" % (np.median(g), len(g)),
                (0.03, 0.83), xycoords="axes fraction", fontsize=7.2, color=S2)
    tidy(ax, ylab="substitutions", xlab="construct over wild type",
         title="c")
    return save(fig, "fig4_set_effect")


def figure5():
    fig = plt.figure(figsize=(7.0, 5.2))
    gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.34)

    prof = tload("esm2_step9", "revised_profile.tsv")
    prof.sort(key=lambda r: float(r["p_revised"]))
    p = np.array([float(r["p_revised"]) for r in prof])
    labels = [r["chicken_pos"] for r in prof]
    below = p < THR

    ax = fig.add_subplot(gs[0, :])
    x = np.arange(len(p))
    ax.bar(x[~below], p[~below], color=S1, width=0.72)
    ax.bar(x[below], p[below], color=S2, width=0.72)
    ax.axhline(THR, color=INK2, linewidth=1.0, linestyle="--")
    ax.text(len(p) - 0.4, THR * 1.18, "threshold %g" % THR, fontsize=7.2,
            color=INK2, ha="right")
    for xi in x[below]:
        ax.annotate(labels[xi], (xi, p[xi]), textcoords="offset points",
                    xytext=(0, 4), ha="center", fontsize=7, color=S2)
    ax.set_yscale("log")
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=6.4, rotation=90)
    tidy(ax, ylab="likelihood of the human residue", xlab="chicken position",
         title="a")

    ax = fig.add_subplot(gs[1, 0])
    dropped = N["step9"]["neighbourhood_of_the_swap"]["dropped"][0]
    near = set(N["step9"]["neighbourhood_of_the_swap"]["near"])
    d = np.array([abs(int(r["chicken_pos"]) - dropped) for r in prof])
    fo = np.array([float(r["fold_revised_over_step7"]) for r in prof])
    isnear = np.array([int(r["chicken_pos"]) in near for r in prof])
    ax.scatter(d[~isnear], fo[~isnear], s=22, color=GREY, zorder=3,
               label="elsewhere")
    ax.scatter(d[isnear], fo[isnear], s=34, color=S3, zorder=4,
               label="within %d residues" % N["step9"]["neighbourhood_of_the_swap"]["near_max_distance"])
    ax.axhline(1.0, color=INK2, linewidth=0.9, linestyle=":")
    ax.set_xscale("log")
    ax.set_xlim(0.72, max(d) * 1.7)
    ax.set_ylim(0.6, max(fo) * 1.34)
    tidy(ax, ylab="fold over the previous construct",
         xlab="distance from the dropped rescuer at %d, residues" % dropped,
         title="b")
    ax.legend(loc="upper right", handletextpad=0.3, borderpad=0.1)

    ax = fig.add_subplot(gs[1, 1])
    rng = N["af3"]["plddt_range_at_watch"]
    arms = ["wt", "primaries", "construct"]
    cols = {"wt": GREY, "primaries": S1, "construct": S2}
    watch = [str(w) for w in N["af3"]["watched_positions"]]
    for j, q in enumerate(watch):
        for k, a in enumerate(arms):
            lo, hi = rng[a][q]
            xp = j + (k - 1) * 0.24
            ax.plot([xp, xp], [lo, hi], color=cols[a], linewidth=5,
                    solid_capstyle="butt",
                    label=a if j == 0 else None)
    lo = min(v[0] for a in arms for v in rng[a].values())
    hi = max(v[1] for a in arms for v in rng[a].values())
    ax.set_ylim(lo - 0.22 * (hi - lo), hi + 0.06 * (hi - lo))
    ax.set_xlim(-0.55, len(watch) - 0.45)
    ax.set_xticks(range(len(watch)))
    ax.set_xticklabels(["chicken " + q for q in watch])
    tidy(ax, ylab="pLDDT over %d models" % N["af3"]["models_per_arm"]["wt"],
         title="c")
    ax.legend(loc="lower center", handletextpad=0.4, borderpad=0.1, ncol=3,
              columnspacing=0.9)
    return save(fig, "fig5_construct_and_models")


def figure7():
    arms = ("wt", "primaries", "construct")
    cols = {a: tload("dynamics_aware", "%s.tsv" % a) for a in arms}
    dyn = {a: [r for r in cols[a] if r["dcs_0_10"]] for a in arms}
    if not all(dyn[a] for a in arms):
        die("no scored positions in the dynamics-aware tables")
    DN = N["dynamics"]
    watch = [str(w) for w in DN["watched"]]
    targets = sorted({str(k) for k in DN["arm"]["wt"]["watched"]}, key=int)

    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.2),
                             gridspec_kw={"width_ratios": [1.3, 1.0],
                                          "wspace": 0.3})

    ax = axes[0]
    rows = dyn["wt"]
    nd = np.array([float(r["norm_dynamics"]) for r in rows])
    cn = np.array([float(r["conservation"]) for r in rows])
    pos = np.array([int(r["chicken_pos"]) for r in rows])
    inset = np.array([r["in_position_set"] == "1" for r in rows])
    resc = np.array([r["is_rescuer"] == "1" for r in rows])
    raw = cn * nd
    k = int(DN["top_n"])
    bound = np.sort(raw)[-k]
    ax.scatter(nd[~inset & ~resc], cn[~inset & ~resc], s=9, color=GREY,
               zorder=2, linewidths=0, label="other scored positions")
    ax.scatter(nd[inset], cn[inset], s=22, color=S1, zorder=4, linewidths=0,
               label="pharmacophore set")
    ax.scatter(nd[resc], cn[resc], s=30, color=S4, zorder=5, marker="^",
               linewidths=0, label="rescuer")
    xs = np.linspace(max(nd.min(), 1e-3), nd.max(), 300)
    with np.errstate(divide="ignore", invalid="ignore"):
        ax.plot(xs, bound / xs, color=INK2, linewidth=0.8, linestyle=(0, (4, 2)),
                zorder=3, label="top %s per cent boundary"
                % ("%.0f" % (100 * float(DN["top_fraction"]))))
    for w in watch:
        hit = [i for i, p in enumerate(pos) if p == int(w)]
        if not hit or raw[hit[0]] < bound:
            continue
        i = hit[0]
        ax.annotate("chicken %s" % w, (nd[i], cn[i]),
                    textcoords="offset points", xytext=(8, -1),
                    fontsize=7.0, color=INK, va="center")
    ax.set_ylim(0, max(cn.max() * 1.12, 0.1))
    ax.set_xlim(0, nd.max() * 1.05)
    tidy(ax, ylab="column conservation", xlab="normalised fluctuation",
         title="a",
         grid="both")
    ax.legend(loc="upper right", handletextpad=0.4, borderpad=0.25,
              labelspacing=0.3, fontsize=6.8)

    ax = axes[1]
    byarm = {a: {int(r["chicken_pos"]): r for r in dyn[a]} for a in arms}
    hues = [S2, S1, S3, S4]
    x = np.arange(len(arms))
    for j, w in enumerate(targets):
        p = int(w)
        if not all(p in byarm[a] for a in arms):
            continue
        y = [float(byarm[a][p]["dcs_0_10"]) for a in arms]
        c = hues[j % len(hues)]
        ax.plot(x, y, color=c, linewidth=1.3, marker="o", markersize=4,
                markeredgewidth=0, zorder=4)
        ax.annotate(" %s" % w, (x[-1], y[-1]), textcoords="offset points",
                    xytext=(4, -1), fontsize=7.0, color=c, va="center")
    ax.set_xticks(x)
    ax.set_xticklabels(["wild type", "primaries", "construct"])
    ax.set_xlim(-0.25, len(arms) - 0.45)
    tidy(ax, ylab="dynamic conserved score",
         title="b")
    return save(fig, "fig7_dynamics_aware")


if __name__ == "__main__":
    print("figures from %s" % os.path.relpath(DATA, ROOT))
    built = [figure1(), figure2(), figure3(), figure4(), figure5(),
             figure7()]
    missing = [b for b in built if not os.path.isfile(b)]
    if missing:
        die("these did not write: %s" % missing)
    print("%d figures in %s" % (len(built), os.path.relpath(FIGS, ROOT)))
