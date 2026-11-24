import csv
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figfont

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
_sf = __import__("33_structure_figures")
trace, kabsch, core = _sf.trace, _sf.kabsch, _sf.core
S1, S2, S3, S4 = _sf.S1, _sf.S2, _sf.S3, _sf.S4
GREY, WHITE, INK, INK2, MUTED = _sf.GREY, _sf.WHITE, _sf.INK, _sf.INK2, _sf.MUTED

ROOT = os.path.join(HERE, "..")
DATA = os.path.join(ROOT, "data")
FIGS = os.path.join(ROOT, "figures")
DPI = 400

figfont.use_arial()
plt.rcParams.update({
    "figure.facecolor": WHITE, "axes.facecolor": WHITE, "savefig.facecolor": WHITE,
    "font.family": "Arial", "font.size": 8.5,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.7, "axes.labelcolor": INK,
    "axes.titlesize": 9.5, "axes.titleweight": "bold",
    "xtick.color": INK2, "ytick.color": INK2,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.frameon": False, "legend.fontsize": 8,
    "grid.color": "#e9e8e4", "grid.linewidth": 0.5, "lines.linewidth": 1.4,
})


def die(msg):
    sys.exit("STOP: %s" % msg)


def jload(*p):
    f = os.path.join(DATA, *p)
    if not os.path.isfile(f):
        die("missing %s" % f)
    return json.load(open(f))


def tload(*p):
    f = os.path.join(DATA, *p)
    if not os.path.isfile(f):
        die("missing %s" % f)
    return list(csv.DictReader(open(f), delimiter="\t"))


def tidy(ax, ylab=None, xlab=None, title=None, grid="y"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if grid:
        ax.set_axisbelow(True); ax.grid(axis=grid)
    if ylab:
        ax.set_ylabel(ylab)
    if xlab:
        ax.set_xlabel(xlab)
    if title:
        ax.set_title(title, loc="left", color=INK, pad=6)


def save(fig, name):
    p = os.path.join(FIGS, name + ".png")
    fig.savefig(p, dpi=DPI, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(os.path.join(FIGS, name + ".svg"), bbox_inches="tight",
                pad_inches=0.04)
    plt.close(fig)
    print("  %-16s %7.0f kB" % (name, os.path.getsize(p) / 1024.0))


N = jload("derived_values.json")
ARMS = ("wt", "primaries", "construct")
COLS = {"wt": GREY, "primaries": S3, "construct": S1}
NAME = {"wt": "wild type", "primaries": "humanising only", "construct": "construct"}
MODEL = {a: os.path.join(DATA, "af3", a, "fold_%s_model_0.cif" % a) for a in ARMS}
PSETROWS = tload("pharmacophore_position_set.tsv")
PSET = sorted(int(r["chicken_pos"]) for r in PSETROWS)
PRIMARIES = sorted(int(r["chicken_pos"]) for r in PSETROWS if r["differs"] == "1")
RESCUERS = sorted(int(k) for k in N["construct"]["rescuer_positions"])
D = {}


def figure_s1():
    fig, ax = plt.subplots(figsize=(7.0, 2.5))
    for a in ARMS:
        n, x, b, _ = trace(MODEL[a])
        n, x, b = core(n, x, b)
        ax.plot(n, b, color=COLS[a], label=NAME[a], linewidth=1.1)
    for p in PRIMARIES:
        ax.axvline(p, color=S2, linewidth=0.4, alpha=0.30, zorder=0)
    ax.set_ylim(40, 100)
    tidy(ax, ylab="pLDDT", xlab="chicken position")
    ax.legend(loc="lower center", ncol=3, columnspacing=1.2)
    save(fig, "figS1_plddt_profile")


def figure_s2():
    ref_n, ref_x, _, _ = trace(MODEL["wt"])
    ref_n, ref_x = core(ref_n, ref_x)
    fig, ax = plt.subplots(figsize=(7.0, 2.5))
    out = {}
    for a in ("primaries", "construct"):
        n, x, _, _ = trace(MODEL[a])
        n, x = core(n, x)
        shared = sorted(set(n.tolist()) & set(ref_n.tolist()))
        ia = {int(k): i for i, k in enumerate(n)}
        ir = {int(k): i for i, k in enumerate(ref_n)}
        P = x[[ia[s] for s in shared]]
        Q = ref_x[[ir[s] for s in shared]]
        R, t = kabsch(P, Q)
        d = np.linalg.norm((P @ R.T + t) - Q, axis=1)
        ax.plot(shared, d, color=COLS[a], label=NAME[a] + " against wild type",
                linewidth=1.0)
        out[a] = {"median": float(np.median(d)), "p95": float(np.percentile(d, 95)),
                  "max": float(d.max()),
                  "max_at": int(shared[int(np.argmax(d))]), "n": len(d)}
    for p in PRIMARIES:
        ax.axvline(p, color=S2, linewidth=0.4, alpha=0.30, zorder=0)
    tidy(ax, ylab="deviation, " + u"Å", xlab="chicken position")
    ax.legend(loc="upper left", ncol=2, columnspacing=1.2)
    save(fig, "figS2_per_residue_deviation")
    D["deviation"] = out


def figure_s3():
    sweep = jload("rescue_sweep_report.json")["summary"]
    pols = sorted({r["policy"] for r in sweep})
    radii = [r["radius"] for r in sweep if r["policy"] == pols[0]]
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.5))
    style = {"published": S1, "strict": S2, "annotated": S4}
    for pol in pols:
        rows = [r for r in sweep if r["policy"] == pol]
        axes[0].plot(range(len(rows)), [r["rescuers_crossing"] for r in rows],
                     "o-", color=style[pol], label=pol, markersize=3.6)
        axes[1].plot(range(len(rows)), [r["targets_rescued"] for r in rows],
                     "o-", color=style[pol], label=pol, markersize=3.6)
    for ax, ylab, ttl in ((axes[0], "substitutions crossing",
                           "a"),
                          (axes[1], "targets rescued", "b")):
        ax.set_xticks(range(len(radii)))
        ax.set_xticklabels([str(r) for r in radii])
        tidy(ax, ylab=ylab, xlab="search radius, residues", title=ttl)
    axes[1].set_ylim(0, max(r["targets_rescued"] for r in sweep) + 1)
    axes[0].legend(loc="upper left")
    fig.suptitle("Figure S3   The window and policy sweep in full", x=0.0,
                 ha="left", fontsize=9.5, fontweight="bold", y=1.04)
    save(fig, "figS3_policy_sweep")
    D["sweep_rows"] = [{k: r[k] for k in ("radius", "policy", "targets_rescued",
                                          "rescuers_crossing",
                                          "nearest_rescuer_separation")}
                       for r in sweep]


figure_s1()
figure_s2()
figure_s3()

D["set_rows"] = [{"chicken_pos": int(r["chicken_pos"]),
                  "human_pos": int(r["human_pos"]),
                  "chicken_aa": r["chicken_aa"], "human_aa": r["human_aa"],
                  "differs": r["differs"] == "1",
                  "parker": r["in_parker_named"] == "1",
                  "pemetrexed": r["min_pemetrexed_dist_A"],
                  "mutation": r["in_mutation_literature"] == "1"}
                 for r in PSETROWS]
D["ladder"] = N["coupled"]["ladder"]
D["survival"] = N["set_effect"]["targets"]
D["per_arm"] = N["af3"]["per_arm"]
_AM = jload("af3_models_report.json")["all_models"]
_KEEP = ("model_index", "plddt_chain", "ptm", "ranking_score", "rmsd_to_7BC6",
         "fraction_disordered")
D["all_models"] = {a: [{k: m[k] for k in _KEEP if k in m} for m in _AM[a]]
                   for a in ARMS}
_have = sorted({k for a in ARMS for m in _AM[a] for k in m})
if "plddt_chain" not in _have:
    die("the per-model report carries no chain confidence, keys are %s" % _have)
D["counts"] = {"notes": 4, "tables": 5, "figures": 3}
DYN = {a: {int(r["chicken_pos"]): r
             for r in tload("dynamics_aware", "%s.tsv" % a)}
       for a in ("wt", "primaries", "construct")}
for a, rows in DYN.items():
    if not rows:
        die("no dynamics-aware table for the %s arm. Run "
            "codes/35_dynamics_aware.py first" % a)
PSET_DYN = [p for p in PSET if p in DYN["wt"]]
if not PSET_DYN:
    die("the dynamics-aware table covers none of the set positions")
D["dynamics_set_absent"] = [p for p in PSET if p not in DYN["wt"]]


def _dyn(a, p, key):
    v = DYN[a][p].get(key, "")
    return float(v) if v not in ("", None) else None


D["dynamics_rows"] = [{"chicken_pos": p,
                       "human_pos": DYN["wt"][p]["human_pos"],
                       "chicken_aa": DYN["wt"][p]["chicken_aa"],
                       "conservation": _dyn("wt", p, "conservation"),
                       "dcs_wt": _dyn("wt", p, "dcs_0_10"),
                       "dcs_construct": _dyn("construct", p, "dcs_0_10"),
                       "rcs_wt": _dyn("wt", p, "rcs_0_10"),
                       "rcs_construct": _dyn("construct", p, "rcs_0_10"),
                       "is_rescuer": DYN["wt"][p]["is_rescuer"] == "1"}
                      for p in PSET_DYN]
D["dynamics_rescuer_rows"] = [{"chicken_pos": p,
                               "human_pos": DYN["wt"][p]["human_pos"],
                               "chicken_aa": DYN["wt"][p]["chicken_aa"],
                               "conservation": _dyn("wt", p, "conservation"),
                               "dcs_wt": _dyn("wt", p, "dcs_0_10"),
                               "dcs_construct": _dyn("construct", p, "dcs_0_10"),
                               "rcs_wt": _dyn("wt", p, "rcs_0_10"),
                               "rcs_construct": _dyn("construct", p, "rcs_0_10")}
                              for p in sorted(DYN["wt"])
                              if DYN["wt"][p]["is_rescuer"] == "1"]
if not D["dynamics_rescuer_rows"]:
    die("no rescuer position appears in the dynamics-aware table")

out = os.path.join(DATA, "supplementary_derived.json")
json.dump(D, open(out, "w"), indent=1, sort_keys=True)
print("wrote %s (%d bytes)" % (os.path.relpath(out, ROOT), os.path.getsize(out)))
