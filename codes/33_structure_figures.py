import json
import os
import sys

import gemmi
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figfont

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
DATA = os.path.join(ROOT, "data")
SEQ = os.path.join(ROOT, "seq")
FIGS = os.path.join(ROOT, "figures")

S1, S2, S3, S4 = "#146b4f", "#b0257a", "#8a6d0f", "#5c3a7d"
GREY, WHITE = "#b9b8b3", "#ffffff"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
SLAB, RULE, LIG = "#f1f0ec", "#d9d8d3", "#3f3e3b"
DPI = 400

HALF_THICKNESS = 15.0
PLDDT_FLOOR = 70.0
LW_HELIX, LW_LOOP, LW_FAINT = 6.4, 2.2, 1.5
SITE_SHELL, SITE_CONTEXT = 16.0, 21.0

KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
      "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
      "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}

figfont.use_arial()
plt.rcParams.update({
    "figure.facecolor": WHITE, "axes.facecolor": WHITE, "savefig.facecolor": WHITE,
    "font.family": "Arial", "font.size": 8.5,
    "axes.titlesize": 9.5, "axes.titleweight": "bold",
    "legend.frameon": False, "legend.fontsize": 7.6,
})


def die(msg):
    sys.exit("STOP: %s" % msg)


def jload(*parts):
    p = os.path.join(DATA, *parts)
    if not os.path.isfile(p):
        die("missing %s" % p)
    with open(p) as fh:
        return json.load(fh)


def rgb(c):
    return np.array(matplotlib.colors.to_rgb(c))


def pale(c, f):
    return 1.0 + (rgb(c) - 1.0) * f


def trace(path, chain=None):
    if not os.path.isfile(path):
        die("missing %s" % path)
    st = gemmi.read_structure(path)
    st.setup_entities()
    best = None
    for ch in st[0]:
        if chain is not None and ch.name != chain:
            continue
        rs = [r for r in ch if r.find_atom("CA", "*") is not None]
        if best is None or len(rs) > len(best):
            best = rs
    if best is None or len(best) < 50:
        die("%s has no chain with a backbone trace" % path)
    nums = np.array([r.seqid.num for r in best])
    xyz = np.array([[r.find_atom("CA", "*").pos.x, r.find_atom("CA", "*").pos.y,
                     r.find_atom("CA", "*").pos.z] for r in best])
    b = np.array([r.find_atom("CA", "*").b_iso for r in best])
    aa = np.array([gemmi.find_tabulated_residue(r.name).one_letter_code.upper()
                   for r in best])
    return nums, xyz, b, aa


def helix_mask(xyz):
    d3 = np.linalg.norm(xyz[3:] - xyz[:-3], axis=1)
    d4 = np.linalg.norm(xyz[4:] - xyz[:-4], axis=1)
    m = np.zeros(len(xyz), bool)
    for i in range(len(d4)):
        if 5.0 < d4[i] < 7.2 and 4.4 < d3[i] < 6.1:
            m[i:i + 5] = True
    return m


def runs_of(mask, minlen=1):
    out, s = [], None
    for i, v in enumerate(mask):
        if v and s is None:
            s = i
        if not v and s is not None:
            if i - s >= minlen:
                out.append((s, i - 1))
            s = None
    if s is not None and len(mask) - s >= minlen:
        out.append((s, len(mask) - 1))
    return out


def axis_path(xyz, mask, win=3):
    out = xyz.copy()
    for s, e in runs_of(mask, minlen=2 * win + 1):
        for i in range(s, e + 1):
            w = min(win, i - s, e - i)
            if w:
                out[i] = xyz[i - w:i + w + 1].mean(axis=0)
    return out


def smooth(xyz, per=6):
    p = np.vstack([xyz[0], xyz, xyz[-1]])
    t = np.linspace(0, 1, per, endpoint=False)[:, None]
    out = []
    for i in range(len(p) - 3):
        p0, p1, p2, p3 = p[i], p[i + 1], p[i + 2], p[i + 3]
        out.append(0.5 * ((2 * p1) + (-p0 + p2) * t
                          + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t ** 2
                          + (-p0 + 3 * p1 - 3 * p2 + p3) * t ** 3))
    out.append(p[-2][None, :])
    return np.vstack(out)


def resample(attr, n_in, n_out, per=6):
    u = np.arange(n_out) / float(per)
    return np.array([np.interp(u, np.arange(n_in), a) for a in np.atleast_2d(attr.T)]).T


def ribbon(ax, xyz, cols, lw, zorder=2, fog=0.46, alpha=1.0):
    if len(xyz) < 2:
        return
    seg = np.stack([xyz[:-1, :2], xyz[1:, :2]], axis=1)
    depth = 0.5 * (xyz[:-1, 2] + xyz[1:, 2])
    w = 0.5 * (np.asarray(lw)[:-1] + np.asarray(lw)[1:])
    c = 0.5 * (np.asarray(cols)[:-1] + np.asarray(cols)[1:])
    o = np.argsort(depth)
    d = depth[o]
    span = max(d.max() - d.min(), 1e-6)
    rel = (d - d.min()) / span
    f = (fog + (1.0 - fog) * rel)[:, None]
    ax.add_collection(LineCollection(
        seg[o], colors=np.clip(1.0 + (c[o] - 1.0) * f, 0, 1),
        linewidths=w[o] * (0.72 + 0.28 * rel), alpha=alpha,
        capstyle="round", joinstyle="round", zorder=zorder))


def draw_chain(ax, xyz, mask, cols, plddt=None, zorder=2, fog=0.46,
               lw_helix=LW_HELIX, lw_loop=LW_LOOP, alpha=1.0):
    if len(xyz) < 4:
        return
    path = axis_path(xyz, mask)
    width = np.where(mask, lw_helix, lw_loop).astype(float)
    col = np.asarray(cols, float)
    if plddt is not None:
        low = plddt < PLDDT_FLOOR
        width = np.where(low, LW_FAINT, width)
        col = np.where(low[:, None], rgb(GREY)[None, :], col)
    sm = smooth(path)
    ribbon(ax, sm, resample(col, len(path), len(sm)),
           resample(width[:, None], len(path), len(sm))[:, 0],
           zorder=zorder, fog=fog, alpha=alpha)


def draw_breaks(ax, nums, xyz, mask, cols, plddt=None, **kw):
    cuts = np.where(np.diff(nums) != 1)[0]
    starts = np.concatenate([[0], cuts + 1])
    ends = np.concatenate([cuts, [len(nums) - 1]])
    for s, e in zip(starts, ends):
        sl = slice(s, e + 1)
        draw_chain(ax, xyz[sl], mask[sl], np.asarray(cols)[sl],
                   None if plddt is None else plddt[sl], **kw)


def marker(ax, xyz, idx, color, size, label=None, zorder=6, edge=WHITE, lwx=0.7):
    if not len(idx):
        return
    p = xyz[list(idx)]
    o = np.argsort(p[:, 2])
    ax.scatter(p[o, 0], p[o, 1], s=size, facecolor=color, edgecolor=edge,
               linewidth=lwx, zorder=zorder, label=label)


def frame(ax, title=None):
    ax.set_aspect("equal")
    ax.axis("off")
    if title:
        ax.set_title(title, loc="left", color=INK, pad=5)


def fit(ax, pts, pad=0.07, aspect=1.0, bottom=0.0):
    x, y = pts[:, 0], pts[:, 1]
    cx = 0.5 * (x.min() + x.max())
    w = (x.max() - x.min()) * (1.0 + 2.0 * pad)
    h = (y.max() - y.min()) * (1.0 + 2.0 * pad)
    hi = y.max() + 0.5 * (h - (y.max() - y.min()))
    lo = hi - h * (1.0 + bottom)
    h = hi - lo
    if w / h < aspect:
        w = h * aspect
    else:
        cy = 0.5 * (lo + hi)
        h = w / aspect
        lo = cy - 0.5 * h
    ax.set_xlim(cx - 0.5 * w, cx + 0.5 * w)
    ax.set_ylim(lo, lo + h)


def bilayer(xyz, mask, aa):
    S = np.zeros((3, 3))
    for s, e in runs_of(mask, minlen=12):
        a = xyz[e] - xyz[s]
        if np.linalg.norm(a) < 22.0:
            continue
        S += np.outer(a / np.linalg.norm(a), a / np.linalg.norm(a)) * (e - s + 1)
    if not S.any():
        die("no helix long enough to set the membrane normal")
    w, v = np.linalg.eigh(S)
    n = v[:, np.argmax(w)]
    t = xyz @ n
    h = np.array([KD.get(c, 0.0) for c in aa])
    best = None
    for c in np.arange(t.min() + HALF_THICKNESS, t.max() - HALF_THICKNESS, 0.25):
        sel = np.abs(t - c) < HALF_THICKNESS
        if sel.sum() < 60:
            continue
        if best is None or h[sel].mean() > best[0]:
            best = (h[sel].mean(), c)
    if best is None:
        die("no thirty angstrom window along the normal holds enough residues")
    centre = best[1]
    out = t > centre
    seq = "".join(aa)
    sequon = [i for i in range(len(seq) - 2)
              if seq[i] == "N" and seq[i + 1] != "P" and seq[i + 2] in "ST"]
    glyco = sum(1 for i in sequon if out[i]) - sum(1 for i in sequon if not out[i])
    loop = ~mask
    basic = np.isin(aa, ["R", "K"])
    up = loop & (t > centre + HALF_THICKNESS)
    dn = loop & (t < centre - HALF_THICKNESS)
    charge = 0.0
    if up.sum() and dn.sum():
        charge = basic[up].mean() - basic[dn].mean()
    if glyco > 0 or (glyco == 0 and charge < 0):
        n, centre = -n, -centre
    return n, centre


def frame_from(n_cyto, xyz):
    y = -n_cyto
    c = xyz - xyz.mean(axis=0)
    flat = c - np.outer(c @ y, y)
    w, v = np.linalg.eigh(np.cov(flat.T))
    x = v[:, np.argmax(w)]
    x = x - np.dot(x, y) * y
    x /= np.linalg.norm(x)
    z = np.cross(x, y)
    R = np.column_stack([x, y, z])
    if np.linalg.det(R) < 0:
        R[:, 2] *= -1
    return R


def down_the_normal(R):
    x, y, z = R[:, 0], R[:, 1], R[:, 2]
    return np.column_stack([x, -z, y])


AC = jload("admissible_construct.json")
AF3 = jload("af3_models_report.json")
SPEC_P = os.path.join(ROOT, "md", "inputs", "build_spec.json")
if not os.path.isfile(SPEC_P):
    die("missing %s" % SPEC_P)
with open(SPEC_P) as fh:
    SPEC = json.load(fh)

PSET = [l.rstrip("\n").split("\t") for l in
        open(os.path.join(DATA, "pharmacophore_position_set.tsv"))]
PH = {h: i for i, h in enumerate(PSET[0])}
PRIMARIES = [int(r[PH["chicken_pos"]]) for r in PSET[1:] if r[PH["differs"]] == "1"]
RESCUERS = sorted(int(k) for k in AC["rescuer_positions"])
CONTACT = sorted(int(r[PH["chicken_pos"]]) for r in PSET[1:]
                 if r[PH["min_pemetrexed_dist_A"]] not in ("", None))
CORE_LO, CORE_HI = SPEC["residue_range"]
N_CHANGED = AC["n_positions_changed"]

RG = jload("rescue_geometry.json")
MD = os.path.join(ROOT, "data", "md")
NM_TO_A = 10.0

ARMS = ("wt", "primaries", "construct")
MODEL = {a: os.path.join(DATA, "af3", a, "fold_%s_model_0.cif" % a) for a in ARMS}

C_N, C_C = pale(S1, 0.74), pale(S4, 0.70)


def fluctuation(arm):
    import glob
    files = sorted(glob.glob(os.path.join(MD, "%s_r*_rmsf_chain.xvg" % arm)))
    if not files:
        die("no fluctuation file for the %s arm" % arm)
    per = []
    for f in files:
        rows = [l.split() for l in open(f)
                if l.strip() and l[0] not in "@#"]
        per.append({int(a): float(b) * NM_TO_A for a, b in rows})
    keys = sorted(set.intersection(*(set(d) for d in per)))
    if len(keys) < 100:
        die("the %s replicates share only %d residues" % (arm, len(keys)))
    nums = np.array(keys)
    val = np.array([np.mean([d[k] for d in per]) for k in keys])
    return nums, val, len(files)


def core(nums, *arrays):
    keep = (nums >= CORE_LO) & (nums <= CORE_HI)
    if keep.sum() < 100:
        die("the core range %d to %d keeps only %d residues"
            % (CORE_LO, CORE_HI, keep.sum()))
    return (nums[keep],) + tuple(a[keep] for a in arrays)


def halves(nums, mask):
    mid = (nums > CORE_LO + 120) & (nums < CORE_HI - 120)
    gaps = [(e - s, s, e) for s, e in runs_of((~mask) & mid, minlen=4)]
    if not gaps:
        die("no connection long enough to split the two bundles")
    _, s, e = max(gaps)
    cut = int(nums[(s + e) // 2])
    return cut, np.where(nums[:, None] < cut, C_N[None, :], C_C[None, :])


def kabsch(P, Q):
    pc, qc = P.mean(axis=0), Q.mean(axis=0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def close_gaps(mask, gap=4):
    m = mask.copy()
    for s, e in runs_of(~m):
        if e - s + 1 <= gap and s > 0 and e < len(m) - 1:
            m[s:e + 1] = True
    return m


def running_mean(a, w=5):
    k = np.ones(w) / float(w)
    return np.convolve(np.pad(a, (w // 2, w // 2), mode="edge"), k, mode="valid")


def tm_helices(xyz, mask, n_cyto, centre, minlen=12):
    t = xyz @ n_cyto
    out = []
    for s, e in runs_of(close_gaps(mask), minlen=10):
        d = np.sign(np.diff(running_mean(t[s:e + 1])))
        cuts = ([0] + [i + 1 for i in range(1, len(d))
                       if d[i] and d[i - 1] and d[i] != d[i - 1]] + [e - s + 1])
        for a, b in zip(cuts[:-1], cuts[1:]):
            if b - a < minlen:
                continue
            seg = t[s + a:s + b]
            if (seg.min() < centre - 0.6 * HALF_THICKNESS
                    and seg.max() > centre + 0.6 * HALF_THICKNESS):
                out.append((s + a, s + b - 1))
    return sorted(out)


def crossing(path2, depth, level):
    d = depth - level
    hit = np.where(np.sign(d[:-1]) * np.sign(d[1:]) < 0)[0]
    if not len(hit):
        return path2[np.argmin(np.abs(d))]
    i = hit[len(hit) // 2]
    f = abs(d[i]) / max(abs(d[i]) + abs(d[i + 1]), 1e-9)
    return path2[i] + f * (path2[i + 1] - path2[i])


def figure6():
    fig = plt.figure(figsize=(7.1, 6.6))
    gs = fig.add_gridspec(2, 2, hspace=0.17, wspace=0.10)

    nums, xyz, plddt, aa = trace(MODEL["construct"])
    nums, xyz, plddt, aa = core(nums, xyz, plddt, aa)
    mask = helix_mask(xyz)
    n_cyto, centre = bilayer(xyz, mask, aa)
    R = frame_from(n_cyto, xyz)
    origin = xyz.mean(axis=0)
    base = (xyz - origin) @ R
    mid = float(((centre * n_cyto) - origin) @ R[:, 1])
    cut, cols = halves(nums, mask)
    at = {int(n): i for i, n in enumerate(nums)}
    prim_i = [at[p] for p in PRIMARIES if p in at]
    resc_i = [at[p] for p in RESCUERS if p in at]
    if len(prim_i) != len(PRIMARIES) or len(resc_i) != len(RESCUERS):
        die("the model core does not carry every substituted position")
    tms = tm_helices(xyz, mask, n_cyto, centre)
    if len(tms) != 12:
        die("the model gives %d membrane crossing helices and not twelve" % len(tms))

    ax = fig.add_subplot(gs[0, 0])
    draw_breaks(ax, nums, base, mask, cols, plddt)
    marker(ax, base, prim_i, S2, 26, "humanising")
    marker(ax, base, resc_i, S3, 34, "rescuer", zorder=7)
    frame(ax, "a")
    fit(ax, base, bottom=0.09)
    x0, x1 = ax.get_xlim()
    ax.fill_between([x0, x1], mid - HALF_THICKNESS, mid + HALF_THICKNESS,
                    facecolor=SLAB, edgecolor="none", zorder=0)
    for yy in (mid - HALF_THICKNESS, mid + HALF_THICKNESS):
        ax.plot([x0, x1], [yy, yy], color=RULE, lw=0.7, zorder=1)
    ax.set_xlim(x0, x1)
    ax.annotate("extracellular", (x1 - 0.015 * (x1 - x0),
                                  mid + HALF_THICKNESS + 1.5), ha="right",
                va="bottom", fontsize=6.6, color=MUTED, zorder=8)
    ax.annotate("cytoplasm", (x1 - 0.015 * (x1 - x0),
                              mid - HALF_THICKNESS - 1.5), ha="right",
                va="top", fontsize=6.6, color=MUTED, zorder=8)
    ax.legend(loc="lower left", handletextpad=0.25, borderpad=0.1,
              scatterpoints=1, ncol=2, columnspacing=0.7)

    ax = fig.add_subplot(gs[0, 1])
    Rb = down_the_normal(R)
    vb = (xyz - origin) @ Rb
    path = axis_path(xyz, mask)
    pb = (path - origin) @ Rb
    inside = np.abs(vb[:, 2] - mid) < HALF_THICKNESS
    centres = []
    for k, (i0, i1) in enumerate(tms, 1):
        c = crossing(pb[i0:i1 + 1, :2], pb[i0:i1 + 1, 2], mid)
        centres.append(c)
        inner = nums[i0] < cut
        ax.add_patch(plt.Circle(c, 5.0, facecolor=C_N if inner else C_C,
                                edgecolor=pale(S1 if inner else S4, 0.95),
                                linewidth=0.8, zorder=2))
        ax.annotate(str(k), c, ha="center", va="center", fontsize=6.6,
                    color=INK, zorder=8,
                    path_effects=[pe.withStroke(linewidth=2.1, foreground=WHITE)])
    marker(ax, vb, [i for i in prim_i if inside[i]], S2, 24, zorder=5)
    marker(ax, vb, [i for i in resc_i if inside[i]], S3, 32, zorder=6)
    frame(ax, "b")
    fit(ax, np.vstack([np.array(centres) - 6.0, np.array(centres) + 6.0,
                       vb[inside][:, :2]]))
    ax.annotate("amino terminal half 1 to 6, carboxy terminal half 7 to 12",
                (0.0, 0.015), xycoords="axes fraction", fontsize=6.6, color=MUTED)

    ax = fig.add_subplot(gs[1, 0])
    st = gemmi.read_structure(os.path.join(SEQ, "7BC7.cif"))
    st.setup_entities()
    lig = [a for ch in st[0] for r in ch if r.name == "LYA" for a in r
           if a.element != gemmi.Element("H")]
    if not lig:
        die("no LYA ligand in 7BC7.cif, so the antifolate site cannot be drawn")
    lnums, lxyz, _, laa = trace(os.path.join(SEQ, "7BC7.cif"), chain="A")
    lnums, lxyz, laa = core(lnums, lxyz, laa)
    lmask = helix_mask(lxyz)
    Rl = frame_from(bilayer(lxyz, lmask, laa)[0], lxyz)
    L = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in lig])
    site = L.mean(axis=0)
    v = (lxyz - site) @ Rl
    lv = (L - site) @ Rl
    lat = {int(n): i for i, n in enumerate(lnums)}
    dist = np.linalg.norm(v, axis=1)
    near_i = dist < SITE_SHELL
    shown, drawn = [], set()
    for s, e in runs_of(dist < SITE_CONTEXT, minlen=4):
        draw_breaks(ax, lnums[s:e + 1], v[s:e + 1], lmask[s:e + 1],
                    np.tile(pale(GREY, 0.92), (e - s + 1, 1)),
                    zorder=1, lw_helix=5.2, fog=0.55)
        shown.append(v[s:e + 1])
        drawn.update(range(s, e + 1))
    if not shown:
        die("nothing of the deposited chain falls within the site shell")
    for i in range(len(L)):
        for j in range(i + 1, len(L)):
            if np.linalg.norm(L[i] - L[j]) < 1.75:
                ax.plot(lv[[i, j], 0], lv[[i, j], 1], color=LIG, lw=2.6,
                        solid_capstyle="round", zorder=3)
    o = np.argsort(lv[:, 2])
    ax.scatter(lv[o, 0], lv[o, 1], s=10, facecolor=LIG, edgecolor=LIG,
               linewidth=0.5, zorder=4, label="pemetrexed")
    near = [p for p in CONTACT if p in lat and near_i[lat[p]]]
    adrift = [p for p in near if lat[p] not in drawn]
    if adrift:
        die("the set positions %s fall within the site shell with no backbone "
            "drawn around them" % ", ".join(str(q) for q in adrift))
    marker(ax, v, [lat[p] for p in near if p not in PRIMARIES], WHITE, 24,
           "conserved", zorder=5, edge=MUTED, lwx=0.8)
    marker(ax, v, [lat[p] for p in near if p in PRIMARIES], S2, 32,
           "humanising", zorder=6)
    rin = [lat[p] for p in RESCUERS if p in lat and near_i[lat[p]]]
    marker(ax, v, rin, S3, 32, "rescuer" if rin else None, zorder=7)
    taken, unlabelled = [], 0
    for p in sorted(near):
        if p not in PRIMARIES:
            continue
        q = v[lat[p]]
        put = None
        for dx, dy in ((4.2, 3.4), (-4.2, 3.4), (4.2, -6.6), (-4.2, -6.6)):
            a, b = q[0] + dx * 0.18, q[1] + dy * 0.18
            if not any(abs(a - u) < 2.5 and abs(b - w) < 1.8 for u, w in taken):
                put = (dx, dy, a, b)
                break
        if put is None:
            unlabelled += 1
            continue
        taken.append((put[2], put[3]))
        ax.annotate(str(p), (q[0], q[1]), fontsize=6.3, color=S2, weight="bold",
                    ha="left" if put[0] > 0 else "right",
                    textcoords="offset points", xytext=(put[0], put[1]), zorder=8,
                    path_effects=[pe.withStroke(linewidth=1.9, foreground=WHITE)])
    if unlabelled:
        die("%d humanising positions in the site shell could not be labelled"
            % unlabelled)
    frame(ax, "c")
    fit(ax, np.vstack([v[near_i], lv]), pad=0.10, bottom=0.11)
    ax.legend(loc="lower left", handletextpad=0.25, borderpad=0.1,
              scatterpoints=1, ncol=2, columnspacing=0.7)

    ax = fig.add_subplot(gs[1, 1])
    ref_n, ref_x, _, ref_a = trace(MODEL["wt"])
    ref_n, ref_x, ref_a = core(ref_n, ref_x, ref_a)
    ref_m = helix_mask(ref_x)
    Rw = frame_from(bilayer(ref_x, ref_m, ref_a)[0], ref_x)
    arm_col = {"wt": rgb(GREY), "primaries": rgb(S2), "construct": rgb(S1)}
    allpts = []
    for a in ARMS:
        n, x, _, _ = trace(MODEL[a])
        n, x = core(n, x)
        m = helix_mask(x)
        shared = sorted(set(n.tolist()) & set(ref_n.tolist()))
        ia = {int(k): i for i, k in enumerate(n)}
        ir = {int(k): i for i, k in enumerate(ref_n)}
        Rr, tt = kabsch(x[[ia[s] for s in shared]], ref_x[[ir[s] for s in shared]])
        v = ((x @ Rr.T + tt) - ref_x.mean(axis=0)) @ Rw
        draw_breaks(ax, n, v, m, np.tile(arm_col[a], (len(n), 1)),
                    zorder=2 if a == "wt" else 3,
                    lw_helix=4.6 if a == "wt" else 3.0, lw_loop=1.5,
                    fog=0.60, alpha=0.90)
        allpts.append(v)
    frame(ax, "d")
    fit(ax, np.vstack(allpts), bottom=0.13)
    ax.annotate(u"no pair further apart than %.2f \u00c5\nover the %d shared positions"
                % (AF3["furthest_arm_pair"], AF3["core_positions"]),
                (0.0, 0.085), xycoords="axes fraction", fontsize=6.4, color=MUTED,
                linespacing=1.45)
    ax.legend(handles=[Line2D([], [], color=arm_col[a], lw=2.6,
                              label={"wt": "wild type"}.get(a, a)) for a in ARMS],
              loc="lower left", handletextpad=0.5, borderpad=0.1, ncol=3,
              columnspacing=0.8)

    p = os.path.join(FIGS, "fig6_structures.png")
    fig.savefig(p, dpi=DPI, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(os.path.join(FIGS, "fig6_structures.svg"), bbox_inches="tight",
                pad_inches=0.04)
    plt.close(fig)

    rep = os.path.join(DATA, "structure_figure.json")
    with open(rep, "w") as fh:
        json.dump({
            "model": os.path.relpath(MODEL["construct"], ROOT),
            "residue_range": [CORE_LO, CORE_HI],
            "half_thickness_A": HALF_THICKNESS,
            "bilayer_centre_A": round(float(centre), 2),
            "membrane_crossings": len(tms),
            "crossing_spans": [[int(nums[a]), int(nums[b])] for a, b in tms],
            "crossings_amino_half": sum(1 for a, _ in tms if nums[a] < cut),
            "split_residue": int(cut),
            "plddt_floor": PLDDT_FLOOR,
            "plddt_below_floor": int((plddt < PLDDT_FLOOR).sum()),
            "site_shell_A": SITE_SHELL,
            "site_positions_shown": len(near),
            "site_positions_humanising": len(taken),
            "site_positions_rescuer": len(rin),
            "reference": "7BC7",
        }, fh, indent=1, sort_keys=True)

    print("  bilayer centre %.1f, split at residue %d, %d membrane crossings"
          % (centre, cut, len(tms)))
    print("  wrote %s" % os.path.relpath(rep, ROOT))
    print("  fig6_structures %7.0f kB" % (os.path.getsize(p) / 1024.0))
    return p


def figure7():
    fig = plt.figure(figsize=(7.1, 6.9))
    gs = fig.add_gridspec(2, 2, hspace=0.30, wspace=0.12)

    nums, xyz, plddt, aa = trace(MODEL["construct"])
    nums, xyz, plddt, aa = core(nums, xyz, plddt, aa)
    mask = helix_mask(xyz)
    n_cyto, centre = bilayer(xyz, mask, aa)
    R = frame_from(n_cyto, xyz)
    origin = xyz.mean(axis=0)
    base = (xyz - origin) @ R
    mid = float(((centre * n_cyto) - origin) @ R[:, 1])
    at = {int(n): i for i, n in enumerate(nums)}
    tms = tm_helices(xyz, mask, n_cyto, centre)
    if len(tms) != 12:
        die("the model gives %d membrane crossings and not twelve" % len(tms))
    faint = np.tile(pale(GREY, 1.00), (len(nums), 1))

    ax = fig.add_subplot(gs[0, 0])
    draw_breaks(ax, nums, base, mask, faint, zorder=1, fog=0.46)
    for q in RG["pairs"]:
        t, r = q["target"], q["rescuer"]
        if t not in at or r not in at:
            die("the pair %d and %d is outside the modelled core" % (t, r))
        p0, p1 = base[at[t]], base[at[r]]
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=INK2, lw=0.9,
                zorder=4, solid_capstyle="round")
    marker(ax, base, [at[q["target"]] for q in RG["pairs"]], S2, 34, "rescue target")
    marker(ax, base, [at[q["rescuer"]] for q in RG["pairs"]], S3, 34, "rescuer",
           zorder=7)
    far = max(RG["pairs"], key=lambda q: q["separation"])
    q0, q1 = base[at[far["target"]]], base[at[far["rescuer"]]]
    ax.annotate("%d residues apart" % far["separation"],
                (0.5 * (q0[0] + q1[0]), 0.5 * (q0[1] + q1[1])),
                fontsize=6.3, color=INK2, ha="center",
                textcoords="offset points", xytext=(0, 5), zorder=8,
                path_effects=[pe.withStroke(linewidth=2.1, foreground=WHITE)])
    frame(ax, "a")
    fit(ax, base, bottom=0.09)
    x0, x1 = ax.get_xlim()
    ax.fill_between([x0, x1], mid - HALF_THICKNESS, mid + HALF_THICKNESS,
                    facecolor=SLAB, edgecolor="none", zorder=0)
    ax.set_xlim(x0, x1)
    ax.legend(loc="lower left", handletextpad=0.25, borderpad=0.1,
              scatterpoints=1, ncol=2, columnspacing=0.7)

    ax = fig.add_subplot(gs[0, 1])
    ax.set_aspect("auto")
    sep = np.array([q["separation"] for q in RG["pairs"]], float)
    dis = np.array([q["distance_max"] for q in RG["pairs"]], float)
    out = np.array([q["beyond_window"] for q in RG["pairs"]])
    ax.axvspan(0.7, RG["window"], facecolor=SLAB, edgecolor="none", zorder=0)
    ax.axhline(RG["cutoff_A"], color=MUTED, lw=0.8, ls=(0, (3, 2)), zorder=1)
    ax.scatter(sep[~out], dis[~out], s=42, facecolor=S3, edgecolor=WHITE,
               linewidth=0.7, zorder=3, label="inside the window")
    ax.scatter(sep[out], dis[out], s=52, facecolor=S2, edgecolor=WHITE,
               linewidth=0.7, zorder=4, label="outside it")
    at_point = {}
    for q in RG["pairs"]:
        at_point.setdefault((q["separation"], q["distance_max"]), []).append(
            "%d to %d" % (q["rescuer"], q["target"]))
    for (x, y), names in at_point.items():
        ax.annotate(", ".join(sorted(names)), (x, y), fontsize=6.2, color=INK2,
                    textcoords="offset points", xytext=(6, 3), zorder=5,
                    path_effects=[pe.withStroke(linewidth=2.0, foreground=WHITE)])
    ax.set_xscale("log")
    ax.set_xlim(0.7, 140)
    ax.set_ylim(0, RG["cutoff_A"] * 1.22)
    ax.set_xticks([1, 3, 10, 30, 100])
    ax.set_xticklabels(["1", "3", "10", "30", "100"])
    ax.set_xlabel("separation along the sequence, residues", labelpad=2)
    ax.set_ylabel(u"minimum heavy-atom distance, \u00c5", labelpad=2)
    ax.set_title("b", loc="left", color=INK, pad=5)
    ax.annotate("contact cutoff", (130, RG["cutoff_A"]), fontsize=6.4,
                color=MUTED, ha="right", va="bottom")
    ax.annotate("published window", (RG["window"], 0.06 * RG["cutoff_A"]),
                fontsize=6.4, color=MUTED, ha="left", va="bottom",
                textcoords="offset points", xytext=(3, 0))
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_linewidth(0.7); ax.spines[sp].set_color(MUTED)
    ax.legend(loc="upper left", handletextpad=0.25, borderpad=0.1,
              scatterpoints=1)

    ax = fig.add_subplot(gs[1, 0])
    WATCH = sorted(AF3["watched_positions"])
    sixth, seventh = tms[5], tms[6]
    cols = np.tile(pale(GREY, 0.62), (len(nums), 1))
    cols[sixth[0]:sixth[1] + 1] = rgb(S1)
    link = (nums > nums[sixth[1]]) & (nums < nums[seventh[0]])
    cols[link] = rgb(S4)
    draw_breaks(ax, nums, base, mask, cols, zorder=2, fog=0.58)
    marker(ax, base, [at[w] for w in WATCH if w in at], S2, 40, zorder=7)
    for w in WATCH:
        if w in at:
            ax.annotate(str(w), tuple(base[at[w]][:2]), fontsize=6.6, color=S2,
                        weight="bold", textcoords="offset points",
                        xytext=(5, 3), zorder=8,
                        path_effects=[pe.withStroke(linewidth=2.1,
                                                    foreground=WHITE)])
    frame(ax, "c")
    fit(ax, base, bottom=0.11)
    x0, x1 = ax.get_xlim()
    ax.fill_between([x0, x1], mid - HALF_THICKNESS, mid + HALF_THICKNESS,
                    facecolor=SLAB, edgecolor="none", zorder=0)
    ax.set_xlim(x0, x1)
    ax.legend(handles=[Line2D([], [], color=S1, lw=3.4,
                              label="the sixth membrane crossing"),
                       Line2D([], [], color=S4, lw=3.4,
                              label="the connection between the halves")],
              loc="lower left", handletextpad=0.5, borderpad=0.1, ncol=1,
              labelspacing=0.25)

    ax = fig.add_subplot(gs[1, 1])
    fn, fv, nrep = fluctuation("construct")
    shared = np.array([int(n) in at for n in fn])
    if shared.sum() < 100:
        die("the fluctuation file shares only %d residues with the model"
            % int(shared.sum()))
    order = [at[int(n)] for n in fn[shared]]
    val = fv[shared]
    lo, hi = float(np.percentile(val, 2)), float(np.percentile(val, 98))
    ramp = LinearSegmentedColormap.from_list("flux", [pale(GREY, 0.60), S2])
    fcols = faint.copy()
    fcols[order] = ramp(np.clip((val - lo) / max(hi - lo, 1e-9), 0, 1))[:, :3]
    draw_breaks(ax, nums, base, mask, fcols, zorder=2, fog=0.50)
    frame(ax, "d")
    fit(ax, base, bottom=0.11)
    x0, x1 = ax.get_xlim()
    ax.fill_between([x0, x1], mid - HALF_THICKNESS, mid + HALF_THICKNESS,
                    facecolor=SLAB, edgecolor="none", zorder=0)
    ax.set_xlim(x0, x1)
    cax = ax.inset_axes([0.03, 0.055, 0.42, 0.028])
    cax.imshow(np.linspace(0, 1, 128)[None, :], aspect="auto", cmap=ramp)
    cax.set_xticks([0, 127]); cax.set_yticks([])
    cax.set_xticklabels([u"%.1f \u00c5" % lo, u"%.1f" % hi], fontsize=6.4,
                        color=INK2)
    cax.tick_params(length=0, pad=1.5)
    for sp in cax.spines.values():
        sp.set_linewidth(0.4); sp.set_color(MUTED)
    ax.annotate("construct, mean of %d replicates" % nrep, (0.03, 0.135),
                xycoords="axes fraction", fontsize=6.4, color=MUTED)

    p = os.path.join(FIGS, "fig7_rescue_geometry.png")
    fig.savefig(p, dpi=DPI, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(os.path.join(FIGS, "fig7_rescue_geometry.svg"),
                bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)

    rep = os.path.join(DATA, "junction_figure.json")
    with open(rep, "w") as fh:
        json.dump({
            "watched": [int(w) for w in WATCH],
            "crossing_six": [int(nums[sixth[0]]), int(nums[sixth[1]])],
            "crossing_seven": [int(nums[seventh[0]]), int(nums[seventh[1]])],
            "connection": [int(nums[sixth[1]]) + 1, int(nums[seventh[0]]) - 1],
            "connection_length": int(nums[seventh[0]]) - int(nums[sixth[1]]) - 1,
            "fluctuation_replicates": nrep,
            "fluctuation_low": round(lo, 2),
            "fluctuation_high": round(hi, 2),
        }, fh, indent=1, sort_keys=True)
    print("  connection between the halves %d to %d, fluctuation %.2f to %.2f A"
          % (int(nums[sixth[1]]) + 1, int(nums[seventh[0]]) - 1, lo, hi))
    print("  wrote %s" % os.path.relpath(rep, ROOT))
    print("  fig7_rescue_geometry %7.0f kB" % (os.path.getsize(p) / 1024.0))
    return p


def toc_graphic():
    CM = 1.0 / 2.54
    W, H = 8.9 * CM, 3.6 * CM
    fig = plt.figure(figsize=(W, H))
    fig.patch.set_facecolor(WHITE)

    nums, xyz, plddt, aa = trace(MODEL["construct"])
    nums, xyz, plddt, aa = core(nums, xyz, plddt, aa)
    mask = helix_mask(xyz)
    n_cyto, centre = bilayer(xyz, mask, aa)
    R = frame_from(n_cyto, xyz)
    origin = xyz.mean(axis=0)
    mid = float(((centre * n_cyto) - origin) @ R[:, 1])
    at = {int(n): i for i, n in enumerate(nums)}
    prim_i = [at[q] for q in PRIMARIES if q in at]
    resc_i = [at[q] for q in RESCUERS if q in at]

    ref_n, ref_x, _, ref_a = trace(MODEL["wt"])
    ref_n, ref_x, ref_a = core(ref_n, ref_x, ref_a)
    ref_m = helix_mask(ref_x)
    shared = sorted(set(ref_n.tolist()) & set(nums.tolist()))
    ia = {int(k): i for i, k in enumerate(ref_n)}
    ir = {int(k): i for i, k in enumerate(nums)}
    Rr, tt = kabsch(ref_x[[ia[q] for q in shared]], xyz[[ir[q] for q in shared]])
    ref_v = ((ref_x @ Rr.T + tt) - origin) @ R
    base = (xyz - origin) @ R

    panels = ((0.035, ref_n, ref_v, ref_m, None, None, "chicken scaffold"),
              (0.56, nums, base, mask, prim_i, resc_i, "humanised surrogate"))
    for x0, nn, vv, mm, pi, ri, label in panels:
        ax = fig.add_axes([x0, 0.16, 0.40, 0.80])
        draw_breaks(ax, nn, vv, mm, np.tile(pale(GREY, 0.95), (len(nn), 1)),
                    zorder=2, fog=0.50, lw_helix=2.4, lw_loop=0.9)
        if pi:
            marker(ax, vv, pi, S2, 5.5, zorder=5, lwx=0.25)
            marker(ax, vv, ri, S3, 7.0, zorder=6, lwx=0.25)
        frame(ax)
        fit(ax, np.vstack([base, ref_v]), pad=0.03)
        xl = ax.get_xlim()
        ax.fill_between(xl, mid - HALF_THICKNESS, mid + HALF_THICKNESS,
                        facecolor=SLAB, edgecolor="none", zorder=0)
        ax.set_xlim(xl)
        ax.annotate(label, (0.5, -0.015), xycoords="axes fraction", ha="center",
                    va="top", fontsize=5.4, color=INK2)

    ax = fig.add_axes([0.455, 0.16, 0.09, 0.80])
    ax.axis("off")
    ax.annotate("", (0.92, 0.56), xytext=(0.08, 0.56), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color=INK2, linewidth=0.8,
                                shrinkA=0, shrinkB=0))
    ax.annotate("%d changes" % N_CHANGED, (0.5, 0.63), xycoords="axes fraction",
                ha="center", va="bottom", fontsize=5.4, color=INK2)
    ax.annotate("%d rescuers" % len(RESCUERS), (0.5, 0.49),
                xycoords="axes fraction", ha="center", va="top", fontsize=5.4,
                color=INK2)

    p = os.path.join(FIGS, "toc_graphic.png")
    fig.savefig(p, dpi=600, facecolor=WHITE)
    plt.close(fig)
    from PIL import Image
    w, h = Image.open(p).size
    want = (round(8.9 * CM * 600), round(3.6 * CM * 600))
    if (w, h) != want:
        die("the graphic came out %d by %d pixels and the print size asks for "
            "%d by %d" % (w, h, want[0], want[1]))
    print("  toc_graphic %d by %d px, %.2f by %.2f cm at 600 dpi"
          % (w, h, w / 600.0 * 2.54, h / 600.0 * 2.54))
    return p


if __name__ == "__main__":
    for build in (figure6, figure7, toc_graphic):
        out = build()
        if not os.path.isfile(out):
            die("%s did not write" % build.__name__)
        print("wrote %s" % os.path.relpath(out, ROOT))
