import json
import os
import sys

import gemmi

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
DATA = os.path.join(ROOT, "data")
SEQ = os.path.join(ROOT, "seq")

ENTRIES = ("7BC6", "7BC7")
CUTOFF_A = 6.0
WINDOW = 3


def die(msg):
    sys.exit("STOP: %s" % msg)


def residues(entry):
    p = os.path.join(SEQ, "%s.cif" % entry)
    if not os.path.isfile(p):
        die("missing %s" % p)
    s = gemmi.read_structure(p)
    s.setup_entities()
    s.remove_waters()
    return {r.seqid.num: r for r in s[0]["A"]}


def main():
    ac = os.path.join(DATA, "admissible_construct.json")
    if not os.path.isfile(ac):
        die("missing %s" % ac)
    with open(ac) as fh:
        chosen = json.load(fh)["chosen"]
    if not chosen:
        die("the construct carries no rescue, so there is no geometry to report")

    res = {e: residues(e) for e in ENTRIES}
    pairs = []
    for target in sorted(chosen, key=int):
        t = int(target)
        for r in chosen[target]["positions"]:
            absent = [e for e in ENTRIES if t not in res[e] or r not in res[e]]
            if absent:
                die("positions %d and %d are not both modelled in %s"
                    % (t, r, ", ".join(absent)))
            d = {e: min(x.pos.dist(y.pos)
                        for x in res[e][t] for y in res[e][r]) for e in ENTRIES}
            pairs.append({
                "target": t,
                "rescuer": r,
                "label": chosen[target]["label"],
                "separation": abs(t - r),
                "distance": {e: round(d[e], 2) for e in ENTRIES},
                "distance_max": round(max(d.values()), 2),
                "touching": all(v < CUTOFF_A for v in d.values()),
                "beyond_window": abs(t - r) > WINDOW,
            })

    out = {
        "cutoff_A": CUTOFF_A,
        "window": WINDOW,
        "entries": list(ENTRIES),
        "pairs": pairs,
        "n_pairs": len(pairs),
        "n_targets": len({p["target"] for p in pairs}),
        "n_touching": sum(1 for p in pairs if p["touching"]),
        "n_beyond_window": sum(1 for p in pairs if p["beyond_window"]),
        "separation_min": min(p["separation"] for p in pairs),
        "separation_max": max(p["separation"] for p in pairs),
        "distance_min": round(min(p["distance_max"] for p in pairs), 2),
        "distance_max": round(max(p["distance_max"] for p in pairs), 2),
    }
    p = os.path.join(DATA, "rescue_geometry.json")
    with open(p, "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)

    print("%-7s %-8s %-6s %s" % ("target", "rescuer", "sep", "  ".join(ENTRIES)))
    for q in pairs:
        print("%-7d %-8d %-6d %s%s"
              % (q["target"], q["rescuer"], q["separation"],
                 "  ".join("%5.2f" % q["distance"][e] for e in ENTRIES),
                 "" if q["touching"] else "   not touching"))
    print("\n%d of %d pairs touch in both entries within %.1f A"
          % (out["n_touching"], out["n_pairs"], CUTOFF_A))
    print("%d of %d lie beyond the published window of plus or minus %d"
          % (out["n_beyond_window"], out["n_pairs"], WINDOW))
    print("wrote %s" % os.path.relpath(p, ROOT))


if __name__ == "__main__":
    main()
