import collections
import itertools
import json

import gemmi

CUTOFFS = (4.0, 5.0, 6.0, 8.0)
WINDOW = 3

rows = [l.rstrip('\n').split('\t') for l in open('data/pharmacophore_position_set.tsv')]
col = {k: i for i, k in enumerate(rows[0])}
prim = [(int(r[col['chicken_pos']]), r[col['chicken_aa']],
         int(r[col['human_pos']]), r[col['human_aa']])
        for r in rows[1:] if r[col['differs']] == '1']
prim.sort()


def distances(entry):
    s = gemmi.read_structure(f'seq/{entry}.cif')
    s.setup_entities()
    s.remove_waters()
    res = {r.seqid.num: r for r in s[0]['A']}
    absent = [p for p, _, _, _ in prim if p not in res]
    if absent:
        raise SystemExit(f'FATAL not modelled in {entry}: {absent}')
    return {(p, q): min(x.pos.dist(y.pos) for x in res[p] for y in res[q])
            for (p, _, _, _), (q, _, _, _) in itertools.combinations(prim, 2)}


d6, d7 = distances('7BC6'), distances('7BC7')
out = {'n_primary': len(prim), 'n_pairs': len(d6), 'window': WINDOW, 'cutoffs': {}}

for cut in CUTOFFS:
    edges = sorted(k for k in d6 if d6[k] <= cut and d7[k] <= cut)
    inside = [e for e in edges if abs(e[0] - e[1]) <= WINDOW]
    g = collections.defaultdict(set)
    for p, q in edges:
        g[p].add(q)
        g[q].add(p)
    seen, comps = set(), []
    for p, _, _, _ in prim:
        if p in seen:
            continue
        stack, comp = [p], set()
        while stack:
            x = stack.pop()
            if x in comp:
                continue
            comp.add(x)
            seen.add(x)
            stack.extend(g[x] - comp)
        comps.append(sorted(comp))
    comps.sort(key=len, reverse=True)
    out['cutoffs'][str(cut)] = {
        'coupled_pairs': len(edges),
        'pairs_inside_window': len(inside),
        'pairs_outside_window': len(edges) - len(inside),
        'max_separation': max((abs(p - q) for p, q in edges), default=0),
        'positions_coupled': sum(len(c) for c in comps if len(c) > 1),
        'isolated_positions': [c[0] for c in comps if len(c) == 1],
        'components': [c for c in comps if len(c) > 1],
        'edges': [[p, q, abs(p - q), round(d6[(p, q)], 2), round(d7[(p, q)], 2)]
                  for p, q in edges],
    }

json.dump(out, open('data/primary_coupling.json', 'w'), indent=1)
for cut in CUTOFFS:
    c = out['cutoffs'][str(cut)]
    print(f'{cut:4.1f} A  pairs {c["coupled_pairs"]:3d}  '
          f'outside +/-{WINDOW} {c["pairs_outside_window"]:3d}  '
          f'max sep {c["max_separation"]:3d}  '
          f'coupled positions {c["positions_coupled"]:2d}/{len(prim)}  '
          f'components {[len(x) for x in c["components"]]}  '
          f'isolated {c["isolated_positions"]}')
