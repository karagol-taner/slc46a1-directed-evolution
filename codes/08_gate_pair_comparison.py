import json
import os
import statistics
import sys

SCAN = 'data/esm2_step2/gate_perturbation.tsv'


def die(msg):
    sys.exit('FATAL ' + msg)


if not os.path.exists(SCAN):
    die('%s not found; run codes/SLC46A1_step2_rescue_scan.ipynb on Colab and copy '
        'esm2_step2/gate_perturbation.tsv into data/' % SCAN)

gates = json.load(open('data/zhao_gate_chicken.json'))
z = json.load(open('data/zhao2016_gate.json'))
h2c = {g['human_pos']: g['chicken_pos'] for g in gates}
gate_pos = {g['chicken_pos'] for g in gates}

ent = {}
for rec in open('seq/7BC6.fasta').read().split('>'):
    if rec.strip():
        h, *r = rec.strip().split('\n')
        ent[h.split('|')[0]] = ''.join(r)
chicken = ent['7BC6_1'][:473]
for g in gates:
    if chicken[g['chicken_pos'] - 1] != g['chicken_aa']:
        die('gate chicken %d is %s in sequence, %s in file'
            % (g['chicken_pos'], chicken[g['chicken_pos'] - 1], g['chicken_aa']))

rows = [l.rstrip('\n').split('\t') for l in open(SCAN)]
head, body = rows[0], rows[1:]
c = {k: i for i, k in enumerate(head)}
if not body:
    die('scan file has no data rows')

positions = sorted({int(r[c['variant_pos']]) for r in body})
expected = len(positions) * 19 * len(gates)
if len(body) != expected:
    die('scan has %d rows, expected %d for %d positions and %d gates'
        % (len(body), expected, len(positions), len(gates)))
print(f'gate scan rows {len(body)}, positions {len(positions)}, gates {len(gates)}')

by_gate = {}
for g in gates:
    gp = g['chicken_pos']
    vs = [(float(r[c['p_gate_wt']]), int(r[c['variant_pos']]), r[c['variant_to']])
          for r in body if int(r[c['gate_pos']]) == gp]
    vs.sort()
    rank = {(j, y): i + 1 for i, (p, j, y) in enumerate(vs)}
    per_pos = {}
    for (j, y), rk in rank.items():
        per_pos.setdefault(j, []).append(rk)
    score = {j: statistics.median(v) for j, v in per_pos.items()}
    by_gate[gp] = {'rank': rank, 'score': score, 'n': len(vs)}

PAIRS = [('spontaneous disulfide', p) for p in z['spontaneous_disulfide']] + \
        [('CuPh and cadmium induced', p) for p in z['induced_by_cuph_and_cadmium']]

print('\nmeasured pairs, each member scored against its partner')
print('pair            kind                        sep   partner percentile  best rank')
out = []
for kind, label in PAIRS:
    a, b = label.split('/')
    ha, hb = int(a[1:-1]), int(b[1:-1])
    if ha not in h2c or hb not in h2c:
        die('pair %s has a member outside the gate set' % label)
    ca, cb = h2c[ha], h2c[hb]
    for src, dst in ((ca, cb), (cb, ca)):
        g = by_gate[dst]
        scores = sorted(g['score'].values())
        s = g['score'][src]
        pct = 100.0 * sum(1 for v in scores if v <= s) / len(scores)
        best = min(rk for (j, y), rk in g['rank'].items() if j == src)
        out.append({'pair': label, 'kind': kind,
                    'substituted_chicken': src, 'read_chicken': dst,
                    'separation': abs(src - dst),
                    'median_rank': s, 'percentile': round(pct, 2),
                    'best_rank': best, 'n_variants': g['n']})
        print(f"{label:<15} {kind:<26} {abs(src-dst):4d}  "
              f"{pct:17.2f}  {best:9d}")

print('\nfor reference, the same score for every gate pair including the unmeasured one')
for g1 in gates:
    for g2 in gates:
        if g1['chicken_pos'] >= g2['chicken_pos']:
            continue
        g = by_gate[g2['chicken_pos']]
        scores = sorted(g['score'].values())
        s = g['score'][g1['chicken_pos']]
        pct = 100.0 * sum(1 for v in scores if v <= s) / len(scores)
        measured = any(set(p.split('/')) ==
                       {f"{g1['human_aa']}{g1['human_pos']}C",
                        f"{g2['human_aa']}{g2['human_pos']}C"}
                       for _, p in PAIRS)
        print(f"  chicken {g1['chicken_pos']:3d} on {g2['chicken_pos']:3d}  "
              f"sep {abs(g1['chicken_pos']-g2['chicken_pos']):3d}  "
              f"percentile {pct:6.2f}  measured {measured}")

json.dump({'source': SCAN, 'pairs': out,
           'statistic': 'median rank of a position\'s 19 substitutions, among all '
                        'positions, ranking variants by the wild-type probability at '
                        'the partner gate position, most disruptive first'},
          open('data/gate_pair_comparison.json', 'w'), indent=1)
print('\nwrote data/gate_pair_comparison.json')
