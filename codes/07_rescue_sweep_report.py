import json
import os
import sys

SCAN = 'data/esm2_step2/rescue_scan.tsv'
THRESHOLD = 1e-3
RADII = [3, 5, 8, 12, 20, 35, 60, None]
POLICIES = ['published', 'strict', 'annotated']


def die(msg):
    sys.exit('FATAL ' + msg)


if not os.path.exists(SCAN):
    die('%s not found; run codes/SLC46A1_step2_rescue_scan.ipynb on Colab and '
        'copy esm2_step2/rescue_scan.tsv into data/' % SCAN)

m = json.load(open('structures/7BC6_to_human.json'))
ck2hu = {int(k): v for k, v in m['ck2hu'].items()}
hu2ck = {int(k): v for k, v in m['hu2ck'].items()}

human = ''.join(l.strip() for l in open('seq/Q96NT5.fasta') if not l.startswith('>'))
ent = {}
for rec in open('seq/7BC6.fasta').read().split('>'):
    if rec.strip():
        h, *r = rec.strip().split('\n')
        ent[h.split('|')[0]] = ''.join(r)
chicken = ent['7BC6_1'][:473]

parker = json.load(open('data/parker_residues.json'))['residues']
protected = set()
for k, v in parker.items():
    cp = int(k)
    hp = v['human']
    if chicken[cp - 1] == human[hp - 1]:
        protected.add(cp)
hfm_rows = [l.rstrip('\n').split('\t') for l in open('data/hfm_positions.tsv')][1:]
unmapped = []
for r in hfm_rows:
    hp = int(r[0])
    cp = hu2ck.get(hp)
    if cp is None:
        unmapped.append(hp)
    else:
        protected.add(cp)
if unmapped:
    die('hereditary folate malabsorption positions without a chicken partner: %s'
        % unmapped)
print(f'protected positions for the annotated policy: {len(protected)}')

rows = [l.rstrip('\n').split('\t') for l in open(SCAN)]
head, body = rows[0], rows[1:]
c = {k: i for i, k in enumerate(head)}
if not body:
    die('scan file has no data rows')

targets = sorted({int(r[c['target_pos']]) for r in body})
baseline = {}
for r in body:
    tp = int(r[c['target_pos']])
    baseline.setdefault(tp, float(r[c['p_baseline']]))
for tp, b in baseline.items():
    if b >= THRESHOLD:
        die('target %d baseline %.3e is not below the threshold' % (tp, b))

n_pos = len(set(int(r[c['rescuer_pos']]) for r in body))
expected = n_pos * 19 * len(targets)
if len(body) != expected:
    die('scan has %d rows, expected %d for %d positions and %d targets'
        % (len(body), expected, n_pos, len(targets)))
print(f'scan rows {len(body)}, rescuer positions {n_pos}, targets {len(targets)}')


def allowed(policy, jpos, target_pos):
    if jpos == target_pos:
        return False
    if policy == 'published':
        return True
    hj = ck2hu.get(jpos)
    identical = hj is not None and chicken[jpos - 1] == human[hj - 1]
    if policy == 'strict':
        return not identical
    if policy == 'annotated':
        return jpos not in protected
    die('unknown policy ' + policy)


best = {}
for r in body:
    j = int(r[c['rescuer_pos']])
    tp = int(r[c['target_pos']])
    p = float(r[c['p_target']])
    sep = int(r[c['separation']])
    key = (j, r[c['rescuer_to']], tp)
    best[key] = (p, sep)

print('\nradius  policy      targets rescued  rescuers crossing  nearest rescuer')
summary = []
for rad in RADII:
    for pol in POLICIES:
        rescued, crossing, nearest = set(), 0, None
        for (j, y, tp), (p, sep) in best.items():
            if rad is not None and sep > rad:
                continue
            if not allowed(pol, j, tp):
                continue
            if p >= THRESHOLD:
                crossing += 1
                rescued.add(tp)
                nearest = sep if nearest is None else min(nearest, sep)
        label = 'all' if rad is None else str(rad)
        summary.append({'radius': label, 'policy': pol,
                        'targets_rescued': len(rescued),
                        'rescuers_crossing': crossing,
                        'nearest_rescuer_separation': nearest,
                        'rescued_positions': sorted(rescued)})
        print(f'{label:>6}  {pol:<10}  {len(rescued):>15}  {crossing:>17}  '
              f'{"none" if nearest is None else nearest}')

per_target = {}
for tp in targets:
    for pol in POLICIES:
        cands = sorted(((p, j, y, sep) for (j, y, t), (p, sep) in best.items()
                        if t == tp and allowed(pol, j, tp)), reverse=True)[:5]
        per_target.setdefault(tp, {})[pol] = [
            {'p_target': p, 'rescuer_pos': j, 'rescuer_to': y, 'separation': sep,
             'crosses': p >= THRESHOLD, 'fold_change': p / baseline[tp]}
            for p, j, y, sep in cands]

print('\nstrongest five rescuers per target, by policy, target position excluded')
for tp in targets:
    print(f'  target {chicken[tp-1]}{tp}, baseline {baseline[tp]:.2e}')
    for pol in POLICIES:
        rows = per_target[tp][pol]
        if not rows:
            print(f'    {pol:<10} no candidate positions'); continue
        head = rows[0]
        extra = sum(1 for e in rows if e['crosses'])
        print(f"    {pol:<10} best {chicken[head['rescuer_pos']-1]}"
              f"{head['rescuer_pos']}{head['rescuer_to']}"
              f"  sep {head['separation']:3d}  p {head['p_target']:.3e}"
              f"  x{head['fold_change']:.1f}  crosses {head['crosses']}"
              f"  ({extra} of the top 5 cross)")

json.dump({'threshold': THRESHOLD, 'protected_count': len(protected),
           'protected': sorted(protected), 'summary': summary,
           'per_target': {str(k): v for k, v in per_target.items()},
           'baseline': {str(k): v for k, v in baseline.items()}},
          open('data/rescue_sweep_report.json', 'w'), indent=1)
print('\nwrote data/rescue_sweep_report.json')
