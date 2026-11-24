import json
import os
import sys

POSSET = 'data/pharmacophore_position_set.tsv'
SWEEP = 'data/rescue_sweep_report.json'
D4 = 'data/double_rescue_report.json'
C5 = 'data/construct_report.json'
C6 = 'data/construct_rescue_report.json'
S8 = 'data/step8_report.json'
PARKER = 'data/parker_residues.json'
HFM = 'data/hfm_positions.tsv'
MAP = 'structures/7BC6_to_human.json'
OUT = 'data/admissible_construct.json'
LOG = 'data/20_admissible_construct_log.txt'
THRESHOLD = 1e-3

lines = []


def log(m):
    print(m, flush=True)
    lines.append(str(m))


def die(m):
    with open(LOG, 'w') as f:
        f.write('\n'.join(lines + ['FATAL ' + m]) + '\n')
    sys.exit('FATAL ' + m)


for p in (POSSET, SWEEP, D4, C5, C6, PARKER, HFM, MAP):
    if not os.path.exists(p):
        die(p + ' not found')

ent = {}
for rec in open('seq/7BC6.fasta').read().split('>'):
    if rec.strip():
        q = rec.strip().split('\n')
        ent[q[0].split('|')[0]] = ''.join(q[1:])
chicken = ent['7BC6_1'][:473]
human = ''.join(l.strip() for l in open('seq/Q96NT5.fasta') if not l.startswith('>'))
if len(chicken) != 473 or len(human) != 459:
    die('sequence lengths %d and %d, expected 473 and 459' % (len(chicken), len(human)))

rows = [l.rstrip('\n').split('\t') for l in open(POSSET)]
pc = {k: i for i, k in enumerate(rows[0])}
SET_AT = {}
for r in rows[1:]:
    cp = int(r[pc['chicken_pos']])
    SET_AT[cp] = {'chicken_aa': r[pc['chicken_aa']], 'human_aa': r[pc['human_aa']],
                  'human_pos': int(r[pc['human_pos']]),
                  'differs': r[pc['differs']] == '1',
                  'parker': r[pc['in_parker_named']] == '1'}
PRIM = sorted(p for p, v in SET_AT.items() if v['differs'])

hu2ck = {int(k): v for k, v in json.load(open(MAP))['hu2ck'].items()}
LOOSE = set()
for k, v in json.load(open(PARKER))['residues'].items():
    if chicken[int(k) - 1] == human[v['human'] - 1]:
        LOOSE.add(int(k))
for r in [l.rstrip('\n').split('\t') for l in open(HFM)][1:]:
    cp = hu2ck.get(int(r[0]))
    if cp is None:
        die('hereditary folate malabsorption position %s has no chicken partner' % r[0])
    LOOSE.add(cp)
STRICT = set(SET_AT)
log('positions off limits under the annotated policy of codes/07: %d' % len(LOOSE))
log('positions off limits here, the whole pharmacophore set: %d' % len(STRICT))
log('primary humanising substitutions: %d' % len(PRIM))


def clean(positions):
    return all(p not in STRICT for p in positions)


def check(pos, frm):
    if chicken[pos - 1] != frm:
        die('a scan states %s at chicken %d, the sequence has %s'
            % (frm, pos, chicken[pos - 1]))


CAND = {}


def offer(target, label, positions, residues, p, source, background):
    for q, frm in zip(positions, residues[0]):
        check(q, frm)
    CAND.setdefault(target, []).append(
        {'label': label, 'positions': list(positions), 'to': list(residues[1]),
         'from_chicken': list(residues[0]), 'p': p, 'source': source,
         'background': background, 'clean': clean(positions),
         'loose_clean': all(q not in LOOSE and q not in PRIM for q in positions),
         'touches': [{'pos': q, 'in_set': q in SET_AT, 'primary': q in PRIM,
                      'parker': bool(SET_AT.get(q, {}).get('parker'))}
                     for q in positions]})


sweep = json.load(open(SWEEP))
RADIUS = sweep.get('published_radius', 3)
for t, blocks in sweep['per_target'].items():
    t = int(t)
    for r in blocks.get('published', []):
        if not r['crosses'] or r['separation'] > RADIUS:
            continue
        j = int(r['rescuer_pos'])
        offer(t, '%s%d%s' % (chicken[j - 1], j, r['rescuer_to']), [j],
              ([chicken[j - 1]], [r['rescuer_to']]), r['p_target'],
              'step 2 sweep, window of %d' % RADIUS, 'wild-type chicken')

d4 = json.load(open(D4))
for r in d4['top_doubles']:
    if not r['crosses']:
        continue
    offer(int(d4['target_chicken_pos']),
          '%s%d%s with %s%d%s' % (r['from_i'], r['pos_i'], r['to_i'],
                                  r['from_j'], r['pos_j'], r['to_j']),
          [r['pos_i'], r['pos_j']],
          ([r['from_i'], r['from_j']], [r['to_i'], r['to_j']]),
          r['p_target'], 'step 4 pair scan, coupled positions', 'wild-type chicken')

c5 = json.load(open(C5))
A = c5['part_a']
for r in A['top_doubles']:
    if not r['crosses']:
        continue
    offer(int(A['target']),
          '%s%d%s with %s%d%s' % (r['from_i'], r['pos_i'], r['to_i'],
                                  r['from_j'], r['pos_j'], r['to_j']),
          [r['pos_i'], r['pos_j']],
          ([r['from_i'], r['from_j']], [r['to_i'], r['to_j']]),
          r['p_target'], 'step 5 pair scan, nearest positions', 'wild-type chicken')

c6 = json.load(open(C6))
for t, blk in c6['per_target'].items():
    t = int(t)
    for r in blk.get('crossing_top', []) or []:
        j = r['pos']
        offer(t, '%s%d%s' % (chicken[j - 1], j, r['to']), [j],
              ([chicken[j - 1]], [r['to']]), r['p_construct'],
              'step 6 sweep, construct background', 'the construct')
s8 = json.load(open(S8))
for r in s8['top_doubles']:
    if not r['crosses']:
        continue
    offer(int(s8['target']),
          '%s%d%s with %s%d%s' % (r['from_construct_i'], r['pos_i'], r['to_i'],
                                  r['from_construct_j'], r['pos_j'], r['to_j']),
          [r['pos_i'], r['pos_j']],
          ([chicken[r['pos_i'] - 1], chicken[r['pos_j'] - 1]], [r['to_i'], r['to_j']]),
          r['p_target'], 'step 8 pair scan, construct background', 'the construct')

for policy, picks in c6['chosen'].items():
    for r in picks:
        if not r.get('found'):
            continue
        j = r['pos']
        offer(int(r['target']), '%s%d%s' % (chicken[j - 1], j, r['to']), [j],
              ([chicken[j - 1]], [r['to']]), r['p_construct'],
              'step 6 sweep, construct background', 'the construct')

TARGETS = sorted(CAND)
log('')
log('rescues on record, by target')
log('  target  total  crossing and clean  crossing and clean under the looser policy')
for t in TARGETS:
    c = CAND[t]
    log('  %6d %6d %19d %43d'
        % (t, len(c), sum(1 for r in c if r['clean']),
           sum(1 for r in c if r['loose_clean'])))

CHOSEN, unres = {}, []
log('')
log('the strongest rescue touching no position in the pharmacophore set')
log('  target  rescue                  p            source                               '
    'background')
overruled = []
for t in sorted(set(PRIM) & set(list(CAND) + PRIM)):
    pool = [r for r in CAND.get(t, []) if r['clean'] and r['p'] >= THRESHOLD]
    if not pool:
        if t in CAND:
            unres.append(t)
            log('  %6d  %s' % (t, 'nothing crossing that leaves the set untouched'))
        continue
    in_construct = [r for r in pool if r['background'] == 'the construct']
    b = max(in_construct or pool, key=lambda r: r['p'])
    by_p = max(pool, key=lambda r: r['p'])
    if b is not by_p:
        overruled.append({'target': t, 'chosen': b['label'], 'chosen_p': b['p'],
                          'chosen_background': b['background'],
                          'chosen_positions': list(b['positions']),
                          'chosen_to': list(b['to']),
                          'strongest': by_p['label'], 'strongest_p': by_p['p'],
                          'strongest_background': by_p['background'],
                          'strongest_positions': list(by_p['positions']),
                          'strongest_to': list(by_p['to']),
                          'positions_swapped': sorted(set(b['positions'])
                                                      ^ set(by_p['positions']))})
    CHOSEN[t] = b
    log('  %6d  %-22s %.6e  %-36s %s'
        % (t, b['label'], b['p'], b['source'], b['background']))
for q in CHOSEN.values():
    if q['background'] not in ('the construct', 'wild-type chicken'):
        die('a chosen rescue carries an unrecognised background, %s, so the preference '
            'between backgrounds cannot be applied' % q['background'])
log('')
if overruled:
    log('rescues where a construct-background measurement displaced a larger wild-type '
        'one, because a likelihood measured against wild-type chicken is not on the same '
        'scale as one measured in the background the residue will have')
    for r in overruled:
        log('  chicken %d: %s at %.6e in %s, over %s at %.6e in %s'
            % (r['target'], r['chosen'], r['chosen_p'], r['chosen_background'],
               r['strongest'], r['strongest_p'], r['strongest_background']))
else:
    log('no target had a construct-background rescue displace a larger wild-type one')

seen = {}
for t, b in CHOSEN.items():
    for q, y in zip(b['positions'], b['to']):
        if q in seen and seen[q] != y:
            die('two chosen rescuers want different residues at chicken %d, %s and %s'
                % (q, seen[q], y))
        seen[q] = y
log('')
log('the construct this gives: %d primaries plus %d rescuer positions, %d positions '
    'changed in all, and no position in the pharmacophore set overwritten'
    % (len(PRIM), len(seen), len(set(PRIM) | set(seen))))
log('rescued with the set untouched: %s' % sorted(CHOSEN))
log('still without an admissible rescue: %s' % (unres or 'none'))

seq = list(chicken)
for p in PRIM:
    seq[p - 1] = SET_AT[p]['human_aa']
for q, y in seen.items():
    seq[q - 1] = y
seq = ''.join(seq)
nd = sum(1 for a, b in zip(seq, chicken) if a != b)
if nd != len(set(PRIM) | set(seen)):
    die('the assembled sequence differs at %d positions, %d expected'
        % (nd, len(set(PRIM) | set(seen))))
for p in PRIM:
    if seq[p - 1] != SET_AT[p]['human_aa']:
        die('the assembled sequence is not humanised at chicken %d' % p)
log('assembled sequence checked: %d positions differ from chicken and every primary '
    'carries its human residue' % nd)

json.dump({'threshold': THRESHOLD,
           'off_limits_strict': sorted(STRICT), 'off_limits_loose': sorted(LOOSE),
           'primaries': PRIM,
           'chosen': {str(t): CHOSEN[t] for t in CHOSEN},
           'unrescued': unres,
           'choice_rule': 'among admissible rescues crossing the threshold, one measured '
                          'in the construct background is preferred to one measured '
                          'against wild-type chicken, and within a background the '
                          'largest likelihood wins',
           'overruled_by_background': overruled,
           'rescuer_positions': {str(q): y for q, y in sorted(seen.items())},
           'n_positions_changed': nd,
           'sequence': seq,
           'candidates': {str(t): CAND[t] for t in CAND}},
          open(OUT, 'w'), indent=1)
log('')
log('wrote ' + OUT)
with open(LOG, 'w') as f:
    f.write('\n'.join(lines) + '\n')
