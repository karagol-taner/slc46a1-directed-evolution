import glob
import json
import os
import statistics
import sys

POSSET = 'data/pharmacophore_position_set.tsv'
ADM = 'data/admissible_construct.json'
RISK = 'data/primary_substitution_risk.tsv'
C5 = 'data/construct_report.json'
S7 = 'data/esm2_step7/step7_summary.json'
S8 = 'data/step8_report.json'
SEQF = 'seq/7BC6.fasta'
OUT = 'data/step9_report.json'
LOG = 'data/23_step9_read_log.txt'
THRESHOLD = 1e-3
TOL = 5e-2

lines = []


def log(m):
    print(m, flush=True)
    lines.append(str(m))


def die(m):
    with open(LOG, 'w') as f:
        f.write('\n'.join(lines + ['FATAL ' + m]) + '\n')
    sys.exit('FATAL ' + m)


def pick():
    if len(sys.argv) > 1:
        return sys.argv[1]
    c = [d for d in sorted(glob.glob('data/esm2_step9*')) if os.path.isdir(d)]
    if not c:
        die('no esm2_step9 folder under data/. Run '
            'codes/SLC46A1_step9_revised_construct.ipynb on Colab, then copy '
            'data/esm2_step9/ from Drive into the local data/.')
    return c[-1]


SRC = pick()
for f in ('revised_profile.tsv', 'step9_summary.json'):
    if not os.path.exists(os.path.join(SRC, f)):
        die('%s holds no %s, so the run did not finish' % (SRC, f))
for f in (POSSET, ADM, RISK, C5, S7, S8, SEQF):
    if not os.path.exists(f):
        die('missing input %s' % f)

summ = json.load(open(os.path.join(SRC, 'step9_summary.json')))
if abs(summ.get('threshold', 0) - THRESHOLD) > 0:
    die('the run used a threshold of %s, this reader expects %.0e'
        % (summ.get('threshold'), THRESHOLD))

rows = [l.rstrip('\n').split('\t') for l in open(POSSET)]
pc = {k: i for i, k in enumerate(rows[0])}
PSET = {int(r[pc['chicken_pos']]): {'chicken_aa': r[pc['chicken_aa']],
                                    'human_aa': r[pc['human_aa']],
                                    'human_pos': int(r[pc['human_pos']]),
                                    'differs': r[pc['differs']] == '1'}
        for r in rows[1:]}
PRIM = sorted(p for p, v in PSET.items() if v['differs'])

ent = {}
for rec in open(SEQF).read().split('>'):
    if rec.strip():
        p = rec.strip().split('\n')
        ent[p[0].split('|')[0]] = ''.join(p[1:])
if '7BC6_1' not in ent:
    die('%s holds no 7BC6_1 record' % SEQF)
chicken = ent['7BC6_1'][:473]
for p in PRIM:
    if chicken[p - 1] != PSET[p]['chicken_aa']:
        die('the chicken sequence has %s at chicken %d, the position set says %s'
            % (chicken[p - 1], p, PSET[p]['chicken_aa']))

adm = json.load(open(ADM))
if sorted(adm['primaries']) != PRIM:
    die('the construct file lists different primaries from the position set')
RESC = {int(k): v for k, v in adm['rescuer_positions'].items()}
built = list(chicken)
for p in PRIM:
    built[p - 1] = PSET[p]['human_aa']
for q, y in RESC.items():
    if q in PSET:
        die('rescuer position %d is inside the pharmacophore set' % q)
    built[q - 1] = y
built = ''.join(built)
if built != adm['sequence']:
    die('the construct file sequence does not match the one rebuilt from its own parts')
if summ['construct']['sequence'] != built:
    n = sum(1 for a, b in zip(summ['construct']['sequence'], built) if a != b)
    die('the run scored a sequence differing from the rebuilt construct at %d positions'
        % n)
if summ['construct']['n_positions_changed'] != adm['n_positions_changed']:
    die('the run records %d changed positions, the construct file %d'
        % (summ['construct']['n_positions_changed'], adm['n_positions_changed']))
log('construct rebuilt from %s and matched against the run, %d primaries and %d rescuer '
    'positions %s, %d changes from chicken'
    % (ADM, len(PRIM), len(RESC), sorted(RESC), adm['n_positions_changed']))

DISP = {int(d['target']): d for d in adm.get('overruled_by_background', [])}
if not DISP:
    die('the construct file records no displaced rescue, so this is not a revised '
        'construct and there was nothing for step 9 to re-score')
SWAP = sorted({q for d in DISP.values() for q in d['positions_swapped']})
log('rescuer positions changed since the earlier set-level run: %s' % SWAP)

prows = [l.rstrip('\n').split('\t') for l in open(os.path.join(SRC, 'revised_profile.tsv'))]
ph = {k: i for i, k in enumerate(prows[0])}
P = {}
for r in prows[1:]:
    t = int(r[ph['chicken_pos']])
    if r[ph['human_aa']] != PSET[t]['human_aa']:
        die('the profile says the human residue at chicken %d is %s, the position set '
            'says %s' % (t, r[ph['human_aa']], PSET[t]['human_aa']))
    P[t] = {k: float(r[ph[k]]) for k in ('p_wt', 'p_primaries', 'p_revised', 'p_step7')}
if sorted(P) != PRIM:
    die('the profile covers %d positions, the %d primaries expected: %s'
        % (len(P), len(PRIM), sorted(set(PRIM) ^ set(P))))

rrows = [l.rstrip('\n').split('\t') for l in open(RISK)]
rc = {k: i for i, k in enumerate(rrows[0])}
WT = {int(r[rc['chicken_pos']]): float(r[rc['p_human_wt_background']]) for r in rrows[1:]}
PR5 = {x['chicken_pos']: x for x in json.load(open(C5))['part_b']['profile']}
worst_wt = max(((abs(P[t]['p_wt'] - WT[t]) / max(P[t]['p_wt'], WT[t], 1e-30)), t)
               for t in PRIM)
worst_pr = max(((abs(P[t]['p_primaries'] - PR5[t]['primaries'])
                 / max(P[t]['p_primaries'], PR5[t]['primaries'], 1e-30)), t)
               for t in PRIM)
if worst_wt[0] > TOL:
    die('the wild-type column does not reproduce step 1 at chicken %d, %.6e against %.6e'
        % (worst_wt[1], P[worst_wt[1]]['p_wt'], WT[worst_wt[1]]))
if worst_pr[0] > TOL:
    die('the primaries column does not reproduce step 5 at chicken %d, %.6e against %.6e'
        % (worst_pr[1], P[worst_pr[1]]['p_primaries'], PR5[worst_pr[1]]['primaries']))
log('scale checks: wild type against step 1 worst %.2e, primaries against step 5 worst '
    '%.2e' % (worst_wt[0], worst_pr[0]))

s7 = json.load(open(S7))
S7P = {x['chicken_pos']: x['p_admissible'] for x in s7['part_a']['profile']}
for t in PRIM:
    if abs(P[t]['p_step7'] - S7P[t]) / max(P[t]['p_step7'], S7P[t], 1e-30) > 1e-9:
        die('the profile quotes %.6e for chicken %d in the earlier construct, the earlier '
            'run recorded %.6e' % (P[t]['p_step7'], t, S7P[t]))

s8 = json.load(open(S8))
T8 = int(s8['target'])
B8 = s8['best_sparing_double']
if B8 is None:
    die('step 8 recorded no admissible crossing pair, so the revised construct cannot '
        'have come from it')
rel8 = abs(P[T8]['p_revised'] - B8['p_target']) / max(P[T8]['p_revised'],
                                                      B8['p_target'], 1e-30)
if rel8 > TOL:
    die('chicken %d scores %.6e in this construct and %.6e in step 8, a relative '
        'difference of %.2e, so the construct scored here is not the one step 8 chose'
        % (T8, P[T8]['p_revised'], B8['p_target'], rel8))
log('step 8 tie: chicken %d at %.6e here against %.6e there, relative difference %.2e'
    % (T8, P[T8]['p_revised'], B8['p_target'], rel8))

DROPPED, ADDED = [], []
for d in DISP.values():
    DROPPED += [q for q in d['strongest_positions'] if q not in d['chosen_positions']]
    ADDED += [q for q in d['chosen_positions'] if q not in d['strongest_positions']]
DROPPED, ADDED = sorted(set(DROPPED)), sorted(set(ADDED))
if not DROPPED:
    die('the displaced rescue vacated no position, so the construct did not change where '
        'the earlier set-level run measured it')
dist = {t: min(abs(t - q) for q in DROPPED) for t in PRIM}
by_dist = sorted(PRIM, key=lambda t: dist[t])
by_fold = sorted(PRIM, key=lambda t: -P[t]['p_revised'] / P[t]['p_step7'])
NEAR = len(DROPPED) + 3
near, far = by_dist[:NEAR], by_dist[NEAR:]
gap = dist[by_dist[NEAR]] - dist[by_dist[NEAR - 1]]
nearest_are_largest = set(near) == set(by_fold[:NEAR])
log('')
log('the swap dropped a rescuer at chicken %s and wrote one at chicken %s'
    % (DROPPED, ADDED))
log('  the %d primaries nearest the vacated position, chicken %s, lie within %d residues '
    'of it, and the next is %d away' % (NEAR, sorted(near), dist[by_dist[NEAR - 1]],
                                        dist[by_dist[NEAR]]))
log('  their folds over the earlier construct run %.2f to %.2f, the other %d run %.2f to '
    '%.2f' % (min(P[t]['p_revised'] / P[t]['p_step7'] for t in near),
              max(P[t]['p_revised'] / P[t]['p_step7'] for t in near), len(far),
              min(P[t]['p_revised'] / P[t]['p_step7'] for t in far),
              max(P[t]['p_revised'] / P[t]['p_step7'] for t in far)))
log('  the nearest %d are exactly the %d largest risers: %s'
    % (NEAR, NEAR, nearest_are_largest))

below = [t for t in PRIM if P[t]['p_revised'] < THRESHOLD]
below7 = [t for t in PRIM if P[t]['p_step7'] < THRESHOLD]
folds = [P[t]['p_revised'] / P[t]['p_step7'] for t in PRIM]
surv = []
for t in sorted(int(k) for k in adm['chosen']):
    ch = adm['chosen'][str(t)]
    surv.append({'target': t, 'label': ch['label'], 'chosen_at': ch['p'],
                 'chosen_background': ch['background'],
                 'in_construct': P[t]['p_revised'],
                 'fold': P[t]['p_revised'] / ch['p'],
                 'holds': P[t]['p_revised'] >= THRESHOLD})
hold = [h['target'] for h in surv if h['holds']]
lost = [h['target'] for h in surv if not h['holds']]

log('')
log('the revised construct, every primary with the other %d and all %d rescuers in place'
    % (len(PRIM) - 1, len(RESC)))
log('chicken  human    wild type   primaries     revised      earlier   fold over earlier')
for t in PRIM:
    log('%7d %6s %12.3e %11.3e %11.3e %11.3e %18.2f'
        % (t, PSET[t]['human_aa'], P[t]['p_wt'], P[t]['p_primaries'],
           P[t]['p_revised'], P[t]['p_step7'], P[t]['p_revised'] / P[t]['p_step7']))
log('')
log('  below %.0e : %d of %d primaries%s'
    % (THRESHOLD, len(below), len(PRIM), (', chicken %s' % below) if below else ''))
log('  the earlier construct was below it at %d%s'
    % (len(below7), (', chicken %s' % below7) if below7 else ''))
log('  fold over the earlier construct: median %.2f, lowest %.2f, highest %.2f'
    % (statistics.median(folds), min(folds), max(folds)))
log('')
log('does each rescue hold with the others present')
for h in surv:
    log('  chicken %d: %-20s chosen at %.6e in %-18s here %.6e, %s'
        % (h['target'], h['label'], h['chosen_at'], h['chosen_background'],
           h['in_construct'], 'holds' if h['holds'] else 'does not hold'))
log('  holding %s of %s, not holding %s'
    % (hold, [h['target'] for h in surv], lost or 'none'))
log('  still without an admissible rescue: %s' % (adm['unrescued'] or 'none'))
log('  verdict         : %d of five high-risk substitutions rescued under the strict '
    'policy, %d primaries below the threshold'
    % (len(hold), len(below)))

json.dump({'source': SRC, 'threshold': THRESHOLD,
           'construct': {'primaries': PRIM,
                         'rescuer_positions': {str(k): v for k, v in RESC.items()},
                         'n_positions_changed': adm['n_positions_changed'],
                         'rebuild_matches_run': True},
           'rescuer_positions_swapped': SWAP,
           'displaced_rescues': adm['overruled_by_background'],
           'checks': {'wt_vs_step1_worst': worst_wt[0],
                      'primaries_vs_step5_worst': worst_pr[0],
                      'step8_target': T8, 'step8_here': P[T8]['p_revised'],
                      'step8_there': B8['p_target'],
                      'step8_relative_difference': rel8},
           'profile': [dict(chicken_pos=t, human_aa=PSET[t]['human_aa'],
                            human_pos=PSET[t]['human_pos'], **P[t],
                            fold_over_step7=P[t]['p_revised'] / P[t]['p_step7'])
                       for t in PRIM],
           'below_threshold': below, 'below_threshold_step7': below7,
           'fold_over_step7': {'median': statistics.median(folds),
                               'min': min(folds), 'max': max(folds)},
           'neighbourhood_of_the_swap': {
               'dropped': DROPPED, 'added': ADDED,
               'n_near': NEAR,
               'near': sorted(near), 'near_max_distance': dist[by_dist[NEAR - 1]],
               'next_distance': dist[by_dist[NEAR]], 'distance_gap': gap,
               'near_fold_min': min(P[t]['p_revised'] / P[t]['p_step7'] for t in near),
               'near_fold_max': max(P[t]['p_revised'] / P[t]['p_step7'] for t in near),
               'far_fold_min': min(P[t]['p_revised'] / P[t]['p_step7'] for t in far),
               'far_fold_max': max(P[t]['p_revised'] / P[t]['p_step7'] for t in far),
               'nearest_are_the_largest_risers': nearest_are_largest},
           'rescue_survival': surv,
           'rescues_holding': hold, 'rescues_lost': lost,
           'unrescued': adm['unrescued']},
          open(OUT, 'w'), indent=1)
log('')
log('wrote ' + OUT)
with open(LOG, 'w') as f:
    f.write('\n'.join(lines) + '\n')
