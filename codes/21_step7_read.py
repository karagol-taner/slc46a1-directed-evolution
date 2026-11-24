import glob
import json
import os
import statistics
import sys

POSSET = 'data/pharmacophore_position_set.tsv'
RISK = 'data/primary_substitution_risk.tsv'
C5 = 'data/construct_report.json'
ADM = 'data/admissible_construct.json'
CG = 'data/coupling_guided_rescue.json'
ST = 'data/structural_candidates.json'
SCAN = 'data/esm2_step2/rescue_scan.tsv'
OUT = 'data/step7_report.json'
LOG = 'data/21_step7_read_log.txt'
THRESHOLD = 1e-3
TARGET_B = 322
AA = 'ACDEFGHIKLMNPQRSTVWY'
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
    c = [d for d in sorted(glob.glob('data/esm2_step7*')) if os.path.isdir(d)]
    if not c:
        die('no esm2_step7 folder under data/. Run '
            'codes/SLC46A1_step7_admissible_and_322.ipynb on Colab and copy esm2_step7/ '
            'from Drive into data/, or pass its path.')
    if len(c) > 1:
        die('several candidates under data/, so name the one to read: %s' % c)
    return c[0]


SRC = pick()
for p in (POSSET, RISK, C5, ADM, CG, ST, SCAN):
    if not os.path.exists(p):
        die(p + ' not found')
need = ['admissible_profile.tsv', 'b322_double_rescue.tsv', 'b322_singles.tsv']
missing = [f for f in need if not os.path.exists(os.path.join(SRC, f))]
if missing:
    die('%s is missing %s, so the run has not finished its last cell' % (SRC, missing))

ent = {}
for rec in open('seq/7BC6.fasta').read().split('>'):
    if rec.strip():
        q = rec.strip().split('\n')
        ent[q[0].split('|')[0]] = ''.join(q[1:])
chicken = ent['7BC6_1'][:473]
human = ''.join(l.strip() for l in open('seq/Q96NT5.fasta') if not l.startswith('>'))
if len(chicken) != 473 or len(human) != 459:
    die('sequence lengths %d and %d' % (len(chicken), len(human)))

rows = [l.rstrip('\n').split('\t') for l in open(POSSET)]
pc = {k: i for i, k in enumerate(rows[0])}
PSET = {}
for r in rows[1:]:
    cp = int(r[pc['chicken_pos']])
    PSET[cp] = {'chicken_aa': r[pc['chicken_aa']], 'human_aa': r[pc['human_aa']],
                'human_pos': int(r[pc['human_pos']]),
                'differs': r[pc['differs']] == '1'}
PRIM = sorted(p for p, v in PSET.items() if v['differs'])
adm = json.load(open(ADM))
RESC = {int(k): v for k, v in adm['rescuer_positions'].items()}
inset = [q for q in RESC if q in PSET]
if inset:
    die('the construct claims to spare the pharmacophore set but writes at %s' % inset)
c5 = json.load(open(C5))
PROF5 = {x['chicken_pos']: x for x in c5['part_b']['profile']}
rr = [l.rstrip('\n').split('\t') for l in open(RISK)]
rcl = {k: i for i, k in enumerate(rr[0])}
WT = {int(r[rcl['chicken_pos']]): float(r[rcl['p_human_wt_background']]) for r in rr[1:]}


def build(primaries_in, rescuers_in, leave_out=None):
    s = list(chicken)
    if primaries_in:
        for p in PRIM:
            if p != leave_out:
                s[p - 1] = PSET[p]['human_aa']
    if rescuers_in:
        for q, y in RESC.items():
            if q != leave_out:
                s[q - 1] = y
    return ''.join(s)


if build(True, True) != adm['sequence']:
    die('the construct rebuilt here differs from the sequence stored in ' + ADM)
log('result             : %s' % SRC)
log('construct          : %d primaries plus %d rescuer positions %s, %d changes, and the '
    'rebuild matches the stored sequence'
    % (len(PRIM), len(RESC), sorted(RESC), adm['n_positions_changed']))

arows = [l.rstrip('\n').split('\t')
         for l in open(os.path.join(SRC, 'admissible_profile.tsv'))]
ah = {k: i for i, k in enumerate(arows[0])}
for k in ('chicken_pos', 'p_wt', 'p_primaries', 'p_admissible', 'p_step5_final'):
    if k not in ah:
        die('the profile lacks the %s column; it has %s' % (k, arows[0]))
prof, worst_wt, worst_pr, worst_5 = [], (0.0, None), (0.0, None), (0.0, None)
for r in arows[1:]:
    t = int(r[ah['chicken_pos']])
    if t not in PRIM:
        die('the profile holds chicken %d, which is not a primary' % t)
    a = float(r[ah['p_wt']])
    rel = abs(a - WT[t]) / max(a, WT[t], 1e-30)
    if rel > worst_wt[0]:
        worst_wt = (rel, (t, a, WT[t]))
    a = float(r[ah['p_primaries']])
    rel = abs(a - PROF5[t]['primaries']) / max(a, PROF5[t]['primaries'], 1e-30)
    if rel > worst_pr[0]:
        worst_pr = (rel, (t, a, PROF5[t]['primaries']))
    a = float(r[ah['p_step5_final']])
    b5 = PROF5[t]['primaries_rescuers_primary_wins']
    rel = abs(a - b5) / max(a, b5, 1e-30)
    if rel > worst_5[0]:
        worst_5 = (rel, (t, a, b5))
    prof.append({'chicken_pos': t, 'human_aa': PSET[t]['human_aa'],
                 'p_wt': float(r[ah['p_wt']]),
                 'p_primaries': float(r[ah['p_primaries']]),
                 'p_admissible': float(r[ah['p_admissible']]),
                 'p_step5_final': float(r[ah['p_step5_final']]),
                 'fold_over_wt': float(r[ah['p_admissible']]) / float(r[ah['p_wt']]),
                 'fold_over_step5': (float(r[ah['p_admissible']])
                                     / float(r[ah['p_step5_final']]))})
if sorted(x['chicken_pos'] for x in prof) != PRIM:
    die('the profile covers %d primaries, the set has %d' % (len(prof), len(PRIM)))
for name, w, src in (('wild type against step 1', worst_wt, RISK),
                     ('primaries against step 5', worst_pr, C5),
                     ('step 5 column against step 5', worst_5, C5)):
    if w[0] > TOL:
        t, a, b = w[1]
        die('the %s column disagrees at chicken %d, %.6e here and %.6e in %s'
            % (name, t, a, b, src))
prof.sort(key=lambda x: x['chicken_pos'])
log('part A scale checks: wild type against step 1 worst %.2e, primaries against step 5 '
    'worst %.2e, step 5 column worst %.2e' % (worst_wt[0], worst_pr[0], worst_5[0]))
log('')
log('chicken  human   wild type   primaries   admissible   step 5 final  fold over wt  '
    'fold over step 5')
for x in prof:
    log('%7d %6s %11.3e %11.3e %12.3e %14.3e %13.2f %17.2f'
        % (x['chicken_pos'], x['human_aa'], x['p_wt'], x['p_primaries'],
           x['p_admissible'], x['p_step5_final'],
           x['fold_over_wt'], x['fold_over_step5']))
belowA = [x['chicken_pos'] for x in prof if x['p_admissible'] < THRESHOLD]
below5 = [x['chicken_pos'] for x in prof if x['p_step5_final'] < THRESHOLD]
log('')
log('admissible construct below %.0e at %d of %d primaries%s'
    % (THRESHOLD, len(belowA), len(prof), (': %s' % belowA) if belowA else ''))
log('step 5 construct  below %.0e at %d of %d primaries%s'
    % (THRESHOLD, len(below5), len(prof), (': %s' % below5) if below5 else ''))
log('median fold of the admissible construct over the step 5 one: %.2f'
    % statistics.median([x['fold_over_step5'] for x in prof]))
HOLD = []
for t in sorted(int(k) for k in adm['chosen']):
    x = [y for y in prof if y['chicken_pos'] == t][0]
    HOLD.append({'target': t, 'label': adm['chosen'][str(t)]['label'],
                 'chosen_at': adm['chosen'][str(t)]['p'],
                 'in_construct': x['p_admissible'],
                 'holds': x['p_admissible'] >= THRESHOLD,
                 'fold': x['p_admissible'] / adm['chosen'][str(t)]['p']})
    log('  rescue at chicken %d, %s: chosen at %.6e, here %.6e, %s'
        % (t, HOLD[-1]['label'], HOLD[-1]['chosen_at'], HOLD[-1]['in_construct'],
           'holds' if HOLD[-1]['holds'] else 'does not hold'))

COUP = [int(r['chicken_pos'])
        for r in json.load(open(CG))['targets'][str(TARGET_B)]['top_coupled']]
NEAR = [int(c['chicken_pos'])
        for c in json.load(open(ST))['targets'][str(TARGET_B)]['candidates']]
CAND = sorted(set(COUP) | set(NEAR))
EXP_PAIRS = len(CAND) * (len(CAND) - 1) // 2
EXP_DOUBLES = EXP_PAIRS * 19 * 19
EXP_SINGLES = len(CAND) * 19
BGB = build(True, False, TARGET_B)
BASE_B = PROF5[TARGET_B]['primaries']

srows = [l.rstrip('\n').split('\t') for l in open(os.path.join(SRC, 'b322_singles.tsv'))]
shh = {k: i for i, k in enumerate(srows[0])}
sing = {}
for r in srows[1:]:
    j, y = int(r[shh['pos']]), r[shh['to']]
    if r[shh['from_construct']] != BGB[j - 1]:
        die('the singles file says the construct has %s at chicken %d, the rebuild has %s'
            % (r[shh['from_construct']], j, BGB[j - 1]))
    sing[(j, y)] = float(r[shh['p_construct']])
if len(sing) != EXP_SINGLES:
    die('%d singles read, %d expected over %d candidate positions'
        % (len(sing), EXP_SINGLES, len(CAND)))

drows = [l.rstrip('\n').split('\t')
         for l in open(os.path.join(SRC, 'b322_double_rescue.tsv'))]
dh = {k: i for i, k in enumerate(drows[0])}
D, seen = [], set()
for r in drows[1:]:
    i, x = int(r[dh['pos_i']]), r[dh['to_i']]
    j, y = int(r[dh['pos_j']]), r[dh['to_j']]
    if i not in CAND or j not in CAND:
        die('the doubles file holds a pair outside the candidate union, %d with %d'
            % (i, j))
    if r[dh['from_construct_i']] != BGB[i - 1] or r[dh['from_construct_j']] != BGB[j - 1]:
        die('the doubles file disagrees with the rebuilt background at %d or %d' % (i, j))
    if (i, x, j, y) in seen:
        die('the doubles file repeats %s%d%s with %s%d%s' % (BGB[i - 1], i, x,
                                                             BGB[j - 1], j, y))
    seen.add((i, x, j, y))
    if (i, x) not in sing or (j, y) not in sing:
        die('a half of %s%d%s with %s%d%s has no single score'
            % (BGB[i - 1], i, x, BGB[j - 1], j, y))
    p = float(r[dh['p_target']])
    pi, pj = sing[(i, x)], sing[(j, y)]
    D.append({'pos_i': i, 'to_i': x, 'from_construct_i': BGB[i - 1],
              'pos_j': j, 'to_j': y, 'from_construct_j': BGB[j - 1],
              'p_target': p, 'best_single_of_pair': max(pi, pj),
              'fold_vs_baseline': p / BASE_B,
              'fold_vs_best_single_of_pair': p / max(pi, pj),
              'crosses': p >= THRESHOLD,
              'admissible': (i not in PSET) and (j not in PSET)})
if len(D) != EXP_DOUBLES:
    die('%d doubles read, %d expected over %d pairs' % (len(D), EXP_DOUBLES, EXP_PAIRS))
D.sort(key=lambda r: -r['p_target'])
cross = [r for r in D if r['crosses']]
adm_cross = [r for r in cross if r['admissible']]
folds = [r['fold_vs_best_single_of_pair'] for r in D]
best_single = max(sing.items(), key=lambda kv: kv[1])
log('')
log('part B, chicken %d by pairs in the construct background' % TARGET_B)
log('  candidate union  : %d positions, %d coupled and %d nearest, %d in common'
    % (len(CAND), len(COUP), len(NEAR), len(set(COUP) & set(NEAR))))
log('  baseline         : %.6e from the step 5 primaries background' % BASE_B)
log('  singles          : %d, strongest %s%d%s at %.6e'
    % (len(sing), BGB[best_single[0][0] - 1], best_single[0][0], best_single[0][1],
       best_single[1]))
log('  doubles          : %d over %d pairs, %d crossing %.0e, %d of those leaving the '
    'set untouched' % (len(D), EXP_PAIRS, len(cross), THRESHOLD, len(adm_cross)))
log('  fold over the better half: median %.2f, max %.2f'
    % (statistics.median(folds), max(folds)))
log('  the five strongest')
for r in D[:5]:
    log('    %s%d%s with %s%d%s at %.6e, %.1f-fold over baseline, %.2f-fold over the '
        'better half, spares the set %s'
        % (r['from_construct_i'], r['pos_i'], r['to_i'], r['from_construct_j'],
           r['pos_j'], r['to_j'], r['p_target'], r['fold_vs_baseline'],
           r['fold_vs_best_single_of_pair'], r['admissible']))
if adm_cross:
    r = adm_cross[0]
    log('  strongest crossing pair that spares the set: %s%d%s with %s%d%s at %.6e'
        % (r['from_construct_i'], r['pos_i'], r['to_i'], r['from_construct_j'],
           r['pos_j'], r['to_j'], r['p_target']))
log('  verdict          : %s'
    % ('chicken %d is rescued by a pair in the construct background, %d crossing and %d '
       'of those sparing the set' % (TARGET_B, len(cross), len(adm_cross)) if cross else
       'no pair at these %d positions lifts chicken %d above %.0e in the construct '
       'background, and the measurement bounds those pairs and not the whole double space'
       % (len(CAND), TARGET_B, THRESHOLD)))

json.dump({'source': SRC, 'threshold': THRESHOLD,
           'construct': {'primaries': PRIM,
                         'rescuer_positions': {str(k): v for k, v in RESC.items()},
                         'n_changes': adm['n_positions_changed'],
                         'rebuild_matches_stored': True},
           'part_a': {'wt_vs_step1_worst': worst_wt[0],
                      'primaries_vs_step5_worst': worst_pr[0],
                      'step5_column_worst': worst_5[0],
                      'below_threshold_admissible': belowA,
                      'below_threshold_step5': below5,
                      'median_fold_over_step5': statistics.median(
                          [x['fold_over_step5'] for x in prof]),
                      'rescue_survival': HOLD, 'profile': prof},
           'part_b': {'target': TARGET_B, 'coupled': COUP, 'nearest': NEAR,
                      'candidates': CAND, 'baseline': BASE_B,
                      'n_pairs': EXP_PAIRS, 'n_doubles': len(D), 'n_singles': len(sing),
                      'best_single': {'pos': best_single[0][0],
                                      'to': best_single[0][1], 'p': best_single[1]},
                      'n_crossing': len(cross),
                      'n_crossing_sparing_the_set': len(adm_cross),
                      'best_double': D[0],
                      'best_sparing_double': adm_cross[0] if adm_cross else None,
                      'top_doubles': D[:20],
                      'fold_over_better_half': {'median': statistics.median(folds),
                                                'max': max(folds)}}},
          open(OUT, 'w'), indent=1)
log('')
log('wrote ' + OUT)
with open(LOG, 'w') as f:
    f.write('\n'.join(lines) + '\n')
