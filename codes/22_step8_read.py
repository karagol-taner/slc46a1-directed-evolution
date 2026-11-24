import glob
import json
import os
import statistics
import sys

POSSET = 'data/pharmacophore_position_set.tsv'
ADM = 'data/admissible_construct.json'
CG = 'data/coupling_guided_rescue.json'
ST = 'data/structural_candidates.json'
SEQ = 'seq/7BC6.fasta'
S7 = 'data/esm2_step7/step7_summary.json'
OUT = 'data/step8_report.json'
LOG = 'data/22_step8_read_log.txt'
THRESHOLD = 1e-3
TARGET = 407
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
    c = [d for d in sorted(glob.glob('data/esm2_step8*')) if os.path.isdir(d)]
    if not c:
        die('no esm2_step8 folder under data/. Run '
            'codes/SLC46A1_step8_pairs_at_407.ipynb on Colab, then copy '
            'data/esm2_step8/ from Drive into the local data/.')
    return c[-1]


SRC = pick()
for f in ('c407_double_rescue.tsv', 'c407_singles.tsv', 'step8_summary.json'):
    if not os.path.exists(os.path.join(SRC, f)):
        die('%s holds no %s, so the run did not finish' % (SRC, f))
for f in (POSSET, ADM, CG, ST, S7, SEQ):
    if not os.path.exists(f):
        die('missing input %s' % f)

summ = json.load(open(os.path.join(SRC, 'step8_summary.json')))
if summ.get('target') != TARGET:
    die('the run scored chicken %s, this reader expects %d' % (summ.get('target'), TARGET))
if abs(summ.get('threshold', 0) - THRESHOLD) > 0:
    die('the run used a threshold of %s, this reader expects %.0e'
        % (summ.get('threshold'), THRESHOLD))

rows = [l.rstrip('\n').split('\t') for l in open(POSSET)]
pc = {k: i for i, k in enumerate(rows[0])}
PSET = {}
for r in rows[1:]:
    PSET[int(r[pc['chicken_pos']])] = {'chicken_aa': r[pc['chicken_aa']],
                                       'human_aa': r[pc['human_aa']],
                                       'human_pos': int(r[pc['human_pos']]),
                                       'differs': r[pc['differs']] == '1'}
PRIM = sorted(p for p, v in PSET.items() if v['differs'])
adm = json.load(open(ADM))
if sorted(adm['primaries']) != PRIM:
    die('the admissible construct lists different primaries from the position set')
RESC = {int(k): v for k, v in adm['rescuer_positions'].items()}
if str(TARGET) not in adm['chosen']:
    die('the admissible construct carries no rescue at chicken %d' % TARGET)
CH = adm['chosen'][str(TARGET)]
R_T = [int(p) for p in CH['positions']]
R_TO = list(CH['to'])
KEEP = {q: y for q, y in RESC.items() if q not in R_T}

full = adm['sequence']
if len(full) != 473:
    die('the stored construct sequence is %d residues, expected 473' % len(full))
ent = {}
for rec in open(SEQ).read().split('>'):
    if rec.strip():
        parts = rec.strip().split('\n')
        ent[parts[0].split('|')[0]] = ''.join(parts[1:])
if '7BC6_1' not in ent:
    die('%s holds no 7BC6_1 record' % SEQ)
chicken = ent['7BC6_1'][:473]
if len(chicken) != 473:
    die('the chicken sequence is %d residues, expected 473' % len(chicken))
for p in PRIM:
    if chicken[p - 1] != PSET[p]['chicken_aa']:
        die('the chicken sequence has %s at chicken %d, the position set says %s'
            % (chicken[p - 1], p, PSET[p]['chicken_aa']))
    if full[p - 1] != PSET[p]['human_aa']:
        die('the stored construct is not humanised at chicken %d' % p)
for q, y in RESC.items():
    if full[q - 1] != y:
        die('the stored construct does not carry the rescuer %s at %d' % (y, q))
    if chicken[q - 1] == y:
        die('the rescuer at %d is not a substitution' % q)

BG = list(chicken)
for p in PRIM:
    if p != TARGET:
        BG[p - 1] = PSET[p]['human_aa']
for q, y in KEEP.items():
    BG[q - 1] = y
BG = ''.join(BG)
if summ['background']['sequence'] != BG:
    n = sum(1 for a, b in zip(summ['background']['sequence'], BG) if a != b)
    die('the background the run recorded differs from the rebuilt one at %d positions' % n)
rebuilt = list(BG)
for q, y in zip(R_T, R_TO):
    rebuilt[q - 1] = y
reference = list(full)
reference[TARGET - 1] = PSET[TARGET]['chicken_aa']
if ''.join(rebuilt) != ''.join(reference):
    die('the background plus %s does not rebuild the stored construct with chicken %d '
        'left alone' % (CH['label'], TARGET))
log('background rebuilt from %s and checked against the run, %d positions differ from '
    'chicken' % (ADM, sum(1 for a, b in zip(BG, chicken) if a != b)))
log('  %d primaries written in, rescuers held at %s, released at %s'
    % (len(PRIM) - 1, sorted(KEEP), R_T))

COUP = [int(r['chicken_pos'])
        for r in json.load(open(CG))['targets'][str(TARGET)]['top_coupled']]
NEAR = [int(c['chicken_pos'])
        for c in json.load(open(ST))['targets'][str(TARGET)]['candidates']]
CAND = sorted(set(COUP) | set(NEAR))
if CAND != sorted(summ['candidates']):
    die('the run scored %d candidate positions, the rankings give %d: %s against %s'
        % (len(summ['candidates']), len(CAND), sorted(summ['candidates']), CAND))
EXP_PAIRS = len(CAND) * (len(CAND) - 1) // 2
EXP_DOUBLES = EXP_PAIRS * 19 * 19
EXP_SINGLES = len(CAND) * 19
IN_SET = [p for p in CAND if p in PSET]
log('candidates %d, %d coupled and %d nearest, %d in both, %d inside the position set %s'
    % (len(CAND), len(COUP), len(NEAR), len(set(COUP) & set(NEAR)), len(IN_SET), IN_SET))

srows = [l.rstrip('\n').split('\t') for l in open(os.path.join(SRC, 'c407_singles.tsv'))]
sh = {k: i for i, k in enumerate(srows[0])}
sing, prev = {}, {}
for r in srows[1:]:
    j, y = int(r[sh['pos']]), r[sh['to']]
    if r[sh['from_construct']] != BG[j - 1]:
        die('the singles file says the background has %s at chicken %d, the rebuild has %s'
            % (r[sh['from_construct']], j, BG[j - 1]))
    sing[(j, y)] = float(r[sh['p_construct']])
    if r[sh['p_wt_step2']]:
        prev[(j, y)] = float(r[sh['p_wt_step2']])
if len(sing) != EXP_SINGLES:
    die('%d singles read, %d expected over %d candidate positions'
        % (len(sing), EXP_SINGLES, len(CAND)))

drows = [l.rstrip('\n').split('\t')
         for l in open(os.path.join(SRC, 'c407_double_rescue.tsv'))]
dh = {k: i for i, k in enumerate(drows[0])}
D, seen = [], set()
for r in drows[1:]:
    i, x = int(r[dh['pos_i']]), r[dh['to_i']]
    j, y = int(r[dh['pos_j']]), r[dh['to_j']]
    if i not in CAND or j not in CAND:
        die('the doubles file holds a pair outside the candidate union, %d with %d'
            % (i, j))
    if r[dh['from_construct_i']] != BG[i - 1] or r[dh['from_construct_j']] != BG[j - 1]:
        die('the doubles file disagrees with the rebuilt background at %d or %d' % (i, j))
    if (i, x, j, y) in seen:
        die('the doubles file repeats %s%d%s with %s%d%s'
            % (BG[i - 1], i, x, BG[j - 1], j, y))
    seen.add((i, x, j, y))
    if (i, x) not in sing or (j, y) not in sing:
        die('a half of %s%d%s with %s%d%s has no single score'
            % (BG[i - 1], i, x, BG[j - 1], j, y))
    p = float(r[dh['p_target']])
    pi, pj = sing[(i, x)], sing[(j, y)]
    D.append({'pos_i': i, 'to_i': x, 'from_construct_i': BG[i - 1],
              'pos_j': j, 'to_j': y, 'from_construct_j': BG[j - 1],
              'p_target': p, 'best_single_of_pair': max(pi, pj),
              'fold_vs_baseline': p / summ['checks']['baseline'],
              'fold_vs_best_single_of_pair': p / max(pi, pj),
              'crosses': p >= THRESHOLD,
              'spares_the_set': (i not in PSET) and (j not in PSET)})
if len(D) != EXP_DOUBLES:
    die('%d doubles read, %d expected over %d pairs' % (len(D), EXP_DOUBLES, EXP_PAIRS))
D.sort(key=lambda r: -r['p_target'])

key4 = [r for r in D if r['pos_i'] == R_T[0] and r['to_i'] == R_TO[0]
        and r['pos_j'] == R_T[1] and r['to_j'] == R_TO[1]]
if len(key4) != 1:
    die('the sweep holds %d rows for %s, expected one' % (len(key4), CH['label']))
P4 = key4[0]['p_target']
s7 = json.load(open(S7))
s7row = [x for x in s7['part_a']['profile'] if x['chicken_pos'] == TARGET]
if not s7row:
    die('step 7 holds no profile row for chicken %d, so the two runs cannot be tied '
        'to one scale' % TARGET)
S7P = s7row[0]['p_admissible']
rel = abs(P4 - S7P) / max(P4, S7P, 1e-30)
if rel > TOL:
    die('%s scores %.6e here and %.6e in step 7, a relative difference of %.2e, so the '
        'two runs are not on one scale' % (CH['label'], P4, S7P, rel))
log('scale check: %s at %.6e here against %.6e in step 7, relative difference %.2e'
    % (CH['label'], P4, S7P, rel))

BASE = summ['checks']['baseline']
if BASE >= THRESHOLD:
    die('the background alone reaches %.6e, which is at or above the threshold' % BASE)
cross = [r for r in D if r['crosses']]
spare = [r for r in cross if r['spares_the_set']]
better = [r for r in D if r['p_target'] > P4]
better_spare = [r for r in better if r['spares_the_set']]
folds = [r['fold_vs_best_single_of_pair'] for r in D]
best_single = max(sing.items(), key=lambda kv: kv[1])
comparable = [(k, v) for k, v in sing.items() if k in prev]
bg_fold = sorted(v / prev[k] for k, v in comparable)

log('')
log('chicken %d to human %s at human %d, in the background the construct has'
    % (TARGET, PSET[TARGET]['human_aa'], PSET[TARGET]['human_pos']))
log('  baseline        : %.6e, against %.6e for the same position in step 5 primaries'
    % (BASE, s7row[0]['p_primaries']))
log('  %s      : %.6e in this background, %.6e where step 4 chose it'
    % (CH['label'], P4, CH['p']))
log('  singles         : %d, strongest %s%d%s at %.6e, inside the set %s'
    % (len(sing), BG[best_single[0][0] - 1], best_single[0][0], best_single[0][1],
       best_single[1], best_single[0][0] in PSET))
if comparable:
    log('  against step 2  : %d singles in both, median %.2f of the wild-type value, '
        '5th and 95th percentile %.2f and %.2f'
        % (len(comparable), statistics.median(bg_fold),
           bg_fold[int(0.05 * (len(bg_fold) - 1))],
           bg_fold[int(0.95 * (len(bg_fold) - 1))]))
log('  doubles         : %d over %d pairs, %d crossing %.0e, %d of those sparing the set'
    % (len(D), EXP_PAIRS, len(cross), THRESHOLD, len(spare)))
log('  beating %s: %d, of which %d spare the set'
    % (CH['label'], len(better), len(better_spare)))
log('  fold over the better half: median %.2f, max %.2f'
    % (statistics.median(folds), max(folds)))
log('  the five strongest')
for r in D[:5]:
    log('    %s%d%s with %s%d%s at %.6e, %.1f-fold over baseline, %.2f-fold over the '
        'better half, spares the set %s'
        % (r['from_construct_i'], r['pos_i'], r['to_i'],
           r['from_construct_j'], r['pos_j'], r['to_j'], r['p_target'],
           r['fold_vs_baseline'], r['fold_vs_best_single_of_pair'], r['spares_the_set']))
if spare:
    r = spare[0]
    log('  strongest crossing pair that spares the set: %s%d%s with %s%d%s at %.6e'
        % (r['from_construct_i'], r['pos_i'], r['to_i'],
           r['from_construct_j'], r['pos_j'], r['to_j'], r['p_target']))
else:
    log('  no crossing pair leaves the pharmacophore set untouched')
log('  verdict         : %s'
    % ('chicken %d is rescued in the construct background by a pair that spares the set'
       % TARGET if spare else
       ('chicken %d is rescued only by a pair that overwrites the set' % TARGET if cross
        else 'no pair at these %d positions lifts chicken %d above %.0e in the background '
             'the construct has' % (len(CAND), TARGET, THRESHOLD))))

json.dump({'source': SRC, 'threshold': THRESHOLD, 'target': TARGET,
           'human_pos': PSET[TARGET]['human_pos'], 'human_aa': PSET[TARGET]['human_aa'],
           'background': {'primaries_in': len(PRIM) - 1,
                          'rescuers_held': {str(k): v for k, v in KEEP.items()},
                          'positions_released': R_T,
                          'differs_from_chicken': sum(1 for a, b in zip(BG, chicken)
                                                      if a != b),
                          'rebuild_matches_run': True},
           'baseline': BASE, 'baseline_step5_primaries': s7row[0]['p_primaries'],
           'step4_pair': {'label': CH['label'], 'positions': R_T, 'to': R_TO,
                          'p_wt_background': CH['p'], 'p_construct': P4,
                          'p_construct_step7': S7P, 'relative_difference': rel,
                          'fold_vs_wt_background': P4 / CH['p']},
           'coupled': COUP, 'nearest': NEAR, 'candidates': CAND,
           'candidates_in_set': IN_SET,
           'n_pairs': EXP_PAIRS, 'n_doubles': len(D), 'n_singles': len(sing),
           'background_effect': {'n_comparable': len(comparable),
                                 'fold_median': statistics.median(bg_fold) if bg_fold
                                 else None,
                                 'fold_p05': bg_fold[int(0.05 * (len(bg_fold) - 1))]
                                 if bg_fold else None,
                                 'fold_p95': bg_fold[int(0.95 * (len(bg_fold) - 1))]
                                 if bg_fold else None},
           'best_single': {'pos': best_single[0][0], 'to': best_single[0][1],
                           'p': best_single[1],
                           'in_set': best_single[0][0] in PSET},
           'n_crossing': len(cross), 'n_crossing_sparing_the_set': len(spare),
           'n_better_than_step4_pair': len(better),
           'n_better_sparing_the_set': len(better_spare),
           'best_double': D[0],
           'best_sparing_double': spare[0] if spare else None,
           'top_doubles': D[:20],
           'fold_over_better_half': {'median': statistics.median(folds),
                                     'max': max(folds)}},
          open(OUT, 'w'), indent=1)
log('')
log('wrote ' + OUT)
with open(LOG, 'w') as f:
    f.write('\n'.join(lines) + '\n')
