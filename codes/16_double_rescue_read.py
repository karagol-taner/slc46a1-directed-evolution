import glob
import json
import os
import statistics
import sys

DBL = 'double_rescue.tsv'
SNG = 'singles_recheck.tsv'
CGR = 'data/coupling_guided_rescue.json'
RISK = 'data/primary_substitution_risk.tsv'
SCAN = 'data/esm2_step2/rescue_scan.tsv'
OUT = 'data/double_rescue_report.json'
LOG = 'data/16_double_rescue_read_log.txt'
THRESHOLD = 1e-3
TARGET_CK = 407
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
    cand = [d for d in sorted(glob.glob('data/esm2_step4*')) if os.path.isdir(d)]
    if not cand:
        die('no esm2_step4 folder under data/. Run '
            'codes/SLC46A1_step4_double_rescue.ipynb on Colab and copy esm2_step4/ '
            'from Drive into data/, or pass its path.')
    if len(cand) > 1:
        die('several candidates under data/, so name the one to read: %s' % cand)
    return cand[0]


SRC = pick()
for p in (CGR, RISK, SCAN):
    if not os.path.exists(p):
        die(p + ' not found')
dpath, spath = os.path.join(SRC, DBL), os.path.join(SRC, SNG)
for p in (dpath, spath):
    if not os.path.exists(p):
        die(p + ' not found, so the run has not finished its last cell')

cgr = json.load(open(CGR))
T = cgr['targets'][str(TARGET_CK)]
CAND = [int(r['chicken_pos']) for r in T['top_coupled']]
N_POS = len(CAND)
EXP_PAIRS = N_POS * (N_POS - 1) // 2
EXP_DOUBLES = EXP_PAIRS * 19 * 19
EXP_SINGLES = N_POS * 19
if abs(cgr['threshold'] - THRESHOLD) > 1e-12:
    die('the sweep used a threshold of %g, this script uses %g'
        % (cgr['threshold'], THRESHOLD))

ent = {}
for rec in open('seq/7BC6.fasta').read().split('>'):
    if rec.strip():
        parts = rec.strip().split('\n')
        ent[parts[0].split('|')[0]] = ''.join(parts[1:])
chicken = ent['7BC6_1'][:473]
if len(chicken) != 473:
    die('chicken native length %d, expected 473' % len(chicken))

rows = [l.rstrip('\n').split('\t') for l in open(RISK)]
col = {k: i for i, k in enumerate(rows[0])}
tr = [r for r in rows[1:] if int(r[col['chicken_pos']]) == TARGET_CK]
if len(tr) != 1:
    die('chicken %d appears %d times in %s' % (TARGET_CK, len(tr), RISK))
P_BASELINE = float(tr[0][col['p_human_wt_background']])
if P_BASELINE >= THRESHOLD:
    die('baseline %.3e is not below the threshold' % P_BASELINE)
TARGET_HU_POS = int(tr[0][col['human_pos']])
TARGET_HU_AA = tr[0][col['human_aa']]

prev = {}
f = open(SCAN)
head = f.readline().rstrip('\n').split('\t')
sc = {k: i for i, k in enumerate(head)}
for line in f:
    r = line.rstrip('\n').split('\t')
    if int(r[sc['target_pos']]) != TARGET_CK:
        continue
    prev[(int(r[sc['rescuer_pos']]), r[sc['rescuer_to']])] = float(r[sc['p_target']])
if not prev:
    die('%s holds no rows for target %d' % (SCAN, TARGET_CK))
best_prev = max(((k, v) for k, v in prev.items() if k[0] != TARGET_CK),
                key=lambda kv: kv[1])

log('result             : %s' % SRC)
log('target             : chicken %d to human %s%d, baseline %.3e from %s'
    % (TARGET_CK, TARGET_HU_AA, TARGET_HU_POS, P_BASELINE, RISK))
log('candidate positions: %d, %s' % (N_POS, CAND))
log('expected rows      : %d singles and %d doubles over %d pairs'
    % (EXP_SINGLES, EXP_DOUBLES, EXP_PAIRS))
log('step 2 for target  : %d single substitutions, strongest rescuer %s%d at %.6e, '
    'the target position itself set aside'
    % (len(prev), best_prev[0][1], best_prev[0][0], best_prev[1]))

srows = [l.rstrip('\n').split('\t') for l in open(spath)]
sh = {k: i for i, k in enumerate(srows[0])}
for k in ('pos', 'from', 'to', 'p_target'):
    if k not in sh:
        die('%s lacks the %s column; it has %s' % (spath, k, srows[0]))
sing = {}
for r in srows[1:]:
    j, y = int(r[sh['pos']]), r[sh['to']]
    if chicken[j - 1] != r[sh['from']]:
        die('%s says chicken %d is %s, the sequence says %s'
            % (spath, j, r[sh['from']], chicken[j - 1]))
    sing[(j, y)] = float(r[sh['p_target']])
if len(sing) != EXP_SINGLES:
    die('%s holds %d singles, %d expected for %d candidate positions'
        % (spath, len(sing), EXP_SINGLES, N_POS))
absent = [k for k in sing if k not in prev]
if absent:
    die('%d singles are absent from %s, the first being %s%d'
        % (len(absent), SCAN, absent[0][1], absent[0][0]))
worst = (0.0, None)
for k, v in sing.items():
    rel = abs(v - prev[k]) / max(v, prev[k], 1e-30)
    if rel > worst[0]:
        worst = (rel, (k, v, prev[k]))
if worst[0] > TOL:
    k, a, b = worst[1]
    die('the rescored singles disagree with step 2: %s%d gives %.6e in the new run and '
        '%.6e in step 2, a relative difference of %.2e above the tolerance of %.0e'
        % (k[1], k[0], a, b, worst[0], TOL))
best_single = max(sing.items(), key=lambda kv: kv[1])
log('singles            : %d read, worst relative difference against step 2 %.2e, '
    'strongest %s%d at %.6e'
    % (len(sing), worst[0], best_single[0][1], best_single[0][0], best_single[1]))

drows = [l.rstrip('\n').split('\t') for l in open(dpath)]
dh = {k: i for i, k in enumerate(drows[0])}
need = ('pos_i', 'from_i', 'to_i', 'pos_j', 'from_j', 'to_j', 'p_target')
for k in need:
    if k not in dh:
        die('%s lacks the %s column; it has %s' % (dpath, k, drows[0]))
D = []
seen = set()
for r in drows[1:]:
    i, x = int(r[dh['pos_i']]), r[dh['to_i']]
    j, y = int(r[dh['pos_j']]), r[dh['to_j']]
    if chicken[i - 1] != r[dh['from_i']] or chicken[j - 1] != r[dh['from_j']]:
        die('%s disagrees with the sequence at chicken %d or %d' % (dpath, i, j))
    if i not in CAND or j not in CAND:
        die('%s holds a pair outside the candidate list, %d with %d' % (dpath, i, j))
    key = (i, x, j, y)
    if key in seen:
        die('%s repeats the variant %s%d%s with %s%d%s'
            % (dpath, chicken[i - 1], i, x, chicken[j - 1], j, y))
    seen.add(key)
    if (i, x) not in sing or (j, y) not in sing:
        die('a half of %s%d%s with %s%d%s has no single score'
            % (chicken[i - 1], i, x, chicken[j - 1], j, y))
    p = float(r[dh['p_target']])
    pi, pj = sing[(i, x)], sing[(j, y)]
    D.append({'pos_i': i, 'from_i': chicken[i - 1], 'to_i': x,
              'pos_j': j, 'from_j': chicken[j - 1], 'to_j': y,
              'p_target': p, 'p_single_i': pi, 'p_single_j': pj,
              'best_single_of_pair': max(pi, pj),
              'fold_vs_baseline': p / P_BASELINE,
              'fold_vs_best_single_of_pair': p / max(pi, pj),
              'crosses': p >= THRESHOLD})
if len(D) != EXP_DOUBLES:
    die('%s holds %d doubles, %d expected for %d pairs at 19 by 19'
        % (dpath, len(D), EXP_DOUBLES, EXP_PAIRS))
npairs = len({(r['pos_i'], r['pos_j']) for r in D})
if npairs != EXP_PAIRS:
    die('%s covers %d pairs, %d expected' % (dpath, npairs, EXP_PAIRS))
D.sort(key=lambda r: -r['p_target'])
log('doubles            : %d read over %d pairs, every wild type and every half checked'
    % (len(D), npairs))

crossing = [r for r in D if r['crosses']]
above = [r for r in D if r['p_target'] > best_single[1]]
folds = [r['fold_vs_best_single_of_pair'] for r in D]
log('')
log('crossing %.0e      : %d of %d doubles' % (THRESHOLD, len(crossing), len(D)))
log('above the strongest single: %d of %d' % (len(above), len(D)))
log('fold over the better half of the pair: median %.2f, 95th percentile %.2f, max %.2f'
    % (statistics.median(folds),
       sorted(folds)[int(0.95 * (len(folds) - 1))], max(folds)))
log('')
log('the ten strongest doubles')
log('  variant                      likelihood  fold over baseline  fold over the '
    'better half  crosses')
for r in D[:10]:
    log('  %s%d%s with %s%d%s %18.6e %19.1f %26.2f  %s'
        % (r['from_i'], r['pos_i'], r['to_i'], r['from_j'], r['pos_j'], r['to_j'],
           r['p_target'], r['fold_vs_baseline'], r['fold_vs_best_single_of_pair'],
           'yes' if r['crosses'] else 'no'))

by_pair = {}
for r in D:
    by_pair.setdefault((r['pos_i'], r['pos_j']), []).append(r)
pair_best = sorted(((max(v, key=lambda r: r['p_target'])) for v in by_pair.values()),
                   key=lambda r: -r['p_target'])
log('')
log('the strongest double at each of the five best pairs')
for r in pair_best[:5]:
    log('  chicken %d with %d: %s%d%s with %s%d%s at %.6e, %.2f-fold over the better half'
        % (r['pos_i'], r['pos_j'], r['from_i'], r['pos_i'], r['to_i'],
           r['from_j'], r['pos_j'], r['to_j'], r['p_target'],
           r['fold_vs_best_single_of_pair']))

log('')
if crossing:
    log('verdict            : %d double substitution%s lift%s chicken %d above %.0e, so '
        'the position is rescued by a pair where no single substitution reaches it'
        % (len(crossing), '' if len(crossing) == 1 else 's',
           's' if len(crossing) == 1 else '', TARGET_CK, THRESHOLD))
else:
    log('verdict            : no double substitution at the %d coupled positions lifts '
        'chicken %d above %.0e. The strongest reaches %.6e against the threshold, and '
        'the measurement bounds the pairs at these positions and not the whole double '
        'space, which holds about %.1f million variants.'
        % (N_POS, TARGET_CK, THRESHOLD, D[0]['p_target'],
           (473 * 19) * (473 * 19 - 19) / 2 / 1e6))

json.dump({'source': SRC, 'threshold': THRESHOLD, 'target_chicken_pos': TARGET_CK,
           'target_human_pos': TARGET_HU_POS, 'target_human_aa': TARGET_HU_AA,
           'p_baseline': P_BASELINE, 'candidate_positions': CAND,
           'n_pairs': npairs, 'n_doubles': len(D), 'n_singles': len(sing),
           'worst_relative_difference_vs_step2': worst[0],
           'best_single': {'pos': best_single[0][0], 'to': best_single[0][1],
                           'p': best_single[1]},
           'best_double': D[0],
           'n_doubles_crossing': len(crossing),
           'n_doubles_above_best_single': len(above),
           'fold_over_better_half': {
               'median': statistics.median(folds),
               'p95': sorted(folds)[int(0.95 * (len(folds) - 1))],
               'max': max(folds)},
           'top_doubles': D[:20],
           'best_per_pair': pair_best},
          open(OUT, 'w'), indent=1)
log('')
log('wrote ' + OUT)
with open(LOG, 'w') as f:
    f.write('\n'.join(lines) + '\n')
