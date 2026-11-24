import glob
import json
import os
import statistics
import sys

POSSET = 'data/pharmacophore_position_set.tsv'
RISK = 'data/primary_substitution_risk.tsv'
SCAN = 'data/esm2_step2/rescue_scan.tsv'
STRUCT = 'data/structural_candidates.json'
SWEEP = 'data/rescue_sweep_report.json'
D4 = 'data/double_rescue_report.json'
OUT = 'data/construct_report.json'
LOG = 'data/18_construct_read_log.txt'
THRESHOLD = 1e-3
TARGET_A = 144
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
    cand = [d for d in sorted(glob.glob('data/esm2_step5*')) if os.path.isdir(d)]
    if not cand:
        die('no esm2_step5 folder under data/. Run '
            'codes/SLC46A1_step5_construct.ipynb on Colab and copy esm2_step5/ from '
            'Drive into data/, or pass its path.')
    if len(cand) > 1:
        die('several candidates under data/, so name the one to read: %s' % cand)
    return cand[0]


SRC = pick()
for p in (POSSET, RISK, SCAN, STRUCT, SWEEP, D4):
    if not os.path.exists(p):
        die(p + ' not found')
need = {'a144_double_rescue.tsv', 'a144_singles_recheck.tsv', 'construct_profile.tsv'}
missing = [f for f in sorted(need) if not os.path.exists(os.path.join(SRC, f))]
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
    die('sequence lengths are %d and %d, expected 473 and 459' % (len(chicken), len(human)))

rows = [l.rstrip('\n').split('\t') for l in open(POSSET)]
pcol = {k: i for i, k in enumerate(rows[0])}
PSET = {}
for r in rows[1:]:
    cp = int(r[pcol['chicken_pos']])
    PSET[cp] = {'chicken_aa': r[pcol['chicken_aa']], 'human_aa': r[pcol['human_aa']],
                'human_pos': int(r[pcol['human_pos']]),
                'differs': r[pcol['differs']] == '1',
                'role': r[pcol['role_or_note']]}
PRIM = sorted(p for p, v in PSET.items() if v['differs'])
log('result             : %s' % SRC)
log('position set       : %d positions, %d primary humanising substitutions'
    % (len(PSET), len(PRIM)))

risk = {}
rr = [l.rstrip('\n').split('\t') for l in open(RISK)]
rc = {k: i for i, k in enumerate(rr[0])}
for r in rr[1:]:
    risk[int(r[rc['chicken_pos']])] = float(r[rc['p_human_wt_background']])

sc = json.load(open(STRUCT))
if str(TARGET_A) not in sc['targets']:
    die('%s has no ranking for chicken %d' % (STRUCT, TARGET_A))
CAND = [int(c['chicken_pos']) for c in sc['targets'][str(TARGET_A)]['candidates']]
DIST = {int(c['chicken_pos']): c for c in sc['targets'][str(TARGET_A)]['candidates']}
EXP_PAIRS = len(CAND) * (len(CAND) - 1) // 2
EXP_DOUBLES = EXP_PAIRS * 19 * 19
EXP_SINGLES = len(CAND) * 19

prev = {}
f = open(SCAN)
h = f.readline().rstrip('\n').split('\t')
s2 = {k: i for i, k in enumerate(h)}
for line in f:
    r = line.rstrip('\n').split('\t')
    if int(r[s2['target_pos']]) == TARGET_A:
        prev[(int(r[s2['rescuer_pos']]), r[s2['rescuer_to']])] = float(r[s2['p_target']])
if not prev:
    die('%s holds no rows for target %d' % (SCAN, TARGET_A))

srows = [l.rstrip('\n').split('\t') for l in open(os.path.join(SRC, 'a144_singles_recheck.tsv'))]
sh = {k: i for i, k in enumerate(srows[0])}
sing = {}
for r in srows[1:]:
    j, y = int(r[sh['pos']]), r[sh['to']]
    if chicken[j - 1] != r[sh['from']]:
        die('the singles file says chicken %d is %s, the sequence says %s'
            % (j, r[sh['from']], chicken[j - 1]))
    sing[(j, y)] = float(r[sh['p_target']])
if len(sing) != EXP_SINGLES:
    die('%d singles read, %d expected for %d candidate positions'
        % (len(sing), EXP_SINGLES, len(CAND)))
worst = (0.0, None)
for k, v in sing.items():
    if k not in prev:
        die('%s%d is absent from the step 2 scan' % (k[1], k[0]))
    rel = abs(v - prev[k]) / max(v, prev[k], 1e-30)
    if rel > worst[0]:
        worst = (rel, (k, v, prev[k]))
if worst[0] > TOL:
    k, a, b = worst[1]
    die('the rescored singles disagree with step 2 at %s%d, %.6e here and %.6e there, '
        '%.2e above the tolerance of %.0e' % (k[1], k[0], a, b, worst[0], TOL))
best_single = max(sing.items(), key=lambda kv: kv[1])

drows = [l.rstrip('\n').split('\t') for l in open(os.path.join(SRC, 'a144_double_rescue.tsv'))]
dh = {k: i for i, k in enumerate(drows[0])}
D, seen = [], set()
for r in drows[1:]:
    i, x = int(r[dh['pos_i']]), r[dh['to_i']]
    j, y = int(r[dh['pos_j']]), r[dh['to_j']]
    if chicken[i - 1] != r[dh['from_i']] or chicken[j - 1] != r[dh['from_j']]:
        die('the doubles file disagrees with the sequence at chicken %d or %d' % (i, j))
    if i not in CAND or j not in CAND:
        die('the doubles file holds a pair outside the candidate list, %d with %d' % (i, j))
    if (i, x, j, y) in seen:
        die('the doubles file repeats %s%d%s with %s%d%s'
            % (chicken[i - 1], i, x, chicken[j - 1], j, y))
    seen.add((i, x, j, y))
    p = float(r[dh['p_target']])
    pi, pj = sing[(i, x)], sing[(j, y)]
    D.append({'pos_i': i, 'from_i': chicken[i - 1], 'to_i': x,
              'pos_j': j, 'from_j': chicken[j - 1], 'to_j': y, 'p_target': p,
              'best_single_of_pair': max(pi, pj),
              'fold_vs_baseline': p / risk[TARGET_A],
              'fold_vs_best_single_of_pair': p / max(pi, pj),
              'crosses': p >= THRESHOLD})
if len(D) != EXP_DOUBLES:
    die('%d doubles read, %d expected over %d pairs' % (len(D), EXP_DOUBLES, EXP_PAIRS))
D.sort(key=lambda r: -r['p_target'])
crossing = [r for r in D if r['crosses']]
folds = [r['fold_vs_best_single_of_pair'] for r in D]
log('')
log('part A, chicken %d by pairs at the %d nearest positions in the coordinates'
    % (TARGET_A, len(CAND)))
log('  candidates       : %s' % CAND)
log('  baseline         : %.3e from %s' % (risk[TARGET_A], RISK))
log('  singles          : %d, worst relative difference against step 2 %.2e, strongest '
    '%s%d at %.6e' % (len(sing), worst[0], best_single[0][1], best_single[0][0],
                      best_single[1]))
log('  doubles          : %d over %d pairs, %d crossing %.0e'
    % (len(D), EXP_PAIRS, len(crossing), THRESHOLD))
log('  fold over the better half: median %.2f, max %.2f'
    % (statistics.median(folds), max(folds)))
log('  the five strongest')
for r in D[:5]:
    q = [x for x in (r['pos_i'], r['pos_j'])]
    log('    %s%d%s with %s%d%s at %.6e, %.1f-fold over baseline, %.2f-fold over the '
        'better half, distances %s'
        % (r['from_i'], r['pos_i'], r['to_i'], r['from_j'], r['pos_j'], r['to_j'],
           r['p_target'], r['fold_vs_baseline'], r['fold_vs_best_single_of_pair'],
           ' and '.join('%.2f' % DIST[p]['min_dist'] for p in q)))
log('  verdict          : %s'
    % ('chicken %d is rescued by a pair, %d of %d crossing'
       % (TARGET_A, len(crossing), len(D)) if crossing else
       'no pair at these positions lifts chicken %d above %.0e, and the measurement '
       'bounds those pairs and not the whole double space' % (TARGET_A, THRESHOLD)))

prows = [l.rstrip('\n').split('\t') for l in open(os.path.join(SRC, 'construct_profile.tsv'))]
ph = {k: i for i, k in enumerate(prows[0])}
for k in ('chicken_pos', 'human_aa', 'p_baseline_step1', 'wt', 'primaries'):
    if k not in ph:
        die('the profile lacks the %s column; it has %s' % (k, prows[0]))
BACKS = [c for c in prows[0] if c in
         ('wt', 'primaries', 'primaries_rescuers_primary_wins',
          'primaries_rescuers_rescue_wins')]
prof = []
for r in prows[1:]:
    cp = int(r[ph['chicken_pos']])
    if cp not in PSET or not PSET[cp]['differs']:
        die('the profile holds chicken %d, which is not a primary' % cp)
    if r[ph['human_aa']] != PSET[cp]['human_aa']:
        die('the profile says human %s at chicken %d, the position set says %s'
            % (r[ph['human_aa']], cp, PSET[cp]['human_aa']))
    rec = {'chicken_pos': cp, 'human_aa': r[ph['human_aa']],
           'role': PSET[cp]['role'], 'human_pos': PSET[cp]['human_pos']}
    for b in BACKS:
        rec[b] = float(r[ph[b]])
    prof.append(rec)
if sorted(x['chicken_pos'] for x in prof) != PRIM:
    die('the profile covers %d primaries, the position set has %d'
        % (len(prof), len(PRIM)))
worstw = (0.0, None)
for x in prof:
    a, b = x['wt'], risk[x['chicken_pos']]
    rel = abs(a - b) / max(a, b, 1e-30)
    if rel > worstw[0]:
        worstw = (rel, (x['chicken_pos'], a, b))
if worstw[0] > TOL:
    cp, a, b = worstw[1]
    die('the wild-type column does not reproduce step 1 at chicken %d, %.6e against '
        '%.6e, %.2e above the tolerance of %.0e' % (cp, a, b, worstw[0], TOL))
prof.sort(key=lambda x: x['chicken_pos'])
log('')
log('part B, every primary scored with the others already written in')
log('  wild-type column reproduces step 1 at all %d primaries, worst relative '
    'difference %.2e' % (len(prof), worstw[0]))
log('  chicken  human  %s' % '  '.join('%-13s' % b[:13] for b in BACKS))
for x in prof:
    log('  %7d %6s  %s'
        % (x['chicken_pos'], x['human_aa'],
           '  '.join('%-13.3e' % x[b] for b in BACKS)))
log('')
for b in BACKS:
    below = [x['chicken_pos'] for x in prof if x[b] < THRESHOLD]
    log('  %-34s below %.0e at %d of %d primaries%s'
        % (b, THRESHOLD, len(below), len(prof),
           (': %s' % below) if below else ''))
base = {x['chicken_pos']: x['wt'] for x in prof}
for b in BACKS:
    if b == 'wt':
        continue
    fold = [x[b] / base[x['chicken_pos']] for x in prof]
    worse = [x['chicken_pos'] for x in prof if x[b] < base[x['chicken_pos']]]
    log('  %-34s median fold over wild type %.2f, lower than wild type at %d primaries'
        % (b, statistics.median(fold), len(worse)))

sweep = json.load(open(SWEEP))
d4 = json.load(open(D4))
RADIUS = sweep.get('published_radius', 3)
alone = {}
for t, b in sweep['per_target'].items():
    ok = [r for r in b.get('published', []) if r['crosses'] and r['separation'] <= RADIUS]
    if ok:
        x = max(ok, key=lambda r: r['p_target'])
        alone[int(t)] = {'p': x['p_target'],
                         'variant': '%s%s' % (x['rescuer_pos'], x['rescuer_to']),
                         'kind': 'single within the window'}
bd4 = d4['best_double']
alone[int(d4['target_chicken_pos'])] = {
    'p': bd4['p_target'],
    'variant': '%s%d%s with %s%d%s' % (bd4['from_i'], bd4['pos_i'], bd4['to_i'],
                                       bd4['from_j'], bd4['pos_j'], bd4['to_j']),
    'kind': 'double, coupled positions'}
alone[TARGET_A] = {
    'p': D[0]['p_target'],
    'variant': '%s%d%s with %s%d%s' % (D[0]['from_i'], D[0]['pos_i'], D[0]['to_i'],
                                       D[0]['from_j'], D[0]['pos_j'], D[0]['to_j']),
    'kind': 'double, nearest positions'}
missing = [t for t in alone if t not in {x['chicken_pos'] for x in prof}]
if missing:
    die('a rescued target is absent from the construct profile: %s' % missing)
FINAL, ALT = 'primaries_rescuers_primary_wins', 'primaries_rescuers_rescue_wins'
for b in (FINAL, ALT):
    if b not in BACKS:
        die('the profile lacks the %s background; it has %s' % (b, BACKS))
pmap = {x['chicken_pos']: x for x in prof}
survival = []
for t in sorted(alone):
    x = pmap[t]
    rec = {'target': t, 'alone': alone[t]['p'], 'alone_variant': alone[t]['variant'],
           'alone_kind': alone[t]['kind'],
           FINAL: x[FINAL], ALT: x[ALT],
           'fold_final_over_alone': x[FINAL] / alone[t]['p'],
           'fold_alt_over_alone': x[ALT] / alone[t]['p'],
           'overstatement_final': alone[t]['p'] / x[FINAL],
           'holds_final': x[FINAL] >= THRESHOLD,
           'holds_alt': x[ALT] >= THRESHOLD}
    rec['verdict'] = ('both resolutions' if rec['holds_final'] and rec['holds_alt']
                      else 'only where the rescue wins' if rec['holds_alt']
                      else 'only where the engraftment wins' if rec['holds_final']
                      else 'neither')
    survival.append(rec)
holds = [r['target'] for r in survival if r['verdict'] == 'both resolutions']
lost = [r['target'] for r in survival if r['verdict'] == 'neither']
cond = [r['target'] for r in survival
        if r['verdict'] not in ('both resolutions', 'neither')]
log('')
log('does each rescue hold once every substitution is present')
log('  target  measured alone  engraftment wins  rescue wins  fold, combined over alone  '
    'holds')
for r in survival:
    log('  %6d %15.6e %17.6e %12.6e %26.3f  %s'
        % (r['target'], r['alone'], r[FINAL], r[ALT], r['fold_final_over_alone'],
           r['verdict']))
lose = [r for r in survival if r['verdict'] != 'both resolutions']
log('  holds in both resolutions: %s' % (holds or 'none'))
log('  holds in neither         : %s' % (lost or 'none'))
log('  conditional on the conflict: %s' % (cond or 'none'))
if lose:
    log('  where a rescue does not hold in both, the value measured alone overstates the '
        'combined one by a factor of %.1f to %.1f'
        % (min(r['overstatement_final'] for r in lose),
           max(r['overstatement_final'] for r in lose)))

json.dump({'source': SRC, 'threshold': THRESHOLD,
           'part_a': {'target': TARGET_A, 'candidates': CAND,
                      'n_doubles': len(D), 'n_pairs': EXP_PAIRS,
                      'n_singles': len(sing),
                      'worst_relative_difference_vs_step2': worst[0],
                      'baseline': risk[TARGET_A],
                      'best_single': {'pos': best_single[0][0],
                                      'to': best_single[0][1], 'p': best_single[1]},
                      'n_crossing': len(crossing), 'best_double': D[0],
                      'top_doubles': D[:20],
                      'fold_over_better_half': {'median': statistics.median(folds),
                                                'max': max(folds)}},
           'part_b': {'backgrounds': BACKS, 'n_primaries': len(prof),
                      'wt_vs_step1_worst_relative_difference': worstw[0],
                      'below_threshold': {b: [x['chicken_pos'] for x in prof
                                              if x[b] < THRESHOLD] for b in BACKS},
                      'median_fold_over_wt': {
                          b: statistics.median([x[b] / base[x['chicken_pos']]
                                                for x in prof])
                          for b in BACKS if b != 'wt'},
                      'n_lower_than_wt': {
                          b: sum(1 for x in prof if x[b] < base[x['chicken_pos']])
                          for b in BACKS if b != 'wt'},
                      'profile': prof},
           'survival': {'targets': survival, 'holds_both': holds, 'holds_neither': lost,
                        'conditional': cond,
                        'overstatement_range': (
                            [min(r['overstatement_final'] for r in lose),
                             max(r['overstatement_final'] for r in lose)]
                            if lose else None),
                        'lost_fold_range': (
                            [min(r['fold_final_over_alone'] for r in survival
                                 if r['verdict'] == 'neither'),
                             max(r['fold_final_over_alone'] for r in survival
                                 if r['verdict'] == 'neither')] if lost else None)}},
          open(OUT, 'w'), indent=1)
log('')
log('wrote ' + OUT)
with open(LOG, 'w') as f:
    f.write('\n'.join(lines) + '\n')
