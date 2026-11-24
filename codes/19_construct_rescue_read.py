import glob
import json
import os
import statistics
import sys

POSSET = 'data/pharmacophore_position_set.tsv'
RISK = 'data/primary_substitution_risk.tsv'
SCAN = 'data/esm2_step2/rescue_scan.tsv'
CON = 'data/construct_report.json'
MAP = 'structures/7BC6_to_human.json'
PARKER = 'data/parker_residues.json'
HFM = 'data/hfm_positions.tsv'
OUT = 'data/construct_rescue_report.json'
LOG = 'data/19_construct_rescue_read_log.txt'
THRESHOLD = 1e-3
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
    c = [d for d in sorted(glob.glob('data/esm2_step6*')) if os.path.isdir(d)]
    if not c:
        die('no esm2_step6 folder under data/. Run '
            'codes/SLC46A1_step6_construct_rescue.ipynb on Colab and copy esm2_step6/ '
            'from Drive into data/, or pass its path.')
    if len(c) > 1:
        die('several candidates under data/, so name the one to read: %s' % c)
    return c[0]


SRC = pick()
for p in (POSSET, RISK, SCAN, CON, MAP, PARKER, HFM):
    if not os.path.exists(p):
        die(p + ' not found')
need = ['construct_rescue_scan.tsv', 'construct_rescue_profile.tsv']
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
    die('sequence lengths %d and %d, expected 473 and 459' % (len(chicken), len(human)))

rows = [l.rstrip('\n').split('\t') for l in open(POSSET)]
pc = {k: i for i, k in enumerate(rows[0])}
PSET = {}
for r in rows[1:]:
    cp = int(r[pc['chicken_pos']])
    PSET[cp] = {'chicken_aa': r[pc['chicken_aa']], 'human_aa': r[pc['human_aa']],
                'human_pos': int(r[pc['human_pos']]),
                'differs': r[pc['differs']] == '1',
                'parker': r[pc['in_parker_named']] == '1',
                'role': r[pc['role_or_note']]}
PRIM = sorted(p for p, v in PSET.items() if v['differs'])
if len(PRIM) != 20:
    die('%d primary substitutions, expected 20' % len(PRIM))

hu2ck = {int(k): v for k, v in json.load(open(MAP))['hu2ck'].items()}
PROTECTED = set()
for k, v in json.load(open(PARKER))['residues'].items():
    cp, hp = int(k), v['human']
    if chicken[cp - 1] == human[hp - 1]:
        PROTECTED.add(cp)
for r in [l.rstrip('\n').split('\t') for l in open(HFM)][1:]:
    cp = hu2ck.get(int(r[0]))
    if cp is None:
        die('hereditary folate malabsorption position %s has no chicken partner' % r[0])
    PROTECTED.add(cp)

con = json.load(open(CON))
TARGETS = sorted(int(t) for t in con['part_b']['below_threshold']['primaries'])
PROF5 = {x['chicken_pos']: x for x in con['part_b']['profile']}
rr = [l.rstrip('\n').split('\t') for l in open(RISK)]
rcl = {k: i for i, k in enumerate(rr[0])}
WT = {int(r[rcl['chicken_pos']]): float(r[rcl['p_human_wt_background']]) for r in rr[1:]}


def background(leave_out):
    s = list(chicken)
    for p in PRIM:
        if p == leave_out:
            continue
        s[p - 1] = PSET[p]['human_aa']
    return ''.join(s)


log('result             : %s' % SRC)
log('targets from %s: %s' % (CON, TARGETS))
log('primaries %d, protected positions %d' % (len(PRIM), len(PROTECTED)))

PREV = {}
f = open(SCAN)
h = f.readline().rstrip('\n').split('\t')
s2 = {k: i for i, k in enumerate(h)}
for line in f:
    r = line.rstrip('\n').split('\t')
    PREV[(int(r[s2['target_pos']]), int(r[s2['rescuer_pos']]),
          r[s2['rescuer_to']])] = float(r[s2['p_target']])
if not PREV:
    die(SCAN + ' produced no rows')

srows = [l.rstrip('\n').split('\t')
         for l in open(os.path.join(SRC, 'construct_rescue_scan.tsv'))]
sh = {k: i for i, k in enumerate(srows[0])}
for k in ('target', 'pos', 'from_construct', 'to', 'p_construct', 'p_wt_step2'):
    if k not in sh:
        die('the scan lacks the %s column; it has %s' % (k, srows[0]))
BG = {t: background(t) for t in TARGETS}
R, seen, nowt = [], set(), 0
for r in srows[1:]:
    t, j, y = int(r[sh['target']]), int(r[sh['pos']]), r[sh['to']]
    if t not in BG:
        die('the scan holds target %d, which is not in the step 5 target list' % t)
    if r[sh['from_construct']] != BG[t][j - 1]:
        die('the scan says the construct has %s at chicken %d for target %d, the '
            'background rebuilt here has %s'
            % (r[sh['from_construct']], j, t, BG[t][j - 1]))
    if y == BG[t][j - 1]:
        die('the scan holds a non-substitution at chicken %d for target %d' % (j, t))
    if (t, j, y) in seen:
        die('the scan repeats target %d with %s%d%s' % (t, BG[t][j - 1], j, y))
    seen.add((t, j, y))
    pv = PREV.get((t, j, y))
    stated = r[sh['p_wt_step2']].strip()
    if pv is None:
        nowt += 1
        if stated:
            die('the scan states a step 2 value for target %d %s%d%s, which the step 2 '
                'scan does not hold' % (t, BG[t][j - 1], j, y))
    else:
        if not stated:
            die('the scan states no step 2 value for target %d %s%d%s, which the step 2 '
                'scan does hold' % (t, BG[t][j - 1], j, y))
        if abs(float(stated) - pv) / max(float(stated), pv, 1e-30) > TOL:
            die('the scan states a step 2 value of %s for target %d %s%d%s, the step 2 '
                'scan holds %.6e' % (stated, t, BG[t][j - 1], j, y, pv))
    R.append({'target': t, 'pos': j, 'from_construct': BG[t][j - 1],
              'from_chicken': chicken[j - 1], 'to': y,
              'separation': abs(j - t),
              'p_construct': float(r[sh['p_construct']]),
              'p_wt_step2': pv,
              'crosses_construct': float(r[sh['p_construct']]) >= THRESHOLD,
              'is_primary': j in PRIM, 'protected': j in PROTECTED})
EXP = len(TARGETS) * (len(chicken) - 1) * 19
if len(R) != EXP:
    die('%d scan rows, %d expected for %d targets over %d positions'
        % (len(R), EXP, len(TARGETS), len(chicken) - 1))
exp_nowt = len(TARGETS) * (len(PRIM) - 1)
if nowt != exp_nowt:
    die('%d rows have no step 2 counterpart, %d expected, one per other primary per '
        'target' % (nowt, exp_nowt))
log('scan               : %d rows over %d targets, every from-residue agrees with the '
    'rebuilt background' % (len(R), len(TARGETS)))
log('                     %d rows have no step 2 counterpart, which is one per other '
    'primary per target, the chicken residue there' % nowt)


def admissible(policy, j, t):
    if j == t:
        return False
    if policy == 'published':
        return True
    if policy == 'annotated':
        return j not in PROTECTED and j not in PRIM
    die('unknown policy ' + policy)


POLICIES = ['published', 'annotated']
log('')
log('per target, against the construct background')
log('  target  baseline here  baseline wild type  crossing  best substitution  '
    'p here        step 2 value  fold')
PERT = {}
for t in TARGETS:
    here = [r for r in R if r['target'] == t]
    cross = [r for r in here if r['crosses_construct']]
    best = max(here, key=lambda r: r['p_construct'])
    PERT[t] = {'baseline_construct': PROF5[t]['primaries'], 'baseline_wt': WT[t],
               'n_scored': len(here), 'n_crossing': len(cross), 'best': best,
               'crossing_by_policy': {
                   p: len([r for r in cross if admissible(p, r['pos'], t)])
                   for p in POLICIES}}
    log('  %6d %14.6e %19.6e %9d  %s%d%s %18.6e %13s %6s'
        % (t, PROF5[t]['primaries'], WT[t], len(cross), best['from_construct'],
           best['pos'], best['to'], best['p_construct'],
           ('%.6e' % best['p_wt_step2']) if best['p_wt_step2'] else '-',
           ('%.2f' % (best['p_construct'] / best['p_wt_step2']))
           if best['p_wt_step2'] else '-'))
    log('           crossing by policy: %s' % PERT[t]['crossing_by_policy'])

both = [r for r in R if r['p_wt_step2']]
gained = [r for r in both
          if r['crosses_construct'] and r['p_wt_step2'] < THRESHOLD]
lostc = [r for r in both
         if not r['crosses_construct'] and r['p_wt_step2'] >= THRESHOLD]
folds = sorted(r['p_construct'] / r['p_wt_step2'] for r in both)
log('')
log('the background changes which substitutions cross: %d cross here and not in step 2, '
    '%d cross in step 2 and not here, of %d comparable rows'
    % (len(gained), len(lostc), len(both)))
log('fold of the construct value over the step 2 value: median %.2f, 5th and 95th '
    'percentile %.2f and %.2f'
    % (folds[len(folds) // 2], folds[int(0.05 * len(folds))],
       folds[int(0.95 * len(folds))]))

CHOSEN = {}
for policy in POLICIES:
    picks = []
    for t in TARGETS:
        cand = [r for r in R if r['target'] == t and r['crosses_construct']
                and admissible(policy, r['pos'], t)]
        if not cand:
            picks.append({'target': t, 'found': False})
            continue
        b = max(cand, key=lambda r: r['p_construct'])
        picks.append({'target': t, 'found': True, 'pos': b['pos'], 'to': b['to'],
                      'from_construct': b['from_construct'],
                      'from_chicken': b['from_chicken'],
                      'p_construct': b['p_construct'],
                      'separation': b['separation'], 'is_primary': b['is_primary'],
                      'protected': b['protected']})
    CHOSEN[policy] = picks
    log('')
    log('policy %s, re-derived here' % policy)
    log('  target  rescuer        separation  p in the construct  primary  protected')
    for p in picks:
        if not p['found']:
            log('  %6d  %-14s %10s %19s %8s %10s'
                % (p['target'], 'none crossing', '-', '-', '-', '-'))
            continue
        log('  %6d  %s%d%s %21d %19.6e %8s %10s'
            % (p['target'], p['from_construct'], p['pos'], p['to'], p['separation'],
               p['p_construct'], p['is_primary'], p['protected']))

summ = os.path.join(SRC, 'construct_rescue_summary.json')
if os.path.exists(summ):
    run = json.load(open(summ))
    for policy in POLICIES:
        mine = {p['target']: (p.get('pos'), p.get('to')) for p in CHOSEN[policy]}
        theirs = {p['target']: (p.get('pos'), p.get('to'))
                  for p in run.get('chosen', {}).get(policy, [])}
        if mine != theirs:
            diff = [t for t in mine if mine.get(t) != theirs.get(t)]
            die('the run and this script choose different rescuers under %s at targets '
                '%s: here %s, there %s'
                % (policy, diff, [mine[t] for t in diff],
                   [theirs.get(t) for t in diff]))
    log('')
    log('the run and this script choose the same rescuers under both policies')

prows = [l.rstrip('\n').split('\t')
         for l in open(os.path.join(SRC, 'construct_rescue_profile.tsv'))]
ph = {k: i for i, k in enumerate(prows[0])}
for k in ('chicken_pos', 'p_wt_step1', 'p_primaries_step5', 'p_step5_final'):
    if k not in ph:
        die('the profile lacks the %s column; it has %s' % (k, prows[0]))
PCOLS = ['p_' + p for p in POLICIES if 'p_' + p in ph]
if len(PCOLS) != len(POLICIES):
    die('the profile lacks a column for every policy; it has %s' % prows[0])
prof = []
for r in prows[1:]:
    t = int(r[ph['chicken_pos']])
    if t not in PRIM:
        die('the profile holds chicken %d, which is not a primary' % t)
    for col, ref in (('p_wt_step1', WT[t]),
                     ('p_primaries_step5', PROF5[t]['primaries']),
                     ('p_step5_final', PROF5[t]['primaries_rescuers_primary_wins'])):
        a = float(r[ph[col]])
        if abs(a - ref) / max(a, ref, 1e-30) > TOL:
            die('the profile states %s of %.6e at chicken %d, the earlier run holds %.6e'
                % (col, a, t, ref))
    rec = {'chicken_pos': t, 'human_aa': PSET[t]['human_aa'],
           'p_wt_step1': float(r[ph['p_wt_step1']]),
           'p_primaries_step5': float(r[ph['p_primaries_step5']]),
           'p_step5_final': float(r[ph['p_step5_final']])}
    for c in PCOLS:
        rec[c] = float(r[ph[c]])
    prof.append(rec)
if sorted(x['chicken_pos'] for x in prof) != PRIM:
    die('the profile covers %d primaries, the set has %d' % (len(prof), len(PRIM)))
prof.sort(key=lambda x: x['chicken_pos'])
log('')
log('every primary with the others in place, and the three earlier columns verified '
    'against their own runs')
log('chicken  human   wild type   primaries   step 5 final  %s'
    % '  '.join('%-13s' % c for c in PCOLS))
for x in prof:
    log('%7d %6s %11.3e %11.3e %14.3e  %s'
        % (x['chicken_pos'], x['human_aa'], x['p_wt_step1'], x['p_primaries_step5'],
           x['p_step5_final'], '  '.join('%-13.3e' % x[c] for c in PCOLS)))
log('')
BELOW = {}
for key in ['p_step5_final'] + PCOLS:
    below = [x['chicken_pos'] for x in prof if x[key] < THRESHOLD]
    BELOW[key] = below
    log('%-22s below %.0e at %d of %d primaries%s'
        % (key, THRESHOLD, len(below), len(prof), (': %s' % below) if below else ''))
for c in PCOLS:
    med = statistics.median([x[c] / x['p_step5_final'] for x in prof])
    log('%-22s median fold over the step 5 construct %.2f' % (c, med))
best_pol = min(PCOLS, key=lambda c: len(BELOW[c]))
log('')
log('verdict            : searching in the construct background leaves %d of %d '
    'primaries below %.0e under %s, against %d under the step 5 rescuers'
    % (len(BELOW[best_pol]), len(prof), THRESHOLD, best_pol[2:],
       len(BELOW['p_step5_final'])))

json.dump({'source': SRC, 'threshold': THRESHOLD, 'targets': TARGETS,
           'policies': POLICIES, 'n_scan_rows': len(R),
           'n_without_step2_counterpart': nowt,
           'background_effect': {'n_comparable': len(both),
                                 'crossing_here_not_step2': len(gained),
                                 'crossing_step2_not_here': len(lostc),
                                 'fold_median': folds[len(folds) // 2],
                                 'fold_p05': folds[int(0.05 * len(folds))],
                                 'fold_p95': folds[int(0.95 * len(folds))]},
           'per_target': {str(t): {k: v for k, v in PERT[t].items()} for t in TARGETS},
           'chosen': CHOSEN, 'profile': prof, 'below_threshold': BELOW,
           'median_fold_over_step5': {
               c: statistics.median([x[c] / x['p_step5_final'] for x in prof])
               for c in PCOLS},
           'best_policy': best_pol[2:]},
          open(OUT, 'w'), indent=1)
log('')
log('wrote ' + OUT)
with open(LOG, 'w') as f:
    f.write('\n'.join(lines) + '\n')
