import csv
import glob
import io
import json
import os
import sys
import zipfile

import numpy as np

SCAN = 'data/esm2_step2/rescue_scan.tsv'
MAP = 'structures/7BC6_to_human.json'
THRESHOLD = 1e-3
KS = [3, 6, 10, 20, 40, 80, 150, None]
POLICIES = ['published', 'strict', 'annotated']
WINDOW = 3
CONTACT_A = 5.0
OUT = 'data/coupling_guided_rescue.json'
LOG = 'data/15_coupling_guided_rescue_log.txt'

lines = []


def log(m):
    print(m, flush=True)
    lines.append(str(m))


def die(m):
    with open(LOG, 'w') as f:
        f.write('\n'.join(lines + ['FATAL ' + m]) + '\n')
    sys.exit('FATAL ' + m)


def pick_source():
    if len(sys.argv) > 1:
        return sys.argv[1]
    dirs = [d for d in sorted(glob.glob('data/evcouplings/*')) if os.path.isdir(d)]
    zips = sorted(glob.glob('data/evcouplings/*.zip'))
    cand = dirs + zips
    if not cand:
        die('nothing in data/evcouplings/. Put the EVcouplings download there, or '
            'pass its path.')
    if len(cand) > 1:
        die('several candidates in data/evcouplings/, so name the one to read: %s'
            % [os.path.basename(c) for c in cand])
    return cand[0]


for p in (SCAN, MAP):
    if not os.path.exists(p):
        die(p + ' not found')
SRC = pick_source()

def members(src):
    if os.path.isdir(src):
        out = {}
        for root, _, files in os.walk(src):
            for f in files:
                fp = os.path.join(root, f)
                out[os.path.relpath(fp, src)] = ('file', fp)
        return out
    with zipfile.ZipFile(src) as z:
        return {n: ('zip', n) for n in z.namelist()}


def read_text(src, entry):
    kind, ref = entry
    if kind == 'file':
        return open(ref, encoding='utf-8', errors='replace').read()
    with zipfile.ZipFile(src) as z:
        return z.read(ref).decode('utf-8', errors='replace')


M = members(SRC)
hits = [n for n in M if n.endswith('_CouplingScoresCompared_all.csv') and '__MACOSX' not in n]
if len(hits) != 1:
    die('expected one comparison table in %s, found %d' % (SRC, len(hits)))
reader = csv.DictReader(io.StringIO(read_text(SRC, M[hits[0]])))
if not {'i', 'j', 'cn', 'dist'} <= set(reader.fieldnames or []):
    die('%s lacks i, j, cn and dist; it has %s' % (hits[0], reader.fieldnames))
CS, DS = {}, {}
for r in reader:
    try:
        i, j, v = int(float(r['i'])), int(float(r['j'])), float(r['cn'])
    except (TypeError, ValueError):
        continue
    k = (min(i, j), max(i, j))
    CS[k] = v
    try:
        DS[k] = float(r['dist'])
    except (TypeError, ValueError):
        DS[k] = float('nan')
if not CS:
    die(hits[0] + ' produced no usable rows')
cv = np.array(list(CS.values()))
cs_sep = np.array([abs(a - b) for a, b in CS])
nl_sorted = np.sort(cv[cs_sep > 8])
log('couplings          : %d pairs from %s' % (len(CS), os.path.basename(SRC)))


def pct_nl(v):
    return 100.0 * float(np.searchsorted(nl_sorted, v)) / len(nl_sorted)


m = json.load(open(MAP))
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
    cp, hp = int(k), v['human']
    if chicken[cp - 1] == human[hp - 1]:
        protected.add(cp)
unmapped = []
for r in [l.rstrip('\n').split('\t') for l in open('data/hfm_positions.tsv')][1:]:
    hp = int(r[0])
    cp = hu2ck.get(hp)
    if cp is None:
        unmapped.append(hp)
    else:
        protected.add(cp)
if unmapped:
    die('hereditary folate malabsorption positions without a chicken partner: %s'
        % unmapped)
log('protected positions for the annotated policy: %d' % len(protected))

rows = [l.rstrip('\n').split('\t') for l in open(SCAN)]
head, body = rows[0], rows[1:]
c = {k: i for i, k in enumerate(head)}
if not body:
    die(SCAN + ' has no data rows')
targets = sorted({int(r[c['target_pos']]) for r in body})
baseline = {}
for r in body:
    baseline.setdefault(int(r[c['target_pos']]), float(r[c['p_baseline']]))
for tp, b in baseline.items():
    if b >= THRESHOLD:
        die('target %d baseline %.3e is not below the threshold' % (tp, b))
n_pos = len(set(int(r[c['rescuer_pos']]) for r in body))
if len(body) != n_pos * 19 * len(targets):
    die('scan has %d rows, expected %d for %d positions and %d targets'
        % (len(body), n_pos * 19 * len(targets), n_pos, len(targets)))
log('scan               : %d rows, %d rescuer positions, %d targets'
    % (len(body), n_pos, len(targets)))

variants = {}
for r in body:
    j, tp = int(r[c['rescuer_pos']]), int(r[c['target_pos']])
    variants.setdefault(tp, []).append(
        {'pos': j, 'to': r[c['rescuer_to']], 'from': r[c['rescuer_from']],
         'p': float(r[c['p_target']]), 'sep': int(r[c['separation']]),
         'fold': float(r[c['fold_change']])})


def _best_at(V, j):
    here = [v for v in V if v['pos'] == j]
    if not here:
        return {'best_substitution': None, 'best_p': None}
    b = max(here, key=lambda v: v['p'])
    return {'best_substitution': '%s%d%s' % (b['from'], b['pos'], b['to']),
            'best_p': b['p']}


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


def rank_by_coupling(target_ck):
    ht = ck2hu.get(target_ck)
    if ht is None:
        die('chicken %d has no human partner in %s' % (target_ck, MAP))
    scored, unrankable = [], 0
    for j in range(1, len(chicken) + 1):
        if j == target_ck:
            continue
        hj = ck2hu.get(j)
        k = None if hj is None else (min(ht, hj), max(ht, hj))
        if k is None or k not in CS:
            unrankable += 1
            continue
        scored.append({'pos': j, 'human': hj, 'score': CS[k],
                       'dist': DS[k] if DS[k] == DS[k] else None,
                       'sep': abs(j - target_ck)})
    scored.sort(key=lambda r: -r['score'])
    return ht, scored, unrankable


out = {}
for tp in targets:
    ht, ranked, unrankable = rank_by_coupling(tp)
    if not ranked:
        log('')
        log('chicken %d to human %d: no position can be ranked, since the run does '
            'not fit that column' % (tp, ht))
        out[tp] = {'human_pos': ht, 'rankable': 0, 'unrankable': unrankable,
                   'baseline': baseline[tp], 'window_positions': None,
                   'sweep': [], 'window': {}, 'best_overall': None}
        continue
    log('')
    log('chicken %d %s to human %s at human %d, baseline %.3e'
        % (tp, chicken[tp - 1], human[ht - 1], ht, baseline[tp]))
    log('  positions rankable by coupling %d, unrankable %d' % (len(ranked), unrankable))

    V = variants[tp]
    win_positions = sorted({v['pos'] for v in V
                            if v['sep'] <= WINDOW and v['pos'] != tp})
    K_MATCH = len(win_positions)
    log('  the published window of %d admits %d positions, so K = %d is the matched '
        'comparison' % (WINDOW, K_MATCH, K_MATCH))

    def run(admit, policy):
        ok = [v for v in V if v['pos'] in admit and allowed(policy, v['pos'], tp)]
        crossing = [v for v in ok if v['p'] >= THRESHOLD]
        best = max(ok, key=lambda v: v['p']) if ok else None
        return ok, crossing, best

    sweep = []
    ladder = list(KS) + ([K_MATCH] if K_MATCH not in KS else [])
    for K in ladder:
        admit_list = ranked if K is None else ranked[:K]
        admit = {r['pos'] for r in admit_list}
        label = 'all %d' % len(ranked) if K is None else str(K)
        for pol in POLICIES:
            ok, crossing, best = run(admit, pol)
            sweep.append({'rule': 'coupling', 'K': label, 'policy': pol,
                          'positions_admitted': len(admit),
                          'variants_examined': len(ok),
                          'rescuers_crossing': len(crossing),
                          'best_p': best['p'] if best else None,
                          'best_variant': ('%s%d%s' % (best['from'], best['pos'], best['to']))
                                          if best else None,
                          'best_separation': best['sep'] if best else None,
                          'best_fold': best['fold'] if best else None,
                          'matched_to_window': (K == K_MATCH)})
    win = {}
    for pol in POLICIES:
        ok, crossing, best = run(set(win_positions), pol)
        win[pol] = {'rule': 'window', 'radius': WINDOW, 'policy': pol,
                    'positions_admitted': len(win_positions),
                    'variants_examined': len(ok),
                    'rescuers_crossing': len(crossing),
                    'best_p': best['p'] if best else None,
                    'best_variant': ('%s%d%s' % (best['from'], best['pos'], best['to']))
                                    if best else None,
                    'best_separation': best['sep'] if best else None,
                    'best_fold': best['fold'] if best else None}

    log('')
    log('  rule      candidates  policy      admitted  examined  crossing  '
        'best likelihood  best variant      fold')
    for pol in POLICIES:
        w = win[pol]
        log('  window    radius %d    %-10s %9d %9d %9d %16s  %-16s %7s'
            % (WINDOW, pol, w['positions_admitted'], w['variants_examined'],
               w['rescuers_crossing'],
               '%.3e' % w['best_p'] if w['best_p'] else '-',
               w['best_variant'] or '-',
               '%.0f' % w['best_fold'] if w['best_fold'] else '-'))
    for s in sweep:
        if not s['matched_to_window']:
            continue
        log('  coupling  top %-7s %-10s %9d %9d %9d %16s  %-16s %7s'
            % (s['K'], s['policy'], s['positions_admitted'], s['variants_examined'],
               s['rescuers_crossing'],
               '%.3e' % s['best_p'] if s['best_p'] else '-',
               s['best_variant'] or '-',
               '%.0f' % s['best_fold'] if s['best_fold'] else '-'))
    log('')
    log('  coupling-ranked sweep over K')
    log('  K         policy      admitted  examined  crossing  best likelihood  '
        'best variant      separation')
    for s in sweep:
        log('  %-9s %-10s %9d %9d %9d %16s  %-16s %10s'
            % (s['K'], s['policy'], s['positions_admitted'], s['variants_examined'],
               s['rescuers_crossing'],
               '%.3e' % s['best_p'] if s['best_p'] else '-',
               s['best_variant'] or '-',
               s['best_separation'] if s['best_separation'] is not None else '-'))

    all_pub = [s for s in sweep if s['policy'] == 'published' and s['K'].startswith('all')]
    best_overall = None
    cand = [s for s in sweep if s['best_p'] is not None]
    if cand:
        bs = max(cand, key=lambda s: s['best_p'])
        best_overall = {'K': bs['K'], 'policy': bs['policy'],
                        'variant': bs['best_variant'], 'p': bs['best_p'],
                        'separation': bs['best_separation'], 'fold': bs['best_fold']}
    log('')
    if all_pub and all_pub[0]['rescuers_crossing'] == 0:
        log('  no substitution anywhere in the protein lifts this target above %.0e '
            'under the published policy, so the selection rule is not what is '
            'limiting it' % THRESHOLD)
    elif best_overall:
        log('  best admitted rescuer overall: %s reaching %.3e at separation %s'
            % (best_overall['variant'], best_overall['p'], best_overall['separation']))

    top = ranked[:10]
    log('  the ten strongest coupled positions, with what the scan says of them')
    log('  chicken  human  sep   coupling  pct nonlocal  dist A  best substitution  '
        'likelihood')
    for r in top:
        here = [v for v in V if v['pos'] == r['pos']]
        b = max(here, key=lambda v: v['p']) if here else None
        log('  %7d %6d %4d %10.4f %13.1f %7s  %-17s %.3e'
            % (r['pos'], r['human'], r['sep'], r['score'], pct_nl(r['score']),
               ('%.2f' % r['dist']) if r['dist'] is not None else '-',
               ('%s%d%s' % (b['from'], b['pos'], b['to'])) if b else '-',
               b['p'] if b else float('nan')))

    _n = len(top)
    out[tp] = {'human_pos': ht, 'rankable': len(ranked), 'unrankable': unrankable,
               'top_coupled_positions': _n,
               'double_candidate_pairs': _n * (_n - 1) // 2,
               'double_candidate_variants': _n * (_n - 1) // 2 * 19 * 19,
               'baseline': baseline[tp], 'window_positions': win_positions,
               'k_matched_to_window': K_MATCH, 'window': win, 'sweep': sweep,
               'best_overall': best_overall,
               'top_coupled': [dict(
                   {'chicken_pos': r['pos'], 'human_pos': r['human'],
                    'separation': r['sep'], 'score': r['score'],
                    'percentile_nonlocal': pct_nl(r['score']),
                    'distance_A': r['dist'],
                    'in_contact': (r['dist'] is not None
                                   and r['dist'] <= CONTACT_A)},
                   **_best_at(V, r['pos'])) for r in top]}

log('')
log('summary over the five high-risk targets, published policy, matched candidate count')
log('chicken  human  window crossing  coupling crossing  window best  coupling best')
for tp in targets:
    o = out[tp]
    if not o['sweep']:
        log('%7d %6d %17s %18s %12s %14s' % (tp, o['human_pos'], 'not fitted',
                                             'not fitted', '-', '-'))
        continue
    w = o['window']['published']
    cm = [s for s in o['sweep']
          if s['matched_to_window'] and s['policy'] == 'published'][0]
    log('%7d %6d %17d %18d %12s %14s'
        % (tp, o['human_pos'], w['rescuers_crossing'], cm['rescuers_crossing'],
           '%.3e' % w['best_p'] if w['best_p'] else '-',
           '%.3e' % cm['best_p'] if cm['best_p'] else '-'))

json.dump({'source': os.path.basename(SRC), 'scan': SCAN, 'threshold': THRESHOLD,
           'window_radius': WINDOW, 'k_ladder': [('all' if k is None else k) for k in KS],
           'policies': POLICIES, 'contact_cutoff_A': CONTACT_A,
           'targets': {str(k): v for k, v in out.items()}},
          open(OUT, 'w'), indent=1)
log('')
log('wrote ' + OUT)
with open(LOG, 'w') as f:
    f.write('\n'.join(lines) + '\n')
