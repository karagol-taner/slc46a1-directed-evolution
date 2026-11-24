import csv
import glob
import io
import json
import os
import sys
import zipfile

import numpy as np

POSSET = 'data/pharmacophore_position_set.tsv'
SWEEP = 'data/rescue_sweep_report.json'
OUT = 'data/unrescued_partners.json'
LOG = 'data/14_unrescued_partners_log.txt'
PUBLISHED_RADIUS = 3
CONTACT_A = 5.0
LOCAL = 8
TOPN = 10

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


SRC = pick_source()
for p in (POSSET, SWEEP):
    if not os.path.exists(p):
        die(p + ' not found')


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
need = {'i', 'j', 'cn', 'dist'}
if not need <= set(reader.fieldnames or []):
    die('%s lacks %s; it has %s' % (hits[0], sorted(need), reader.fieldnames))

S, D = {}, {}
for r in reader:
    try:
        i, j, v = int(float(r['i'])), int(float(r['j'])), float(r['cn'])
    except (TypeError, ValueError):
        continue
    k = (min(i, j), max(i, j))
    S[k] = v
    try:
        D[k] = float(r['dist'])
    except (TypeError, ValueError):
        D[k] = float('nan')
if not S:
    die(hits[0] + ' produced no usable rows')
vals = np.array(list(S.values()))
seps = np.array([abs(a - b) for a, b in S])
nl_sorted = np.sort(vals[seps > LOCAL])
covered = sorted({x for k in S for x in k})
lo_cov, hi_cov = covered[0], covered[-1]
log('result             : %s' % SRC)
log('coupled pairs      : %d, positions %d to %d' % (len(S), lo_cov, hi_cov))
log('contact            : at most %.1f Angstrom in the deposited structures' % CONTACT_A)


def pct_nl(v):
    return 100.0 * float(np.searchsorted(nl_sorted, v)) / len(nl_sorted)


rows = [l.rstrip('\n').split('\t') for l in open(POSSET)]
col = {k: i for i, k in enumerate(rows[0])}
by_chicken, by_human = {}, {}
for r in rows[1:]:
    rec = {'human_pos': int(r[col['human_pos']]), 'human_aa': r[col['human_aa']],
           'chicken_pos': int(r[col['chicken_pos']]), 'chicken_aa': r[col['chicken_aa']],
           'differs': r[col['differs']] == '1',
           'parker': r[col['in_parker_named']] == '1',
           'role': r[col['role_or_note']]}
    by_chicken[rec['chicken_pos']] = rec
    by_human[rec['human_pos']] = rec

sweep = json.load(open(SWEEP))
threshold = sweep['threshold']
targets = sorted(int(t) for t in sweep['per_target'])
log('sweep report       : %s, threshold %g, %d high-risk targets'
    % (os.path.basename(SWEEP), threshold, len(targets)))

rescued_any = {}
for s in sweep['summary']:
    for p in s.get('rescued_positions', []):
        rescued_any.setdefault(int(p), []).append((s['radius'], s['policy']))
log('targets rescued at some radius and policy: %s'
    % (sorted(rescued_any) if rescued_any else 'none'))

out = []
for ct in targets:
    if ct not in by_chicken:
        die('chicken %d is a sweep target but is not in %s' % (ct, POSSET))
    rec = by_chicken[ct]
    hp = rec['human_pos']
    was = ct in rescued_any
    log('')
    log('chicken %d %s to human %s, human position %d, baseline likelihood %.3e, %s'
        % (ct, rec['chicken_aa'], rec['human_aa'], hp, sweep['baseline'][str(ct)],
           'rescued' if was else 'not rescued at any radius or policy'))
    if hp < lo_cov or hp > hi_cov:
        log('  human %d lies outside the covered region %d to %d, so the run says '
            'nothing about it' % (hp, lo_cov, hi_cov))
        out.append({'chicken_pos': ct, 'human_pos': hp,
                    'chicken_aa': rec['chicken_aa'], 'human_aa': rec['human_aa'],
                    'baseline_likelihood': sweep['baseline'][str(ct)],
                    'rescued': was, 'rescued_at': rescued_any.get(ct, []),
                    'in_covered_region': False, 'fitted_by_run': False,
                    'n_partners_scored': 0, 'n_partners_in_window': 0,
                    'best_partner_overall': None, 'best_partner_in_window': None,
                    'partners_top': [], 'contacts_beyond_window_above_99th': []})
        continue
    part = []
    for q in range(lo_cov, hi_cov + 1):
        if q == hp:
            continue
        k = (min(hp, q), max(hp, q))
        if k not in S:
            continue
        d = D[k]
        pr = by_human.get(q)
        part.append({'partner': q, 'separation': abs(hp - q), 'score': S[k],
                     'percentile_nonlocal': pct_nl(S[k]),
                     'distance_A': d if d == d else None,
                     'in_contact': bool(d <= CONTACT_A) if d == d else None,
                     'within_published_window': abs(hp - q) <= PUBLISHED_RADIUS,
                     'in_position_set': pr is not None,
                     'is_primary': bool(pr and pr['differs']),
                     'parker_named': bool(pr and pr['parker']),
                     'role': (pr or {}).get('role', '')})
    part.sort(key=lambda r: -r['score'])
    if not part:
        log('  human %d lies inside the covered region but is one of the columns the '
            'run does not fit, so it has no partners to report' % hp)
        out.append({'chicken_pos': ct, 'human_pos': hp,
                    'chicken_aa': rec['chicken_aa'], 'human_aa': rec['human_aa'],
                    'baseline_likelihood': sweep['baseline'][str(ct)],
                    'rescued': was, 'rescued_at': rescued_any.get(ct, []),
                    'in_covered_region': True, 'fitted_by_run': False,
                    'n_partners_scored': 0, 'n_partners_in_window': 0,
                    'best_partner_overall': None, 'best_partner_in_window': None,
                    'partners_top': [], 'contacts_beyond_window_above_99th': []})
        continue
    top = part[:TOPN]
    inside = [r for r in part if r['within_published_window']]
    best_in = max(inside, key=lambda r: r['score']) if inside else None
    best_any = part[0] if part else None
    log('  partners scored: %d, of which %d lie inside the published window of %d'
        % (len(part), len(inside), PUBLISHED_RADIUS))
    log('  partner  sep    score  pct nonlocal  dist A  contact  in window  note')
    for r in top:
        dt = '%6.2f' % r['distance_A'] if r['distance_A'] is not None else '     -'
        note = ''
        if r['is_primary']:
            note = 'primary'
        elif r['in_position_set']:
            note = 'in the set'
        if r['role']:
            note = (note + ', ' if note else '') + r['role'][:34]
        log('  %7d %4d %8.4f %13.1f %s %8s %11s  %s'
            % (r['partner'], r['separation'], r['score'], r['percentile_nonlocal'],
               dt, 'yes' if r['in_contact'] else 'no',
               'yes' if r['within_published_window'] else 'no', note))
    if best_any and best_in:
        log('  strongest partner overall is %d at separation %d, scoring %.4f; the '
            'strongest the window reaches is %d at separation %d, scoring %.4f'
            % (best_any['partner'], best_any['separation'], best_any['score'],
               best_in['partner'], best_in['separation'], best_in['score']))
    contacts_outside = [r for r in part
                        if r['in_contact'] and not r['within_published_window']
                        and r['percentile_nonlocal'] >= 99.0]
    log('  partners above the 99th percentile, in contact, and beyond the window: %d %s'
        % (len(contacts_outside), [r['partner'] for r in contacts_outside]))
    out.append({'chicken_pos': ct, 'human_pos': hp,
                'chicken_aa': rec['chicken_aa'], 'human_aa': rec['human_aa'],
                'baseline_likelihood': sweep['baseline'][str(ct)],
                'rescued': was, 'rescued_at': rescued_any.get(ct, []),
                'in_covered_region': True, 'fitted_by_run': True,
                'n_partners_scored': len(part),
                'n_partners_in_window': len(inside),
                'best_partner_overall': best_any,
                'best_partner_in_window': best_in,
                'partners_top': top,
                'contacts_beyond_window_above_99th': contacts_outside})

log('')
unres = [r for r in out if not r['rescued']]
res = [r for r in out if r['rescued']]
log('the %d high-risk targets the sweep did not rescue:' % len(unres))
for r in unres:
    b, w = r['best_partner_overall'], r['best_partner_in_window']
    if b is None:
        log('  human %d: the run does not fit this column, so it has no partners'
            % r['human_pos'])
        continue
    log('  human %d: strongest partner %d at separation %d and %.2f Angstrom, '
        'against %d at separation %d inside the window of %d'
        % (r['human_pos'], b['partner'], b['separation'],
           b['distance_A'] if b['distance_A'] is not None else float('nan'),
           w['partner'] if w else -1, w['separation'] if w else -1,
           PUBLISHED_RADIUS))
    log('    %d partners above the 99th percentile sit in contact and beyond the '
        'window: %s' % (len(r['contacts_beyond_window_above_99th']),
                        [x['partner'] for x in r['contacts_beyond_window_above_99th']]))
log('')
log('the %d the sweep did rescue, for comparison:' % len(res))
for r in res:
    b, w = r['best_partner_overall'], r['best_partner_in_window']
    if b is None:
        continue
    log('  human %d: strongest partner %d at separation %d, and the window reaches '
        '%d at separation %d scoring %.4f against %.4f'
        % (r['human_pos'], b['partner'], b['separation'],
           w['partner'] if w else -1, w['separation'] if w else -1,
           w['score'] if w else float('nan'), b['score']))
gap = [r for r in out if r['best_partner_overall'] and r['best_partner_in_window']]
if gap:
    log('')
    log('separation of the strongest partner, by whether the sweep rescued the target')
    for r in gap:
        log('  human %3d  %-12s strongest partner at separation %3d, window best at %d'
            % (r['human_pos'], 'rescued' if r['rescued'] else 'not rescued',
               r['best_partner_overall']['separation'],
               r['best_partner_in_window']['separation']))

json.dump({'source': os.path.basename(SRC), 'sweep_threshold': threshold,
           'published_radius': PUBLISHED_RADIUS, 'contact_cutoff_A': CONTACT_A,
           'local_cutoff': LOCAL, 'top_n_reported': TOPN,
           'covered_region': [lo_cov, hi_cov], 'targets': out},
          open(OUT, 'w'), indent=1)
log('')
log('wrote ' + OUT)
with open(LOG, 'w') as f:
    f.write('\n'.join(lines) + '\n')
