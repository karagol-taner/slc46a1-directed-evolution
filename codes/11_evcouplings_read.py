import csv
import glob
import io
import json
import os
import statistics
import sys
import zipfile

import numpy as np

HUMAN_LEN = 459
LOCAL = 8
MEFF_PER_L_USABLE = 1.0
SCORE_COL = 'cn'
CONTACT_A = 5.0
DEPTH = 'data/evcouplings_depth.json'
PREV = 'data/dca_diagnostics.json'


def die(m):
    sys.exit('FATAL ' + m)


def pick_source():
    if len(sys.argv) > 1:
        return sys.argv[1]
    zips = sorted(glob.glob('data/evcouplings/*.zip'))
    if len(zips) == 1:
        return zips[0]
    if not zips:
        die('no result given and no zip in data/evcouplings/. Put the EVcouplings '
            'download there, or pass its path.')
    die('several zips in data/evcouplings/, so name the one to read: %s'
        % [os.path.basename(z) for z in zips])


SRC = pick_source()
if not os.path.exists(SRC):
    die(SRC + ' not found')


def members(src):
    if os.path.isdir(src):
        out = {}
        for root, _, files in os.walk(src):
            for f in files:
                p = os.path.join(root, f)
                out[os.path.relpath(p, src)] = ('file', p)
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
def find(suffix, exclude=()):
    hits = [n for n in M if n.endswith(suffix) and '__MACOSX' not in n
            and not any(e in n for e in exclude)]
    if not hits:
        die('no file ending %s in %s' % (suffix, SRC))
    if len(hits) > 1:
        die('several files ending %s: %s' % (suffix, hits))
    return hits[0]


stats_name = find('_alignment_statistics.csv')
rows = list(csv.DictReader(io.StringIO(read_text(SRC, M[stats_name]))))
if not rows:
    die(stats_name + ' has no rows')
st = rows[0]
num_seqs = int(float(st['num_seqs']))
seqlen = int(float(st['seqlen']))
neff_raw = (st.get('N_eff') or '').strip()
first_uc, last_uc = int(float(st['1st_uc'])), int(float(st['last_uc']))
if seqlen != HUMAN_LEN:
    die('the run targets a %d residue sequence, Q96NT5 is %d' % (seqlen, HUMAN_LEN))

num_cov = int(float(st['num_cov']))
if neff_raw:
    meff, meff_src = float(neff_raw), 'N_eff as EVcouplings reports it'
elif os.path.exists(DEPTH):
    dep = json.load(open(DEPTH))
    if dep.get('n_sequences') != num_seqs or dep.get('match_states') != num_cov:
        die('%s was computed on a different alignment: %s sequences and %s match '
            'states against %d and %d here'
            % (DEPTH, dep.get('n_sequences'), dep.get('match_states'), num_seqs, num_cov))
    meff = float(dep['meff'])
    meff_src = ('computed by codes/12_evcouplings_depth.py at %g identity, '
                'because N_eff is blank in this run' % (1.0 - dep['theta']))
else:
    die('N_eff is blank in %s and %s does not exist, so there is no effective '
        'sequence count to report. Run codes/12_evcouplings_depth.py first.'
        % (stats_name, DEPTH))
ratio = meff / seqlen
ratio_fitted = meff / num_cov
print('run                : %s' % os.path.basename(SRC))
print('sequences          : %d' % num_seqs)
print('effective sequences: %.1f (%s)' % (meff, meff_src))
print('target length      : %d, covered region %d to %d, %d columns fitted'
      % (seqlen, first_uc, last_uc, num_cov))
print('Meff / L           : %.3f over the full length, %.3f over the fitted columns '
      '(usable from about %.1f)' % (ratio, ratio_fitted, MEFF_PER_L_USABLE))
print('verdict            : %s'
      % ('underdetermined' if ratio < MEFF_PER_L_USABLE else 'usable'))

cs_name = find('_CouplingScores.csv', exclude=('_longrange',))
reader = csv.DictReader(io.StringIO(read_text(SRC, M[cs_name])))
need = {'i', 'j', SCORE_COL}
if not need <= set(reader.fieldnames or []):
    die('%s lacks %s; it has %s' % (cs_name, sorted(need), reader.fieldnames))

S = {}
vals = []
seps = []
for r in reader:
    try:
        i, j, v = int(float(r['i'])), int(float(r['j'])), float(r[SCORE_COL])
    except (TypeError, ValueError):
        continue
    key = (min(i, j), max(i, j))
    S[key] = v
    vals.append(v)
    seps.append(abs(i - j))
if not S:
    die(cs_name + ' produced no usable rows')
vals = np.array(vals)
seps = np.array(seps)
srt = np.sort(vals)
print('coupled pairs      : %d' % len(S))

order = np.argsort(-vals)
print('\nseparation profile of the strongest couplings')
print(' top N   median sep   fraction with sep <= %d' % LOCAL)
for n in (20, 50, 100, 200, 500):
    k = order[:min(n, len(vals))]
    print('%6d %12.0f %22.2f' % (n, statistics.median(seps[k]),
                                 float((seps[k] <= LOCAL).mean())))
print('all %d pairs: median sep %.0f, fraction local %.2f'
      % (len(vals), statistics.median(seps), float((seps <= LOCAL).mean())))

cmp_name = find('_CouplingScoresCompared_all.csv')
creader = csv.DictReader(io.StringIO(read_text(SRC, M[cmp_name])))
cneed = {'i', 'j', SCORE_COL, 'dist'}
if not cneed <= set(creader.fieldnames or []):
    die('%s lacks %s; it has %s' % (cmp_name, sorted(cneed), creader.fieldnames))
D = {}
cmp_rows = []
for r in creader:
    try:
        i, j, v = int(float(r['i'])), int(float(r['j'])), float(r[SCORE_COL])
    except (TypeError, ValueError):
        continue
    try:
        dd = float(r['dist'])
    except (TypeError, ValueError):
        dd = float('nan')
    key = (min(i, j), max(i, j))
    D[key] = dd
    cmp_rows.append((v, abs(i - j), dd))
if not cmp_rows:
    die(cmp_name + ' produced no usable rows')

sh_name = find('_structure_hits.csv')
sh = list(csv.DictReader(io.StringIO(read_text(SRC, M[sh_name]))))
pdbs = sorted({(r.get('pdb_id') or '').lower() for r in sh if r.get('pdb_id')})
coord = []
for r in sh:
    try:
        coord.append((int(float(r['coord_start'])), int(float(r['coord_end']))))
    except (TypeError, ValueError, KeyError):
        pass
if len(set(coord)) != 1:
    die('%s gives %d different coordinate ranges, %s, so no single range can be '
        'reported' % (sh_name, len(set(coord)), sorted(set(coord))))
coord_range = list(coord[0])
print('\nstructures matched by the run: %s, over residues %s to %s'
      % (', '.join(pdbs) or 'none',
         sh[0].get('coord_start') if sh else '?', sh[0].get('coord_end') if sh else '?'))
if not pdbs:
    die(sh_name + ' names no structure, so no structural readout can be made')

cmp_rows.sort(key=lambda t: -t[0])
nl = [t for t in cmp_rows if t[1] > LOCAL and t[2] == t[2]]
print('non-local couplings with a distance in the structure: %d of %d'
      % (len(nl), len(cmp_rows)))
print('\nprecision of the strongest non-local couplings against the structure, '
      'contact at %.1f Angstrom' % CONTACT_A)
print('  rank cut        N   precision   median distance')
prec_rows = []
for label, n in (('L/2', seqlen // 2), ('L', seqlen), ('2L', 2 * seqlen)):
    k = nl[:min(n, len(nl))]
    if not k:
        continue
    pr = float(np.mean([t[2] <= CONTACT_A for t in k]))
    md = float(np.median([t[2] for t in k]))
    print('%9s %8d %11.3f %17.2f' % (label, len(k), pr, md))
    prec_rows.append({'cut': label, 'n': len(k), 'precision': pr,
                      'median_distance_A': md})
base = float(np.mean([t[2] <= CONTACT_A for t in nl]))
print('all %d non-local pairs: %.3f are in contact, which is the rate to beat'
      % (len(nl), base))

nonlocal_sorted = np.sort(vals[seps > LOCAL])

def look(a, b):
    lo, hi = min(a, b), max(a, b)
    if lo < first_uc or hi > last_uc:
        return None, None, None, 'outside the covered region %d to %d' % (first_uc, last_uc)
    if (lo, hi) not in S:
        return None, None, None, 'not in the coupling table'
    v = S[(lo, hi)]
    q = 100.0 * float(np.searchsorted(srt, v)) / len(srt)
    qn = 100.0 * float(np.searchsorted(nonlocal_sorted, v)) / len(nonlocal_sorted)
    return v, q, qn, ''


PAIRS = [('cross-link N90C/S407C', 90, 407),
         ('cross-link Q45C/L290C', 45, 290),
         ('cross-link L290C/N411C', 290, 411),
         ('F392V with S196L', 392, 196),
         ('F392D with G305L', 392, 305),
         ('D109A with G305L', 109, 305)]

print('\nmeasured pairs')
print('pair                       sep     score  pct all  pct nonlocal   dist A  note')
out_pairs, usable_q = [], []
for lab, a, b in PAIRS:
    v, q, qn, note = look(a, b)
    sep = abs(a - b)
    dd = D.get((min(a, b), max(a, b)), float('nan'))
    dtxt = '%7.2f' % dd if dd == dd else '      -'
    if v is None:
        print('%-26s %4d %9s %8s %13s %s  %s' % (lab, sep, '-', '-', '-', dtxt, note))
    else:
        print('%-26s %4d %9.3f %8.1f %13.1f %s' % (lab, sep, v, q, qn, dtxt))
        usable_q.append(qn)
    out_pairs.append({'label': lab, 'i': a, 'j': b, 'separation': sep,
                      'score': v, 'percentile_all': q, 'percentile_nonlocal': qn,
                      'distance_A': dd if dd == dd else None,
                      'in_contact': bool(dd <= CONTACT_A) if dd == dd else None,
                      'note': note})

if usable_q:
    print('\nmedian percentile of the measured pairs against non-local pairs: %.1f'
          % statistics.median(usable_q))
    print('a set with no detectable enrichment sits near 50')
else:
    print('\nno measured pair could be scored, so no enrichment statement is made')

pos = [l.rstrip('\n').split('\t') for l in open('data/pharmacophore_position_set.tsv')]
pcol = {k: i for i, k in enumerate(pos[0])}
prim = sorted(int(r[pcol['human_pos']]) for r in pos[1:] if r[pcol['differs']] == '1')
prim_rows = []
for x in range(len(prim)):
    for y in range(x + 1, len(prim)):
        v, q, qn, note = look(prim[x], prim[y])
        dd = D.get((prim[x], prim[y]), float('nan'))
        prim_rows.append({'i': prim[x], 'j': prim[y],
                          'separation': abs(prim[x] - prim[y]),
                          'score': v, 'percentile_nonlocal': qn,
                          'distance_A': dd if dd == dd else None,
                          'in_contact': bool(dd <= CONTACT_A) if dd == dd else None,
                          'note': note})
scored = [r for r in prim_rows if r['score'] is not None]
nonloc = [r for r in scored if r['separation'] > LOCAL]
print('\nprimary pairs: %d total, %d scored, %d of those separated by more than %d'
      % (len(prim_rows), len(scored), len(nonloc), LOCAL))
if nonloc:
    top = sorted(nonloc, key=lambda r: -r['score'])[:5]
    print('strongest non-local primary pairs')
    for r in top:
        dtxt = '%6.2f' % r['distance_A'] if r['distance_A'] is not None else '     -'
        print('  %3d-%3d  sep %3d  score %7.3f  pct nonlocal %5.1f  dist %s A'
              % (r['i'], r['j'], r['separation'], r['score'],
                 r['percentile_nonlocal'], dtxt))
    withd = [r for r in nonloc if r['distance_A'] is not None]
    if withd:
        print('%d of the %d scored non-local primary pairs have a distance, and %d '
              'are within %.1f Angstrom'
              % (len(withd), len(nonloc),
                 sum(1 for r in withd if r['in_contact']), CONTACT_A))

fold = None
if os.path.exists(PREV):
    pv = json.load(open(PREV))
    if pv.get('meff'):
        fold = meff / float(pv['meff'])
        print('\nagainst the step 3 alignment: %d sequences and Meff %.1f, so the '
              'effective count rises %.1f-fold' % (pv['n_sequences'], pv['meff'], fold))

json.dump({'source': os.path.basename(SRC), 'score_column': SCORE_COL,
           'previous_meff': (json.load(open(PREV))['meff']
                             if os.path.exists(PREV) else None),
           'meff_fold_change_over_step3': fold,
           'num_seqs': num_seqs, 'meff': meff, 'meff_source': meff_src,
           'length': seqlen, 'covered_region': [first_uc, last_uc],
           'fitted_columns': num_cov,
           'meff_per_length': ratio, 'meff_per_fitted_column': ratio_fitted,
           'structures_matched': pdbs, 'structure_coord_range': coord_range,
           'contact_cutoff_A': CONTACT_A,
           'nonlocal_with_distance': len(nl),
           'baseline_contact_rate_nonlocal': base,
           'precision_against_structure': prec_rows,
           'underdetermined': bool(ratio < MEFF_PER_L_USABLE),
           'local_cutoff': LOCAL, 'n_coupled_pairs': len(S),
           'measured_pairs': out_pairs,
           'median_percentile_measured_nonlocal':
               statistics.median(usable_q) if usable_q else None,
           'primary_pairs': prim_rows},
          open('data/evcouplings_report.json', 'w'), indent=1)
print('\nwrote data/evcouplings_report.json')
