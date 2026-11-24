import json
import os
import statistics
import sys

import numpy as np

REPORT = 'data/dca_step3/coupling_report.json'
MATRIX = 'data/dca_step3/coupling_matrix.npy'
MEFF_PER_L_USABLE = 1.0
LOCAL = 8


def die(m):
    sys.exit('FATAL ' + m)


for p in (REPORT, MATRIX):
    if not os.path.exists(p):
        die(p + ' not found')

rep = json.load(open(REPORT))
S = np.load(MATRIX)
L = rep['length']
if S.shape != (L, L):
    die('matrix is %s but the report says length %d' % (S.shape, L))

meff, nseq = rep['meff'], rep['n_sequences']
ratio = meff / L
print('alignment: %d sequences, Meff %.1f, length %d' % (nseq, meff, L))
print('Meff / L = %.3f (usable from about %.1f)' % (ratio, MEFF_PER_L_USABLE))
print('verdict: %s' % ('underdetermined' if ratio < MEFF_PER_L_USABLE else 'usable'))

iu = np.triu_indices(L, 1)
vals, sep = S[iu], np.abs(iu[0] - iu[1])
order = np.argsort(-vals)
print('\nseparation profile of the strongest couplings')
print(' top N   median sep   fraction with sep <= %d' % LOCAL)
for n in (20, 50, 100, 200, 500):
    k = order[:n]
    print('%6d %12.0f %22.2f' % (n, statistics.median(sep[k]), float((sep[k] <= LOCAL).mean())))
print('for reference, all %d pairs: median sep %.0f, fraction local %.2f'
      % (len(vals), statistics.median(sep), float((sep <= LOCAL).mean())))

def pct_of(a, b):
    v = float(S[a - 1, b - 1])
    return v, 100.0 * float(np.searchsorted(np.sort(vals), v)) / len(vals)

print('\nmeasured pairs against the whole distribution')
rows = []
for g in rep['gate_pairs']:
    if g['measured']:
        v, q = pct_of(g['i'], g['j'])
        rows.append(('cross-link %d/%d' % (g['i'], g['j']), g['sep'], v, q))
for c in rep['compensatory_pairs']:
    v, q = pct_of(c['i'], c['j'])
    rows.append((c['label'], c['sep'], v, q))
print('pair                      sep    score  percentile')
for lab, s, v, q in rows:
    print('%-24s %5d %8.3f %11.1f' % (lab, s, v, q))

meas_q = [q for _, _, _, q in rows]
nonlocal_mask = sep > LOCAL
nonlocal_vals = np.sort(vals[nonlocal_mask])
def pct_nonlocal(v):
    return 100.0 * float(np.searchsorted(nonlocal_vals, v)) / len(nonlocal_vals)
meas_q_nl = [pct_nonlocal(v) for _, _, v, _ in rows]
print('\nmedian percentile of the measured pairs: %.1f against all pairs, '
      '%.1f against pairs separated by more than %d'
      % (statistics.median(meas_q), statistics.median(meas_q_nl), LOCAL))
print('a set with no detectable enrichment sits near 50')

out = {'n_sequences': nseq, 'meff': meff, 'length': L, 'meff_per_length': ratio,
       'usable_threshold': MEFF_PER_L_USABLE,
       'underdetermined': bool(ratio < MEFF_PER_L_USABLE),
       'local_cutoff': LOCAL,
       'separation_profile': {str(n): {
           'median_separation': float(statistics.median(sep[order[:n]])),
           'fraction_local': float((sep[order[:n]] <= LOCAL).mean())} for n in
           (20, 50, 100, 200, 500)},
       'measured_pairs': [{'label': lab, 'separation': s, 'score': v,
                           'percentile_all': q, 'percentile_nonlocal': qn}
                          for (lab, s, v, q), qn in zip(rows, meas_q_nl)],
       'median_percentile_measured_all': statistics.median(meas_q),
       'median_percentile_measured_nonlocal': statistics.median(meas_q_nl)}
json.dump(out, open('data/dca_diagnostics.json', 'w'), indent=1)
print('\nwrote data/dca_diagnostics.json')
