import csv
import glob
import json
import os
import sys
import time

import numpy as np

HUMAN_LEN = 459
THETA = 0.2
ALPHABET = 'ACDEFGHIKLMNPQRSTVWY-'
Q = len(ALPHABET)
GROUP = 2048
COLBLOCK = 153
OUT = 'data/evcouplings_depth.json'
LOG = 'data/12_evcouplings_depth_log.txt'

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
    hits = sorted(glob.glob('data/evcouplings/**/*.a2m', recursive=True))
    if not hits:
        die('no a2m under data/evcouplings/. Put the EVcouplings download there, '
            'or pass the path to the alignment.')
    if len(hits) > 1:
        die('several a2m files under data/evcouplings/, so name the one to read: %s' % hits)
    return hits[0]


SRC = pick()
if not os.path.exists(SRC):
    die(SRC + ' not found')
log('alignment          : ' + SRC)

stats_hits = sorted(glob.glob(os.path.join(os.path.dirname(SRC), '*_alignment_statistics.csv')))
if len(stats_hits) != 1:
    die('expected one alignment statistics file beside the alignment, found %d: %s'
        % (len(stats_hits), stats_hits))
st_rows = list(csv.DictReader(open(stats_hits[0])))
if not st_rows:
    die(stats_hits[0] + ' has no rows')
ST = st_rows[0]
n_seqs_stated = int(float(ST['num_seqs']))
seqlen_stated = int(float(ST['seqlen']))
num_cov = int(float(ST['num_cov']))
first_uc, last_uc = int(float(ST['1st_uc'])), int(float(ST['last_uc']))
neff_stated = (ST.get('N_eff') or '').strip()
log('statistics         : ' + os.path.basename(stats_hits[0]))
log('  stated          : %d sequences, target length %d, %d covered columns, region %d to %d'
    % (n_seqs_stated, seqlen_stated, num_cov, first_uc, last_uc))
log('  N_eff in the run: %s' % (neff_stated if neff_stated else 'blank, which is why this script exists'))
if seqlen_stated != HUMAN_LEN:
    die('the run targets a %d residue sequence, Q96NT5 is %d' % (seqlen_stated, HUMAN_LEN))

drop = set('abcdefghijklmnopqrstuvwxyz.')
names, seqs = [], []
cur = []
with open(SRC) as f:
    for line in f:
        line = line.rstrip('\n')
        if line.startswith('>'):
            if cur:
                seqs.append(''.join(cur))
            cur = []
            names.append(line[1:])
        elif line:
            cur.append(''.join(ch for ch in line if ch not in drop))
if cur:
    seqs.append(''.join(cur))
if not seqs:
    die(SRC + ' holds no sequences')
if len(seqs) != len(names):
    die('%d headers against %d sequences in %s' % (len(names), len(seqs), SRC))

L = len(seqs[0])
bad = [i for i, s in enumerate(seqs) if len(s) != L]
if bad:
    die('%d sequences differ in length from the first after the insert columns are '
        'removed, the first being row %d at %d against %d'
        % (len(bad), bad[0], len(seqs[bad[0]]), L))
if L != num_cov:
    die('the alignment has %d match states, the statistics file states %d covered '
        'columns' % (L, num_cov))
N = len(seqs)
if N != n_seqs_stated:
    die('the alignment holds %d sequences, the statistics file states %d'
        % (N, n_seqs_stated))
log('sequences          : %d, matching the statistics file' % N)
log('match states       : %d, matching the %d covered columns of %d'
    % (L, num_cov, seqlen_stated))
log('target row         : %s' % names[0][:60])

idx = np.full(256, Q - 1, dtype=np.int8)
for k, ch in enumerate(ALPHABET):
    idx[ord(ch)] = k
raw = np.frombuffer(''.join(seqs).encode('ascii', 'replace'), dtype=np.uint8)
codes = idx[raw].reshape(N, L)
present = sorted({ALPHABET[c] for c in np.unique(codes)})
log('symbols after encoding: %s' % ''.join(present))
gap_frac = float((codes == Q - 1).mean())
log('gap fraction       : %.3f' % gap_frac)

cut = (1.0 - THETA) * L
log('identity cut       : %.1f of %d columns' % (cut, L))
blocks = [(s, min(s + COLBLOCK, L)) for s in range(0, L, COLBLOCK)]
neigh = np.zeros(N, dtype=np.int64)
t0 = time.time()


def onehot(sub):
    b, w = sub.shape
    out = np.zeros((b, w * Q), dtype=np.float32)
    out[np.arange(b)[:, None], np.arange(w)[None, :] * Q + sub] = 1.0
    return out


allhot = [onehot(codes[:, a:b]) for a, b in blocks]
log('one-hot blocks     : %d, %.0f MB' % (len(allhot), sum(x.nbytes for x in allhot) / 1e6))
for g0 in range(0, N, GROUP):
    g1 = min(g0 + GROUP, N)
    counts = np.zeros((g1 - g0, N), dtype=np.float32)
    for (a, b), AH in zip(blocks, allhot):
        counts += onehot(codes[g0:g1, a:b]) @ AH.T
    neigh[g0:g1] = (counts >= cut).sum(axis=1)
    del counts
    log('  rows %6d to %6d done, %.0f s elapsed' % (g0, g1, time.time() - t0))

if (neigh < 1).any():
    die('%d sequences counted fewer than one neighbour, which is impossible since '
        'each is its own neighbour' % int((neigh < 1).sum()))
w = 1.0 / neigh.astype(np.float64)
meff = float(w.sum())
ratio = meff / L
ratio_full = meff / HUMAN_LEN

log('')
log('effective sequences: %.1f' % meff)
log('Meff / L           : %.3f over the %d fitted columns' % (ratio, L))
log('Meff / L           : %.3f over the full %d residues of Q96NT5' % (ratio_full, HUMAN_LEN))
log('neighbour counts   : median %d, max %d, singletons %d'
    % (int(np.median(neigh)), int(neigh.max()), int((neigh == 1).sum())))
log('computed in        : %.0f s' % (time.time() - t0))

prev = None
if os.path.exists('data/dca_diagnostics.json'):
    prev = json.load(open('data/dca_diagnostics.json'))
    log('')
    log('against the shallow alignment of step 3: %d sequences, Meff %.1f, Meff/L %.3f'
        % (prev['n_sequences'], prev['meff'], prev['meff_per_length']))
    log('the step 3 ratio used the full length as its denominator, so %.3f is the '
        'comparable figure' % ratio_full)
    log('fold change in Meff: %.1f' % (meff / prev['meff']))

json.dump({'source': os.path.relpath(SRC), 'n_sequences': N, 'match_states': L,
           'target_length': HUMAN_LEN, 'covered_region': [first_uc, last_uc],
           'n_eff_stated_by_run': neff_stated or None,
           'theta': THETA, 'identity_cut_columns': cut,
           'gap_fraction': gap_frac,
           'meff': meff, 'meff_per_fitted_column': ratio,
           'meff_per_target_length': ratio_full,
           'neighbour_median': int(np.median(neigh)),
           'neighbour_max': int(neigh.max()),
           'singletons': int((neigh == 1).sum()),
           'previous_meff': (prev or {}).get('meff'),
           'previous_meff_per_length': (prev or {}).get('meff_per_length'),
           'previous_n_sequences': (prev or {}).get('n_sequences')},
          open(OUT, 'w'), indent=1)
log('')
log('wrote ' + OUT)
with open(LOG, 'w') as f:
    f.write('\n'.join(lines) + '\n')
