import csv
import glob
import io
import json
import os
import statistics
import sys
import zipfile

import numpy as np

POSSET = 'data/pharmacophore_position_set.tsv'
OUT = 'data/evmutation_primaries.json'
LOG = 'data/13_evmutation_primaries_log.txt'
EPI = 'prediction_epistatic'
IND = 'prediction_independent'

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


def find(suffix):
    hits = [n for n in M if n.endswith(suffix) and '__MACOSX' not in n]
    if not hits:
        die('no file ending %s in %s' % (suffix, SRC))
    if len(hits) > 1:
        die('several files ending %s: %s' % (suffix, hits))
    return hits[0]


log('result             : ' + SRC)

mm_name = find('_single_mutant_matrix.csv')
reader = csv.DictReader(io.StringIO(read_text(SRC, M[mm_name])))
need = {'pos', 'wt', 'subs', EPI, IND}
if not need <= set(reader.fieldnames or []):
    die('%s lacks %s; it has %s' % (mm_name, sorted(need), reader.fieldnames))

WT, SC, IC = {}, {}, {}
all_epi = []
for r in reader:
    try:
        pos = int(float(r['pos']))
        epi = float(r[EPI])
        ind = float(r[IND])
    except (TypeError, ValueError):
        continue
    wt, sub = (r['wt'] or '').strip(), (r['subs'] or '').strip()
    if not wt or not sub:
        continue
    if WT.setdefault(pos, wt) != wt:
        die('%s gives two wild-type residues at position %d, %s and %s'
            % (mm_name, pos, WT[pos], wt))
    SC[(pos, sub)] = epi
    IC[(pos, sub)] = ind
    all_epi.append(epi)
if not SC:
    die(mm_name + ' produced no usable rows')
all_epi = np.array(all_epi)
srt = np.sort(all_epi)
positions = sorted(WT)
log('matrix             : %s' % os.path.basename(mm_name))
log('scored positions   : %d, from %d to %d' % (len(positions), positions[0], positions[-1]))
log('scored substitutions: %d' % len(SC))
log('epistatic score    : median %.3f, 5th and 95th percentile %.3f and %.3f'
    % (float(np.median(all_epi)), float(np.percentile(all_epi, 5)),
       float(np.percentile(all_epi, 95))))

en_name = find('_enrichment.csv')
ereader = csv.DictReader(io.StringIO(read_text(SRC, M[en_name])))
if not {'i', 'enrichment'} <= set(ereader.fieldnames or []):
    die('%s lacks i and enrichment; it has %s' % (en_name, ereader.fieldnames))
EN = {}
for r in ereader:
    try:
        EN[int(float(r['i']))] = float(r['enrichment'])
    except (TypeError, ValueError):
        continue
if not EN:
    die(en_name + ' produced no usable rows')
en_sorted = np.sort(np.array(list(EN.values())))
log('enrichment         : %d positions, median %.3f' % (len(EN), float(np.median(en_sorted))))

rows = [l.rstrip('\n').split('\t') for l in open(POSSET)]
col = {k: i for i, k in enumerate(rows[0])}
pset = [{'human_pos': int(r[col['human_pos']]), 'human_aa': r[col['human_aa']],
         'chicken_pos': int(r[col['chicken_pos']]), 'chicken_aa': r[col['chicken_aa']],
         'differs': r[col['differs']] == '1',
         'parker': r[col['in_parker_named']] == '1',
         'role': r[col['role_or_note']]} for r in rows[1:]]
prim = [p for p in pset if p['differs']]
log('position set       : %d positions, %d of them primary' % (len(pset), len(prim)))

mismatch = [(p['human_pos'], p['human_aa'], WT[p['human_pos']])
            for p in pset if p['human_pos'] in WT and WT[p['human_pos']] != p['human_aa']]
checked = sum(1 for p in pset if p['human_pos'] in WT)
if mismatch:
    die('the matrix wild type disagrees with the position set at %d positions, the '
        'first being human %d where the set has %s and the matrix has %s'
        % (len(mismatch), mismatch[0][0], mismatch[0][1], mismatch[0][2]))
log('numbering check    : %d of %d positions are in the matrix and every wild-type '
    'residue agrees' % (checked, len(pset)))
absent = [p['human_pos'] for p in pset if p['human_pos'] not in WT]
log('positions the run does not score: %d %s' % (len(absent), absent))


def pct(sorted_vals, v):
    return 100.0 * float(np.searchsorted(sorted_vals, v)) / len(sorted_vals)


log('')
log('the chicken residue scored as a substitution in the human background')
log('human  human  chicken  epistatic  pct of all  pct at the  independent  '
    'EC enrichment  pct')
log('  pos     aa       aa      score         subs    position        score'
    '                       ')
out = []
scored = []
for p in sorted(prim, key=lambda d: d['human_pos']):
    hp, ha, ca = p['human_pos'], p['human_aa'], p['chicken_aa']
    key = (hp, ca)
    epi = SC.get(key)
    ind = IC.get(key)
    en = EN.get(hp)
    at_pos = [v for (q, s), v in SC.items() if q == hp]
    rec = {'human_pos': hp, 'human_aa': ha, 'chicken_pos': p['chicken_pos'],
           'chicken_aa': ca, 'parker_named': p['parker'], 'role': p['role'],
           'epistatic': epi, 'independent': ind,
           'percentile_all': pct(srt, epi) if epi is not None else None,
           'percentile_at_position': (100.0 * float(sum(1 for v in at_pos if v < epi))
                                      / len(at_pos)) if (epi is not None and at_pos) else None,
           'ec_enrichment': en,
           'ec_enrichment_percentile': pct(en_sorted, en) if en is not None else None}
    out.append(rec)
    if epi is None:
        log('%5d %6s %8s %10s %12s %11s %12s %14s %4s'
            % (hp, ha, ca, 'not scored', '-', '-', '-',
               ('%.2f' % en) if en is not None else '-', '-'))
        continue
    scored.append(rec)
    log('%5d %6s %8s %10.3f %12.1f %11.1f %12.3f %14s %4s'
        % (hp, ha, ca, epi, rec['percentile_all'], rec['percentile_at_position'], ind,
           ('%.2f' % en) if en is not None else '-',
           ('%.0f' % rec['ec_enrichment_percentile']) if en is not None else '-'))

log('')
if scored:
    es = [r['epistatic'] for r in scored]
    ps = [r['percentile_all'] for r in scored]
    log('%d of the %d primaries are scored by the run' % (len(scored), len(prim)))
    log('median epistatic score of the chicken residue in human: %.3f, against %.3f '
        'for every substitution in the matrix'
        % (statistics.median(es), float(np.median(all_epi))))
    log('median percentile of the primaries among all substitutions: %.1f, where a '
        'set with no detectable difference sits near 50' % statistics.median(ps))
    worst = sorted(scored, key=lambda r: r['epistatic'])[:5]
    log('the five chicken residues human tolerates least')
    for r in worst:
        log('  human %d %s to %s, epistatic %.3f at percentile %.1f%s'
            % (r['human_pos'], r['human_aa'], r['chicken_aa'], r['epistatic'],
               r['percentile_all'], ', named by Parker' if r['parker_named'] else ''))
    ens = [r['ec_enrichment_percentile'] for r in scored
           if r['ec_enrichment_percentile'] is not None]
    if ens:
        log('median coupling enrichment percentile of the primaries: %.1f over %d '
            'positions' % (statistics.median(ens), len(ens)))
else:
    log('no primary could be scored, so no statement is made about the set')

json.dump({'source': os.path.basename(SRC), 'matrix': os.path.basename(mm_name),
           'n_positions_scored': len(positions),
           'n_substitutions_scored': len(SC),
           'matrix_median_epistatic': float(np.median(all_epi)),
           'numbering_checked_positions': checked,
           'numbering_mismatches': len(mismatch),
           'positions_not_scored': absent,
           'n_primaries': len(prim), 'n_primaries_scored': len(scored),
           'median_epistatic_primaries':
               statistics.median([r['epistatic'] for r in scored]) if scored else None,
           'median_percentile_primaries':
               statistics.median([r['percentile_all'] for r in scored]) if scored else None,
           'primaries': out},
          open(OUT, 'w'), indent=1)
log('')
log('wrote ' + OUT)
with open(LOG, 'w') as f:
    f.write('\n'.join(lines) + '\n')
