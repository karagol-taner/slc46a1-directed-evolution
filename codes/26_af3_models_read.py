import ast
import json
import os
import re
import statistics
import sys

import gemmi

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
JOBS = os.path.join(ROOT, 'codes', '25_af3_jobs.py')
FASTA = os.path.join(ROOT, 'data', 'construct_sequences.fasta')
ADM = os.path.join(ROOT, 'data', 'admissible_construct.json')
POSSET = os.path.join(ROOT, 'data', 'pharmacophore_position_set.tsv')
AF3 = os.path.join(ROOT, 'data', 'af3')
REFS = {'7BC6': os.path.join(ROOT, 'seq', '7BC6.cif'),
        '7BC7': os.path.join(ROOT, 'seq', '7BC7.cif')}
OUT = os.path.join(ROOT, 'data', 'af3_models_report.json')
LOG = os.path.join(ROOT, 'data', '26_af3_models_read_log.txt')
BASELINE_ARM = 'wt'

lines = []


def log(m):
    print(m, flush=True)
    lines.append(str(m))


def die(m):
    with open(LOG, 'w') as f:
        f.write('\n'.join(lines + ['FATAL ' + m]) + '\n')
    sys.exit('FATAL ' + m)


for p in [JOBS, FASTA, ADM, POSSET] + list(REFS.values()):
    if not os.path.exists(p):
        die('missing input %s' % p)

src = open(JOBS).read()
SEED = next((ast.literal_eval(n.value) for n in ast.parse(src).body
             if isinstance(n, ast.Assign) and getattr(n.targets[0], 'id', None) == 'SEED'),
            None)
if SEED is None:
    die('%s defines no SEED' % os.path.basename(JOBS))

seqs, name = {}, None
for line in open(FASTA):
    line = line.strip()
    if line.startswith('>'):
        name = line[1:].split()[0]
        seqs[name] = []
    elif line and name:
        seqs[name].append(line)
seqs = {k: ''.join(v) for k, v in seqs.items()}
if BASELINE_ARM not in seqs:
    die('%s holds no arm called %s to compare the others against'
        % (os.path.basename(FASTA), BASELINE_ARM))

rows = [l.rstrip('\n').split('\t') for l in open(POSSET)]
pc = {k: i for i, k in enumerate(rows[0])}
PSET = sorted(int(r[pc['chicken_pos']]) for r in rows[1:])
adm = json.load(open(ADM))
SUBS = sorted(set(adm['primaries']) | {int(k) for k in adm['rescuer_positions']})
RESC = sorted(int(k) for k in adm['rescuer_positions'])


def find_root():
    if len(sys.argv) > 1:
        return sys.argv[1]
    cands = set()
    for dirpath, _dirnames, filenames in os.walk(AF3):
        if any(f.endswith('_model_0.cif') for f in filenames):
            cands.add(os.path.dirname(dirpath))
    if not cands:
        die('no folder under %s holds a file ending in _model_0.cif. Unpack the '
            'AlphaFold Server download there, or pass the models root as an argument.'
            % os.path.relpath(AF3, ROOT))
    if len(cands) > 1:
        die('%d candidate model roots under %s: %s. Report rather than pick one.'
            % (len(cands), os.path.relpath(AF3, ROOT),
               sorted(os.path.relpath(c, ROOT) for c in cands)))
    return cands.pop()


MROOT = find_root()
log('models root %s' % os.path.relpath(MROOT, ROOT))
log('seed %d, read from %s' % (SEED, os.path.basename(JOBS)))

folders = {}
for d in sorted(os.listdir(MROOT)):
    full = os.path.join(MROOT, d)
    if os.path.isdir(full) and any(f.endswith('_model_0.cif') for f in os.listdir(full)):
        folders[d] = full
if not folders:
    die('%s holds no subfolder with a model' % MROOT)

matched = {}
for arm in seqs:
    hits = [d for d in folders if re.fullmatch(r'(fold_)?%s' % re.escape(arm), d)]
    if len(hits) != 1:
        die('%d folders under %s match the arm %s: %s. Report rather than pick one.'
            % (len(hits), os.path.relpath(MROOT, ROOT), arm, hits))
    matched[arm] = folders[hits[0]]
log('arms matched: %s'
    % ', '.join('%s to %s' % (a, os.path.basename(p)) for a, p in sorted(matched.items())))


def best_model(folder):
    best = None
    for fn in sorted(os.listdir(folder)):
        m = re.match(r'.*_summary_confidences_(\d+)\.json$', fn)
        if not m:
            continue
        d = json.load(open(os.path.join(folder, fn)))
        if 'ranking_score' not in d:
            die('%s carries no ranking_score' % os.path.join(folder, fn))
        if best is None or d['ranking_score'] > best[1]['ranking_score']:
            best = (int(m.group(1)), d)
    if best is None:
        die('%s holds no summary_confidences file' % folder)
    k, summary = best
    cif = [f for f in os.listdir(folder) if f.endswith('_model_%d.cif' % k)]
    if len(cif) != 1:
        die('%s holds %d files ending in _model_%d.cif' % (folder, len(cif), k))
    return k, summary, os.path.join(folder, cif[0])


def all_models(folder):
    out = []
    for fn in sorted(os.listdir(folder)):
        m = re.match(r'.*_model_(\d+)\.cif$', fn)
        if m:
            out.append((int(m.group(1)), os.path.join(folder, fn)))
    if not out:
        die('%s holds no model file' % folder)
    return sorted(out)


def read_chain(path, require_polymer):
    st = gemmi.read_structure(path)
    st.setup_entities()
    st.remove_waters()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    chain = None
    for ch in st[0]:
        if not require_polymer or len(ch.get_polymer()) > 300:
            chain = ch
            break
    if chain is None:
        die('%s has no chain with more than 300 residues' % path)
    ca, plddt, one = {}, {}, {}
    for r in chain:
        info = gemmi.find_tabulated_residue(r.name)
        a = r.find_atom('CA', '*')
        if a is None or info is None or not info.is_amino_acid():
            continue
        ca[r.seqid.num] = a.pos
        plddt[r.seqid.num] = a.b_iso
        one[r.seqid.num] = info.one_letter_code.upper()
    return ca, plddt, one, chain.name


def rmsd(a, b, over=None):
    shared = sorted(set(a) & set(b) & set(over)) if over else sorted(set(a) & set(b))
    if len(shared) < 200:
        die('only %d alpha carbons in common, too few to superpose' % len(shared))
    sup = gemmi.superpose_positions([a[i] for i in shared], [b[i] for i in shared])
    return sup.rmsd, len(shared)


REF_CA = {}
for tag, path in REFS.items():
    ca, _p, _o, cname = read_chain(path, True)
    REF_CA[tag] = ca
    log('reference %s chain %s, %d alpha carbons' % (tag, cname, len(ca)))
tags = sorted(REF_CA)
CORE = sorted(set(REF_CA[tags[0]]) & set(REF_CA[tags[1]]))
ref_gap, ref_n = rmsd(REF_CA[tags[0]], REF_CA[tags[1]], CORE)
log('every deviation below is measured over the %d positions both references resolve, so '
    'the numbers are on one residue set and the unresolved tails cannot drive them'
    % len(CORE))
log('%s against %s over %d shared alpha carbons: %.3f angstrom, which is the scale every '
    'deviation below is read against' % (tags[0], tags[1], ref_n, ref_gap))

arms, plddts = {}, {}
for arm in sorted(seqs):
    k, summary, cif = best_model(matched[arm])
    ca, pl, one, cname = read_chain(cif, False)
    want = seqs[arm]
    missing = [i + 1 for i in range(len(want)) if i + 1 not in one]
    if missing:
        die('the %s model has no residue at %d of the %d positions the arm was given, '
            'first %s' % (arm, len(missing), len(want), missing[:5]))
    wrong = [(n, one[n], want[n - 1]) for n in sorted(one)
             if n <= len(want) and one[n] != want[n - 1]]
    over = [n for n in one if n > len(want)]
    if wrong or over:
        die('the %s model disagrees with the sequence that arm was given at %d positions, '
            'first %s, and carries %d residues past position %d'
            % (arm, len(wrong), wrong[:3], len(over), len(want)))
    to_ref = {}
    for tag in tags:
        r, n = rmsd(REF_CA[tag], ca, CORE)
        to_ref[tag] = {'rmsd': r, 'n_shared': n}
    tails = sorted(set(ca) - set(CORE))
    arms[arm] = {'folder': os.path.basename(matched[arm]), 'model_index': k,
                 'model_file': os.path.basename(cif), 'chain': cname,
                 'n_residues': len(one),
                 'ranking_score': summary.get('ranking_score'),
                 'ptm': summary.get('ptm'),
                 'fraction_disordered': summary.get('fraction_disordered'),
                 'has_clash': summary.get('has_clash'),
                 'rmsd_to': to_ref,
                 'plddt_chain': statistics.mean(pl.values()),
                 'plddt_set': statistics.mean(pl[p] for p in PSET),
                 'plddt_substituted': statistics.mean(pl[p] for p in SUBS),
                 'plddt_rescuers': statistics.mean(pl[p] for p in RESC),
                 'n_outside_core': len(tails),
                 'plddt_outside_core': (statistics.mean(pl[p] for p in tails)
                                        if tails else None)}
    plddts[arm] = pl
    arms[arm]['_ca'] = ca

pairs = {}
names = sorted(arms)
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        r, n = rmsd(arms[names[i]]['_ca'], arms[names[j]]['_ca'], CORE)
        rw, nw = rmsd(arms[names[i]]['_ca'], arms[names[j]]['_ca'])
        pairs['%s vs %s' % (names[i], names[j])] = {'rmsd': r, 'n_shared': n,
                                                    'rmsd_whole_chain': rw,
                                                    'n_whole_chain': nw}
for a in arms:
    del arms[a]['_ca']

log('')
log('%-11s %-6s %-8s %-6s %-9s %-11s %-11s %s'
    % ('arm', 'model', 'ranking', 'pTM', 'pLDDT', 'RMSD ' + tags[0], 'RMSD ' + tags[1],
       'residues'))
for arm in names:
    a = arms[arm]
    log('%-11s %-6d %-8s %-6s %-9.2f %-11.3f %-11.3f %d'
        % (arm, a['model_index'],
           '-' if a['ranking_score'] is None else '%.4f' % a['ranking_score'],
           '-' if a['ptm'] is None else '%.3f' % a['ptm'], a['plddt_chain'],
           a['rmsd_to'][tags[0]]['rmsd'], a['rmsd_to'][tags[1]]['rmsd'], a['n_residues']))

log('')
log('arm against arm')
for k2, v in sorted(pairs.items()):
    log('  %-26s %.3f angstrom over the %d core positions, %.3f over all %d'
        % (k2, v['rmsd'], v['n_shared'], v['rmsd_whole_chain'], v['n_whole_chain']))
max_pair = max(pairs.values(), key=lambda v: v['rmsd'])['rmsd']
max_ref = max(a['rmsd_to'][t]['rmsd'] for a in arms.values() for t in tags)
log('  furthest arm pair %.3f, furthest arm from a deposited structure %.3f'
    % (max_pair, max_ref))
if max_pair > max_ref:
    die('the arms are further from each other (%.3f) than the furthest of them is from a '
        'deposited structure (%.3f), so they are not one set and a difference between two '
        'of them in simulation would not be about the substitutions.'
        % (max_pair, max_ref))
log('  the arms are closer to each other than to either deposited structure, so they are '
    'comparable with each other')

log('')
log('per-residue confidence, mean pLDDT')
log('%-11s %-10s %-14s %-16s %s'
    % ('arm', 'chain', 'position set', 'substituted', 'rescuer positions'))
for arm in names:
    a = arms[arm]
    log('%-11s %-10.2f %-14.2f %-16.2f %.2f'
        % (arm, a['plddt_chain'], a['plddt_set'], a['plddt_substituted'],
           a['plddt_rescuers']))

base = plddts[BASELINE_ARM]
deltas = {}
for arm in names:
    if arm == BASELINE_ARM:
        continue
    d = [(p, plddts[arm][p] - base[p]) for p in SUBS]
    worst = min(d, key=lambda x: x[1])
    deltas[arm] = {'mean': statistics.mean(v for _p, v in d),
                   'worst_position': worst[0], 'worst_delta': worst[1],
                   'per_position': {str(p): v for p, v in d}}
    log('')
    log('%s against %s at the %d substituted positions: mean change %+.2f, largest fall '
        '%+.2f at chicken %d'
        % (arm, BASELINE_ARM, len(SUBS), deltas[arm]['mean'], worst[1], worst[0]))
    fell = [(p, v) for p, v in d if v <= -5.0]
    if fell:
        log('  positions losing 5 or more: %s'
            % ', '.join('%d (%+.1f)' % (p, v) for p, v in sorted(fell)))
    else:
        log('  no position loses 5 or more against %s' % BASELINE_ARM)

spread = {}
WATCH = sorted({d['worst_position'] for d in deltas.values()}
               | {q for d in deltas.values()
                  for q, v in ((int(k), v) for k, v in d['per_position'].items())
                  if v <= -5.0})
log('')
log('every model of every arm, not only the highest ranking one')
log('%-11s %-7s %-22s %-22s %s'
    % ('arm', 'models', 'pLDDT chain', 'pLDDT substituted', 'RMSD to ' + tags[0]))
for arm in names:
    rows_m = []
    for k, cif in all_models(matched[arm]):
        ca, pl, one, _c = read_chain(cif, False)
        r, _n = rmsd(REF_CA[tags[0]], ca, CORE)
        rows_m.append({'model_index': k,
                       'plddt_chain': statistics.mean(pl.values()),
                       'plddt_set': statistics.mean(pl[p] for p in PSET),
                       'plddt_substituted': statistics.mean(pl[p] for p in SUBS),
                       'rmsd_to_first_reference': r,
                       'watched': {str(q): pl[q] for q in WATCH}})
    spread[arm] = rows_m
    ch = [x['plddt_chain'] for x in rows_m]
    sub = [x['plddt_substituted'] for x in rows_m]
    rr = [x['rmsd_to_first_reference'] for x in rows_m]
    log('%-11s %-7d %-22s %-22s %s'
        % (arm, len(rows_m),
           '%.2f to %.2f, mean %.2f' % (min(ch), max(ch), statistics.mean(ch)),
           '%.2f to %.2f, mean %.2f' % (min(sub), max(sub), statistics.mean(sub)),
           '%.3f to %.3f' % (min(rr), max(rr))))

overlap = {}
log('')
log('the positions where the highest ranking models disagreed, across every model')
for q in WATCH:
    log('  chicken %d' % q)
    rng = {}
    for arm in names:
        vals = [x['watched'][str(q)] for x in spread[arm]]
        rng[arm] = (min(vals), max(vals), statistics.mean(vals))
        log('    %-11s %.1f to %.1f, mean %.1f over %d models'
            % (arm, min(vals), max(vals), statistics.mean(vals), len(spread[arm])))
    lo_b, hi_b, _m = rng[BASELINE_ARM]
    for arm in names:
        if arm == BASELINE_ARM:
            continue
        lo, hi, _mm = rng[arm]
        sep = lo > hi_b or hi < lo_b
        overlap['%d %s' % (q, arm)] = {'position': q, 'arm': arm,
                                       'arm_range': [lo, hi],
                                       'baseline_range': [lo_b, hi_b],
                                       'ranges_separate': sep}
        log('    %s against %s: ranges %s'
            % (arm, BASELINE_ARM, 'do not overlap' if sep else 'overlap'))
sep_any = [k for k, v in overlap.items() if v['ranges_separate']]
log('')
log('positions whose ranges separate from %s across all models: %s'
    % (BASELINE_ARM, sep_any or 'none'))
log('a range that overlaps means the fall seen in the highest ranking pair is inside '
    'model to model variation and carries no weight on its own')


json.dump({'models_root': os.path.relpath(MROOT, ROOT), 'seed': SEED,
           'references': tags,
           'core_positions': len(CORE),
           'reference_gap': {'rmsd': ref_gap, 'n_shared': ref_n,
                             'note': 'the two deposited entries are one state in two '
                                     'crystal forms, so no state is assigned from them'},
           'arms': arms, 'arm_pairs': pairs,
           'furthest_arm_pair': max_pair, 'furthest_arm_to_reference': max_ref,
           'arms_are_one_set': True,
           'baseline_arm': BASELINE_ARM,
           'plddt_deltas_vs_baseline': deltas,
           'all_models': spread, 'watched_positions': WATCH,
           'watched_range_overlap': overlap,
           'watched_separating': sep_any,
           'position_set': PSET, 'substituted': SUBS, 'rescuers': RESC},
          open(OUT, 'w'), indent=1)

log('')
log('wrote ' + os.path.relpath(OUT, ROOT))
with open(LOG, 'w') as f:
    f.write('\n'.join(lines) + '\n')
