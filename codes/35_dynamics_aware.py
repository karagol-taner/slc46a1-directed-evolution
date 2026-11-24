import csv
import glob
import json
import os
import sys

ARMS = ('wt', 'primaries', 'construct')
MD = 'data/md'
EVDIR = 'data/evcouplings'
POSSET = 'data/pharmacophore_position_set.tsv'
CONSTRUCT = 'data/admissible_construct.json'
MAP = 'structures/7BC6_to_human.json'
HUMAN_FASTA = 'seq/Q96NT5.fasta'
CHICKEN_FASTA = 'seq/7BC6.fasta'
OUTDIR = 'data/dynamics_aware'
OUT = 'data/dynamics_aware_summary.json'
LOG = 'data/35_dynamics_aware_log.txt'

REPLICAS_IN_PUBLISHED_METHOD = 3
TOP_FRACTION = 0.10
PUBLISHED_DCS_CUT = 0.95
PUBLISHED_RCS_CUT = 6.0
NM_TO_A = 10.0
TOL_A = 0.0006

lines = []


def log(m):
    print(m, flush=True)
    lines.append(str(m))


def die(m):
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    with open(LOG, 'w') as f:
        f.write('\n'.join(lines + ['FATAL ' + m]) + '\n')
    sys.exit('FATAL ' + m)


def need(path, what):
    if not os.path.isfile(path):
        die('no %s, so %s' % (path, what))
    return path


def fasta_one(path):
    return ''.join(l.strip() for l in open(path) if not l.startswith('>'))


def fasta_entity(path, key):
    ent = {}
    for rec in open(path).read().split('>'):
        if not rec.strip():
            continue
        head, *rest = rec.strip().split('\n')
        ent[head.split('|')[0]] = ''.join(rest)
    if key not in ent:
        die('%s carries no entity %s, it has %s'
            % (path, key, ', '.join(sorted(ent))))
    return ent[key]


def tsv(path):
    with open(path) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        return [dict(zip(head, L.rstrip('\n').split('\t')))
                for L in fh if L.strip()]


def minmax(vals):
    lo, hi = min(vals), max(vals)
    span = hi - lo
    return [0.0 for _ in vals] if span == 0 else [(v - lo) / span for v in vals]


def pearson(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    return 0.0 if sxx == 0 or syy == 0 else sxy / (sxx * syy) ** 0.5


def quantiles(vals):
    s = sorted(vals)
    n = len(s)

    def q(f):
        i = f * (n - 1)
        lo = int(i)
        hi = min(lo + 1, n - 1)
        return s[lo] + (i - lo) * (s[hi] - s[lo])
    return {'mean': sum(s) / n, 'min': s[0], 'max': s[-1],
            'p10': q(0.10), 'p50': q(0.50), 'p90': q(0.90)}


human = fasta_one(need(HUMAN_FASTA, 'the human sequence the numbering is '
                                    'checked against is not on disk'))
chicken = fasta_entity(need(CHICKEN_FASTA, 'the chicken sequence the numbering '
                                           'is checked against is not on disk'),
                       '7BC6_1')[:473]
if len(human) != 459:
    die('the cached human sequence is %d residues and this project works with 459'
        % len(human))
if len(chicken) != 473:
    die('the native part of the chicken entity is %d residues and this project '
        'works with 473' % len(chicken))

m = json.load(open(need(MAP, 'the chicken to human mapping is not on disk. Run '
                             'codes/01_residue_map.py first')))
CK2HU = {int(k): int(v) for k, v in m['ck2hu'].items()}

posset = tsv(need(POSSET, 'the pharmacophore position set is not on disk'))
bad = []
for r in posset:
    cp, hp = int(r['chicken_pos']), int(r['human_pos'])
    if CK2HU.get(cp) != hp:
        bad.append('chicken %d maps to %s and the set pairs it with human %d'
                   % (cp, CK2HU.get(cp), hp))
    elif chicken[cp - 1] != r['chicken_aa']:
        bad.append('chicken %d is %s in the entity and %s in the set'
                   % (cp, chicken[cp - 1], r['chicken_aa']))
    elif human[hp - 1] != r['human_aa']:
        bad.append('human %d is %s in the sequence and %s in the set'
                   % (hp, human[hp - 1], r['human_aa']))
if bad:
    die('the mapping disagrees with %s at %d of %d positions, the first few are:'
        '\n  %s' % (POSSET, len(bad), len(posset), '\n  '.join(bad[:5])))
PHARM_CK = sorted(int(r['chicken_pos']) for r in posset)
log('mapping agrees with %s at all %d positions of the set'
    % (os.path.basename(POSSET), len(posset)))

ac = json.load(open(need(CONSTRUCT, 'the construct is not on disk')))
RESCUER_CK = sorted(int(k) for k in ac['rescuer_positions'])
TARGET_CK = sorted(int(k) for k in ac['chosen'])
PRIMARY_HU = sorted(int(p) for p in ac['primaries'])
HU2CK = {v: k for k, v in CK2HU.items()}
PRIMARY_CK = sorted(HU2CK[p] for p in PRIMARY_HU if p in HU2CK)
if len(PRIMARY_CK) != len(PRIMARY_HU):
    die('%d of the %d primaries have no chicken position in the mapping'
        % (len(PRIMARY_HU) - len(PRIMARY_CK), len(PRIMARY_HU)))

freq_paths = sorted(glob.glob(os.path.join(EVDIR, '*', '*', 'align',
                                           '*_frequencies.csv')))
enr_paths = sorted(glob.glob(os.path.join(EVDIR, '*', '*', 'couplings',
                                          '*_enrichment.csv')))
stat_paths = sorted(glob.glob(os.path.join(EVDIR, '*', '*', 'align',
                                           '*_alignment_statistics.csv')))
for paths, what in ((freq_paths, 'column frequencies'),
                    (enr_paths, 'coupling enrichment'),
                    (stat_paths, 'alignment statistics')):
    if len(paths) != 1:
        die('%d files of %s under %s, and exactly one run is scored'
            % (len(paths), what, EVDIR))

stats = list(csv.DictReader(open(stat_paths[0])))
if len(stats) != 1:
    die('%s carries %d rows and one alignment is scored'
        % (stat_paths[0], len(stats)))
stats = stats[0]
MIN_COL_COV = float(stats['minimum_column_coverage'])
NUM_COV = int(stats['num_cov'])
N_SEQS = int(stats['num_seqs'])
SEQLEN = int(stats['seqlen'])

freq = list(csv.DictReader(open(freq_paths[0])))
if len(freq) != SEQLEN or SEQLEN != len(human):
    die('%s carries %d columns, the statistics say %d and the human chain is %d'
        % (freq_paths[0], len(freq), SEQLEN, len(human)))
off = [int(r['i']) for r in freq if r['A_i'] != human[int(r['i']) - 1]]
if off:
    die('%s names a different residue from the human sequence at %d columns, '
        'the first few being %s' % (freq_paths[0], len(off), off[:5]))

CONS, GAP = {}, {}
for r in freq:
    i = int(r['i'])
    CONS[i] = float(r['conservation'])
    GAP[i] = float(r['-'])
COVERED = {i for i, g in GAP.items() if g <= 1.0 - MIN_COL_COV}
if len(COVERED) != NUM_COV:
    die('a gap frequency of at most %.2f leaves %d columns and %s declares %d '
        'covered, so the covered region cannot be taken from the frequencies'
        % (1.0 - MIN_COL_COV, len(COVERED), os.path.basename(stat_paths[0]),
           NUM_COV))
log('covered columns %d, which reproduces the %d the alignment statistics '
    'declare at a column coverage of %.2f' % (len(COVERED), NUM_COV,
                                              MIN_COL_COV))

ENR = {}
for r in csv.DictReader(open(enr_paths[0])):
    i = int(r['i'])
    if r['A_i'] != human[i - 1]:
        die('%s names %s at human %d and the sequence has %s there'
            % (enr_paths[0], r['A_i'], i, human[i - 1]))
    ENR[i] = float(r['enrichment'])

runs = {}
WATCHED_CK = set()
for arm in ARMS:
    reports = sorted(glob.glob(os.path.join(MD, '%s_r*_report.json' % arm)))
    chains = sorted(glob.glob(os.path.join(MD, '%s_r*_rmsf_chain.xvg' % arm)))
    if not reports:
        die('no run report for the %s arm under %s' % (arm, MD))
    if len(chains) != len(reports):
        die('the %s arm has %d run reports and %d chain files'
            % (arm, len(reports), len(chains)))
    per = []
    for rp in reports:
        d = json.load(open(rp))
        cf = rp.replace('_report.json', '_rmsf_chain.xvg')
        if not os.path.isfile(cf):
            die('%s has no chain file beside it' % rp)
        rows = []
        for L in open(cf):
            L = L.strip()
            if not L or L[0] in '@#':
                continue
            a, b = L.split()
            rows.append((int(a), float(b) * NM_TO_A))
        if len(rows) != d['residues']:
            die('%s carries %d residues and %s reports %d'
                % (os.path.basename(cf), len(rows), os.path.basename(rp),
                   d['residues']))
        by = dict(rows)
        if len(by) != len(rows):
            die('%s repeats a residue number' % os.path.basename(cf))

        w = d['rmsf_A_whole_chain']
        vals = [v for _, v in rows]
        checks = [('first', w['first'], rows[0][1]),
                  ('last', w['last'], rows[-1][1]),
                  ('mean', w['mean'], sum(vals) / len(vals))]
        for pos, want in d['rmsf_A_by_chicken_position'].items():
            checks.append(('position %s' % pos, want, by.get(int(pos))))
        for pos, want in d['rmsf_A_at_watched'].items():
            checks.append(('watched %s' % pos, want, by.get(int(pos))))
        off = ['%s reports %s and the chain file gives %s'
               % (n, a, 'nothing' if b is None else '%.4f' % b)
               for n, a, b in checks if b is None or abs(a - b) > TOL_A]
        if off:
            die('%s disagrees with %s at %d of %d values, the first few being:'
                '\n  %s' % (os.path.basename(cf), os.path.basename(rp),
                            len(off), len(checks), '\n  '.join(off[:5])))
        WATCHED_CK.update(int(k) for k in d['rmsf_A_at_watched'])
        per.append({'replicate': int(d['replicate']), 'seed': str(d['seed']),
                    'ns': d['ns'], 'shift': int(d['numbering_shift']),
                    'checks': len(checks), 'by': by})
        log('%-28s %d residues, %d values agree with the run report'
            % (os.path.basename(cf), len(rows), len(checks)))
    if len({tuple(sorted(p['by'])) for p in per}) != 1:
        die('the replicates of the %s arm do not cover the same residues' % arm)
    if len({p['shift'] for p in per}) != 1:
        die('the replicates of the %s arm carry different numbering shifts' % arm)
    if len({p['seed'] for p in per}) != len(per):
        die('two replicates of the %s arm share a seed' % arm)
    runs[arm] = per

counts = {len(v) for v in runs.values()}
if len(counts) != 1:
    die('the arms carry %s replicates, and the normalisation has to have the '
        'same number under each' % ', '.join('%s=%d' % (a, len(runs[a]))
                                             for a in ARMS))
N_REP = counts.pop()
log('')
log('%d arms at %d replicates each, %g ns per run'
    % (len(ARMS), N_REP, runs[ARMS[0]][0]['ns']))

os.makedirs(OUTDIR, exist_ok=True)
summary = {'method': {
    'definition': 'DCS = conservation * norm_dynamics, RCS = conservation * '
                  '(1 - norm_dynamics), DCopS and RCopS the same two products '
                  'of the coupling enrichment, each rescaled 0 to 10 over the '
                  'scored positions',
    'dynamics': 'per-residue RMSF of the production run, averaged over the '
                'replicates of the arm',
    'normalisation': 'min-max over the whole simulated chain',
    'replicas': N_REP,
    'replicas_in_published_method': REPLICAS_IN_PUBLISHED_METHOD,
    'conservation': 'EVcouplings column conservation of the human alignment',
    'enrichment': 'EVcouplings coupling enrichment of the human alignment',
    'alignment_sequences': N_SEQS,
    'minimum_column_coverage': MIN_COL_COV,
    'covered_columns': NUM_COV,
    'scored_region_source': 'the columns the EVcouplings run covers, checked '
                            'against its own alignment statistics',
    'top_fraction': TOP_FRACTION,
    'published_cuts': {'dcs': PUBLISHED_DCS_CUT, 'rcs': PUBLISHED_RCS_CUT,
                       'applied': False,
                       'why': 'calibrated over a set of many proteins, and the '
                              'rescaling here runs over one'}},
    'arms': {}}

scored = {}
for arm in ARMS:
    per = runs[arm]
    pos = sorted(per[0]['by'])
    avg = {p: sum(r['by'][p] for r in per) / len(per) for p in pos}
    nd = dict(zip(pos, minmax([avg[p] for p in pos])))

    out = []
    for p in pos:
        hu = CK2HU.get(p)
        covered = hu in COVERED if hu else False
        out.append({'chicken': p, 'human': hu,
                    'chicken_aa': chicken[p - 1],
                    'human_aa': human[hu - 1] if hu else None,
                    'rmsf_A': avg[p],
                    'rmsf_A_by_replicate': [r['by'][p] for r in per],
                    'norm_dynamics': nd[p],
                    'conservation': CONS[hu] if covered else None,
                    'enrichment': ENR.get(hu) if covered else None,
                    'covered': covered,
                    'in_position_set': p in PHARM_CK,
                    'is_primary': p in PRIMARY_CK,
                    'is_rescuer': p in RESCUER_CK,
                    'is_rescue_target': p in TARGET_CK})

    usable = [r for r in out if r['conservation'] is not None]
    if len(usable) < 2:
        die('the %s arm has %d positions carrying both a fluctuation and a '
            'conservation value, which is not enough to score' % (arm, len(usable)))
    raw_d = [r['conservation'] * r['norm_dynamics'] for r in usable]
    raw_r = [r['conservation'] * (1.0 - r['norm_dynamics']) for r in usable]
    for r, a, b in zip(usable, minmax(raw_d), minmax(raw_r)):
        r['dcs'] = a * 10.0
        r['rcs'] = b * 10.0

    coupled = [r for r in usable if r['enrichment'] is not None]
    if len(coupled) >= 2:
        raw_dc = [r['enrichment'] * r['norm_dynamics'] for r in coupled]
        raw_rc = [r['enrichment'] * (1.0 - r['norm_dynamics']) for r in coupled]
        for r, a, b in zip(coupled, minmax(raw_dc), minmax(raw_rc)):
            r['dcops'] = a * 10.0
            r['rcops'] = b * 10.0

    for name, rows, hi, lo in (('conservation', usable, 'dcs', 'rcs'),
                               ('enrichment', coupled, 'dcops', 'rcops')):
        if not rows or hi not in rows[0]:
            continue
        x = [r['norm_dynamics'] for r in rows]
        ch, cl = pearson(x, [r[hi] for r in rows]), pearson(x, [r[lo] for r in rows])
        if not (ch > 0 and cl < 0):
            die('in the %s arm the dynamic %s score correlates with the '
                'fluctuation at %.3f and the rigid one at %.3f, and the '
                'definitions require the first positive and the second negative'
                % (arm, name, ch, cl))

    path = os.path.join(OUTDIR, '%s.tsv' % arm)
    with open(path, 'w') as fh:
        fh.write('chicken_pos\thuman_pos\tchicken_aa\thuman_aa\trmsf_A\t'
                 'norm_dynamics\tconservation\tenrichment\tdcs_0_10\trcs_0_10\t'
                 'dcops_0_10\trcops_0_10\tin_position_set\tis_primary\t'
                 'is_rescuer\tis_rescue_target\n')
        for r in out:
            def f(k, fmt='%.4f'):
                return '' if r.get(k) is None else fmt % r[k]
            fh.write('%d\t%s\t%s\t%s\t%.4f\t%.6f\t%s\t%s\t%s\t%s\t%s\t%s\t'
                     '%d\t%d\t%d\t%d\n'
                     % (r['chicken'], '' if r['human'] is None else r['human'],
                        r['chicken_aa'], r['human_aa'] or '', r['rmsf_A'],
                        r['norm_dynamics'], f('conservation'), f('enrichment'),
                        f('dcs', '%.3f'), f('rcs', '%.3f'),
                        f('dcops', '%.3f'), f('rcops', '%.3f'),
                        r['in_position_set'], r['is_primary'], r['is_rescuer'],
                        r['is_rescue_target']))

    k = max(1, int(round(TOP_FRACTION * len(usable))))
    top_d = sorted(usable, key=lambda r: -r['dcs'])[:k]
    top_r = sorted(usable, key=lambda r: -r['rcs'])[:k]
    top_d_set = {r['chicken'] for r in top_d}
    top_r_set = {r['chicken'] for r in top_r}

    def share(rows, which):
        if not rows:
            return None
        return sum(1 for r in rows if r['chicken'] in which) / float(len(rows))

    def brief(rows, key):
        return [{'chicken': r['chicken'], 'human': r['human'],
                 'chicken_aa': r['chicken_aa'], 'human_aa': r['human_aa'],
                 key: r[key], 'rmsf_A': r['rmsf_A'],
                 'conservation': r['conservation'],
                 'in_position_set': r['in_position_set'],
                 'is_rescuer': r['is_rescuer']} for r in rows]

    groups = {}
    for label, members in (('position_set', PHARM_CK), ('primaries', PRIMARY_CK),
                           ('rescuers', RESCUER_CK), ('rescue_targets', TARGET_CK)):
        rows = [r for r in usable if r['chicken'] in members]
        groups[label] = {
            'positions': len(rows),
            'not_scored': sorted(p for p in members
                                 if p not in {r['chicken'] for r in usable}),
            'dcs_mean': (sum(r['dcs'] for r in rows) / len(rows)) if rows else None,
            'rcs_mean': (sum(r['rcs'] for r in rows) / len(rows)) if rows else None,
            'rmsf_A_mean': (sum(r['rmsf_A'] for r in rows) / len(rows)) if rows else None,
            'conservation_mean': (sum(r['conservation'] for r in rows) / len(rows)) if rows else None,
            'share_in_top_dcs': share(rows, top_d_set),
            'share_in_top_rcs': share(rows, top_r_set)}

    rank_d = {r['chicken']: i + 1 for i, r in
              enumerate(sorted(usable, key=lambda x: -x['dcs']))}
    rank_r = {r['chicken']: i + 1 for i, r in
              enumerate(sorted(usable, key=lambda x: -x['rcs']))}
    watch = {}
    for p in sorted(set(TARGET_CK) | WATCHED_CK):
        hit = [r for r in usable if r['chicken'] == p]
        if not hit:
            continue
        r = hit[0]
        watch[str(p)] = {'rmsf_A': r['rmsf_A'], 'norm_dynamics': r['norm_dynamics'],
                         'dcs': r['dcs'], 'rcs': r['rcs'],
                         'conservation': r['conservation'],
                         'rank_dcs': rank_d[p], 'rank_rcs': rank_r[p],
                         'in_top_dcs': p in top_d_set, 'in_top_rcs': p in top_r_set}

    summary['arms'][arm] = {
        'arm': arm, 'replicates': len(per),
        'seeds': [r['seed'] for r in per], 'ns_per_run': per[0]['ns'],
        'residues': len(out), 'scored': len(usable), 'with_enrichment': len(coupled),
        'values_checked_against_run_reports': sum(r['checks'] for r in per),
        'rmsf_A': quantiles([r['rmsf_A'] for r in out]),
        'dcs_0_10': quantiles([r['dcs'] for r in usable]),
        'rcs_0_10': quantiles([r['rcs'] for r in usable]),
        'top_n': k,
        'dynamic_conserved': brief(top_d, 'dcs'),
        'rigid_conserved': brief(top_r, 'rcs'),
        'above_published_dcs_cut': sum(1 for r in usable if r['dcs'] > PUBLISHED_DCS_CUT),
        'above_published_rcs_cut': sum(1 for r in usable if r['rcs'] > PUBLISHED_RCS_CUT),
        'share_of_chain_in_top_dcs': k / float(len(usable)),
        'groups': groups,
        'watched': watch}
    scored[arm] = {r['chicken']: r for r in usable}
    log('%-10s %3d residues, %3d scored, %3d with enrichment, RMSF %.2f A mean, '
        'top %d by each score' % (arm, len(out), len(usable), len(coupled),
                                  summary['arms'][arm]['rmsf_A']['mean'], k))

summary['between_arms'] = {
    'what': 'each arm is rescaled over itself, so the comparison is on the '
            'ranking within an arm and on the fluctuation, not on a difference '
            'of rescaled scores'}
base = scored[ARMS[0]]
for arm in ARMS[1:]:
    other = scored[arm]
    both = sorted(set(base) & set(other))
    if len(both) < 2:
        die('the %s arm and the %s arm share %d scored positions, which is not '
            'enough to compare' % (ARMS[0], arm, len(both)))

    def ranks(d):
        order = sorted(both, key=lambda p: d[p]['dcs'])
        return {p: i + 1 for i, p in enumerate(order)}
    ra, rb = ranks(base), ranks(other)
    d_rmsf = [other[p]['rmsf_A'] - base[p]['rmsf_A'] for p in both]
    moved = sorted(both, key=lambda p: -abs(rb[p] - ra[p]))[:10]
    summary['between_arms']['%s_against_%s' % (arm, ARMS[0])] = {
        'positions': len(both),
        'dcs_rank_correlation': pearson([ra[p] for p in both], [rb[p] for p in both]),
        'rmsf_correlation': pearson([base[p]['rmsf_A'] for p in both],
                                    [other[p]['rmsf_A'] for p in both]),
        'delta_rmsf_A_mean': sum(d_rmsf) / len(d_rmsf),
        'delta_rmsf_A_mean_abs': sum(abs(v) for v in d_rmsf) / len(d_rmsf),
        'rank_shift_mean': sum(abs(rb[p] - ra[p]) for p in both) / float(len(both)),
        'moved_most_in_dcs_rank': [
            {'chicken': p, 'human': base[p]['human'],
             'chicken_aa': base[p]['chicken_aa'],
             'rank_in_%s' % ARMS[0]: ra[p], 'rank_in_%s' % arm: rb[p],
             'rmsf_A_%s' % ARMS[0]: base[p]['rmsf_A'],
             'rmsf_A_%s' % arm: other[p]['rmsf_A'],
             'is_rescuer': base[p]['is_rescuer'],
             'in_position_set': base[p]['in_position_set']} for p in moved]}

USABLE_CK = sorted(scored[ARMS[0]])


def score_from(by):
    pos = sorted(by)
    nd = dict(zip(pos, minmax([by[p] for p in pos])))
    use = [p for p in USABLE_CK if p in nd]
    cons = [CONS[CK2HU[p]] for p in use]
    d = minmax([c * nd[p] for c, p in zip(cons, use)])
    r = minmax([c * (1.0 - nd[p]) for c, p in zip(cons, use)])
    return ({p: v * 10.0 for p, v in zip(use, d)},
            {p: v * 10.0 for p, v in zip(use, r)})


single = {}
for arm in ARMS:
    for rep in runs[arm]:
        single[(arm, rep['replicate'])] = score_from(rep['by'])
REPS = sorted({k[1] for k in single})
direction = {}
for p in sorted(set(TARGET_CK) | WATCHED_CK):
    if p not in USABLE_CK:
        continue
    row = {}
    for key, which in (('dcs', 0), ('rcs', 1)):
        pooled = (scored[ARMS[-1]][p][key] - scored[ARMS[0]][p][key])
        each = [single[(ARMS[-1], r)][which][p] - single[(ARMS[0], r)][which][p]
                for r in REPS]
        row[key] = {
            'pooled_change': pooled,
            'change_by_replicate': each,
            'every_replicate_agrees': all((v > 0) == (pooled > 0) for v in each)}
    direction[str(p)] = row
MD_CK = set(runs[ARMS[0]][0]['by'])
summary['not_scored'] = {}
for label, members in (('position_set', PHARM_CK), ('rescue_targets', TARGET_CK),
                       ('rescuers', RESCUER_CK), ('primaries', PRIMARY_CK)):
    rows = []
    for p in members:
        if p in USABLE_CK:
            continue
        hu = CK2HU.get(p)
        if p not in MD_CK:
            why = 'outside the simulated span'
        elif hu is None:
            why = 'no human counterpart in the alignment'
        else:
            why = 'human column gapped in %.1f per cent of the alignment' % (
                100.0 * GAP[hu])
        rows.append({'chicken': p, 'human': hu, 'simulated': p in MD_CK,
                     'gap_frequency': None if hu is None else GAP[hu],
                     'why': why})
    summary['not_scored'][label] = rows
    for r in rows:
        log('chicken %-4d of the %s carries no score, %s'
            % (r['chicken'], label.replace('_', ' '), r['why']))

summary['single_replicate_direction'] = {
    'what': 'the change in each score from the %s arm to the %s arm, recomputed '
            'with one replicate under the normalisation, beside the change the '
            'pooled score gives' % (ARMS[0], ARMS[-1]),
    'replicates': REPS,
    'positions': direction}
for p, row in sorted(direction.items(), key=lambda kv: int(kv[0])):
    log('chicken %-4s dynamic %+.2f pooled, %s by replicate, every replicate '
        'agrees %s' % (p, row['dcs']['pooled_change'],
                       ', '.join('%+.2f' % v for v in row['dcs']['change_by_replicate']),
                       row['dcs']['every_replicate_agrees']))

json.dump(summary, open(OUT, 'w'), indent=1, sort_keys=True)
log('')
for arm in ARMS:
    s = summary['arms'][arm]
    g = s['groups']['position_set']
    log('%-10s position set %.2f of the top dynamic-conserved and %.2f of the '
        'top rigid-conserved, against %.2f of the chain in each'
        % (arm, g['share_in_top_dcs'], g['share_in_top_rcs'],
           s['share_of_chain_in_top_dcs']))
for key, c in sorted(summary['between_arms'].items()):
    if not isinstance(c, dict) or 'positions' not in c:
        continue
    log('%-26s %d shared positions, dynamic ranking correlates at %.3f, '
        'fluctuation at %.3f, mean change %+.3f A'
        % (key, c['positions'], c['dcs_rank_correlation'],
           c['rmsf_correlation'], c['delta_rmsf_A_mean']))
log('')
log('wrote %s and %d per arm tables in %s' % (OUT, len(scored), OUTDIR))
with open(LOG, 'w') as f:
    f.write('\n'.join(lines) + '\n')
