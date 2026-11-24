import glob
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, 'data')
MD = os.path.join(ROOT, 'md')


def die(msg):
    sys.exit('STOP: %s' % msg)


def load(name, where=D):
    p = os.path.join(where, name)
    if not os.path.isfile(p):
        die('missing %s' % p)
    with open(p) as fh:
        return json.load(fh)


def need(d, *path):
    cur, seen = d, []
    for k in path:
        seen.append(str(k))
        try:
            cur = cur[k]
        except (KeyError, IndexError, TypeError):
            die('no value at %s' % ' -> '.join(seen))
    return cur


def logtext(name):
    p = os.path.join(D, name)
    if not os.path.isfile(p):
        die('missing log %s' % p)
    return open(p, errors='replace').read()


def grab(text, pattern, what, cast=float):
    m = re.search(pattern, text)
    if not m:
        die('%s not found by %r' % (what, pattern))
    return cast(m.group(1))


def tsv(name):
    p = os.path.join(D, name)
    if not os.path.isfile(p):
        die('missing %s' % p)
    rows = [l.rstrip('\n').split('\t') for l in open(p)]
    head = rows[0]
    return [dict(zip(head, r)) for r in rows[1:] if any(r)]


N = {}

L1 = logtext('01_residue_map_log.txt')
N['map'] = {
    'human_len': grab(L1, r'human\s+Q96NT5\s*:\s*(\d+) aa', 'human length', int),
    'chicken_len': grab(L1, r'chicken 7BC6 chain A\s*:\s*(\d+) aa', 'chicken length', int),
    'aligned_columns': grab(L1, r'aligned columns (\d+)', 'aligned columns', int),
    'identical': grab(L1, r'identical (\d+)', 'identical count', int),
    'identity_pct': grab(L1, r'identical \d+ \(([\d.]+)%', 'identity per cent'),
    'optimal_paths': grab(L1, r'alignment score [\d.]+, (\d+) optimal path', 'optimal paths', int),
    'stated_pairs': grab(L1, r'(\d+)/\d+ stated pairs reproduced', 'stated pairs', int),
    'stated_pairs_total': grab(L1, r'\d+/(\d+) stated pairs reproduced', 'stated pairs total', int),
}
L2 = logtext('02_map_parameter_sweep_log.txt')
N['map']['settings'] = grab(L2, r'settings swept: (\d+)', 'parameter settings', int)
N['map']['cooptimal_paths'] = grab(L2, r'alignment paths examined: (\d+)',
                                   'co-optimal paths', int)
N['map']['probed_positions'] = grab(L2, r'unambiguous positions: (\d+)/',
                                    'probed positions', int)
N['map']['positions_changing'] = grab(
    L2, r'positions whose human partner changes with the parameters: (\d+)',
    'positions changing', int)

MAP = json.load(open(os.path.join(ROOT, 'structures', '7BC6_to_human.json')))
pairs = sorted((int(k), v) for k, v in MAP['ck2hu'].items() if v is not None)
offs = [(c, c - h) for c, h in pairs]
runs, prev = [], None
for c, o in offs:
    if o != prev:
        runs.append([o, c, c])
        prev = o
    else:
        runs[-1][2] = c
N['map']['offset'] = {
    'values_chain': sorted({o for _, o in offs}),
    'n_values_chain': len({o for _, o in offs}),
    'change_positions': [r[1] for r in runs[1:]],
    'mapped_first': offs[0][0],
    'mapped_last': offs[-1][0],
}
soff = [(int(r['chicken_pos']), int(r['chicken_pos']) - int(r['human_pos']))
        for r in tsv('pharmacophore_position_set.tsv')]
sruns, prev = [], None
for c, o in sorted(soff):
    if o != prev:
        sruns.append([o, c, c])
        prev = o
    else:
        sruns[-1][2] = c
counts = {}
for _, o in soff:
    counts[o] = counts.get(o, 0) + 1
modal = max(counts, key=counts.get)
N['map']['offset'].update({
    'n_values_set': len(counts),
    'set_runs': sruns,
    'modal_offset': modal,
    'modal_count': counts[modal],
})

pset = tsv('pharmacophore_position_set.tsv')
N['set'] = {
    'n_positions': len(pset),
    'n_primaries': sum(1 for r in pset if r['differs'] == '1'),
    'n_parker_named': sum(1 for r in pset if r['in_parker_named'] == '1'),
    'n_mutation_literature': sum(1 for r in pset if r['in_mutation_literature'] == '1'),
    'n_pemetrexed_contact': sum(1 for r in pset
                                if r['min_pemetrexed_dist_A'] not in ('', None)),
}
hfm = tsv('hfm_parker_table.tsv')
N['set']['n_hfm_rows'] = len(hfm)
_d = sorted(float(r['min_pemetrexed_dist_A']) for r in pset
            if r['min_pemetrexed_dist_A'] not in ('', None))
if not _d:
    die('no position of the set carries a distance to the bound antifolate')
N['set']['pemetrexed_dist_min'] = _d[0]
N['set']['pemetrexed_dist_max'] = _d[-1]

pc = load('primary_coupling.json')
N['coupled'] = {
    'n_primary': need(pc, 'n_primary'),
    'n_pairs': need(pc, 'n_pairs'),
    'window': need(pc, 'window'),
    'cutoffs': need(pc, 'cutoffs'),
}
L6 = logtext('06_primary_coupling_log.txt')
ladder = {}
for line in L6.splitlines():
    m = re.match(r'\s*([\d.]+) A\s+pairs\s+(\d+)\s+outside \+/-(\d+)\s+(\d+)'
                 r'\s+max sep\s+(\d+)\s+coupled positions\s+(\d+)/(\d+)', line)
    if m:
        cut, pairs, win, out, sep, cp, tot = m.groups()
        ladder[float(cut)] = dict(pairs=int(pairs), window=int(win),
                                  beyond_window=int(out), max_separation=int(sep),
                                  coupled=int(cp), of=int(tot))
if not ladder:
    die('no cutoff ladder rows in the primary coupling log')
N['coupled']['ladder'] = {str(k): v for k, v in sorted(ladder.items())}
REPORTED_CUT = 6.0
if REPORTED_CUT not in ladder:
    die('the %.1f Angstrom cutoff reported is not in the ladder, '
        'which has %s' % (REPORTED_CUT, sorted(ladder)))
N['coupled']['reported_cutoff_A'] = REPORTED_CUT
N['coupled'].update({'n_contacting': ladder[REPORTED_CUT]['coupled'],
                     'n_coupled_pairs': ladder[REPORTED_CUT]['pairs'],
                     'n_pairs_beyond_window': ladder[REPORTED_CUT]['beyond_window'],
                     'largest_separation': ladder[REPORTED_CUT]['max_separation']})

ev = load('evcouplings_report.json')
dep = load('evcouplings_depth.json')
dca = load('dca_diagnostics.json')
N['dca'] = {
    'n_sequences': need(ev, 'num_seqs'),
    'meff': need(ev, 'meff'),
    'meff_per_length': need(ev, 'meff_per_length'),
    'fitted_columns': need(ev, 'fitted_columns'),
    'length': need(ev, 'length'),
    'own_alignment_meff_per_length': need(dca, 'meff_per_length'),
    'own_alignment_n_sequences': need(dca, 'n_sequences'),
    'gap_fraction': need(dep, 'gap_fraction'),
    'theta': need(dep, 'theta'),
}
L11 = logtext('11_evcouplings_read_log.txt')
N['dca']['contact_cutoff_A'] = grab(
    L11, r'precision of the strongest non-local couplings against the structure, '
         r'contact at ([\d.]+) Angstrom', 'contact cutoff')
N['dca']['precision_L2'] = grab(L11, r'\s+L/2\s+\d+\s+([\d.]+)',
                                'precision at L/2')
N['dca']['n_at_L2'] = grab(L11, r'\s+L/2\s+(\d+)\s+[\d.]+', 'N at L/2', int)
N['dca']['median_distance_L2'] = grab(
    L11, r'\s+L/2\s+\d+\s+[\d.]+\s+([\d.]+)', 'median distance at L/2')
N['dca']['n_nonlocal_pairs'] = grab(
    L11, r'all (\d+) non-local pairs', 'non-local pair count', int)
N['dca']['baseline_contact_rate'] = grab(
    L11, r'all \d+ non-local pairs: ([\d.]+) are in contact',
    'baseline contact rate')
N['dca']['n_nonlocal_with_distance'] = grab(
    L11, r'non-local couplings with a distance in the structure: (\d+)',
    'non-local couplings with a distance', int)

rs = load('rescue_sweep_report.json')
pub = [r for r in need(rs, 'summary') if r['policy'] == 'published']
if not pub:
    die('no published-policy rows in rescue_sweep_report.json')
N['sweep'] = {
    'threshold': need(rs, 'threshold'),
    'radii': [r['radius'] for r in pub],
    'crossings_min': min(r['rescuers_crossing'] for r in pub),
    'crossings_max': max(r['rescuers_crossing'] for r in pub),
    'targets_rescued': sorted({r['targets_rescued'] for r in pub}),
    'rescued_positions': sorted({p for r in pub for p in r['rescued_positions']}),
    'nearest_separations': sorted({r['nearest_rescuer_separation'] for r in pub}),
    'n_rows': len(need(rs, 'summary')),
    'protected_count': need(rs, 'protected_count'),
}

N['sweep']['by_policy'] = {}
for pol in sorted({r['policy'] for r in need(rs, 'summary')}):
    rows = [r for r in need(rs, 'summary') if r['policy'] == pol]
    N['sweep']['by_policy'][pol] = {
        'crossings_min': min(r['rescuers_crossing'] for r in rows),
        'crossings_max': max(r['rescuers_crossing'] for r in rows),
        'targets_rescued': sorted({r['targets_rescued'] for r in rows}),
        'rescued_positions': sorted({p for r in rows for p in r['rescued_positions']}),
    }

pr = load('parker_residues.json')
N['parker'] = {
    'source': need(pr, 'source'),
    'n': len(need(pr, 'residues')),
    'residues': need(pr, 'residues'),
}

cg = load('coupling_guided_rescue.json')
N['coupling_guided'] = {
    'window_radius': need(cg, 'window_radius'),
    'k_ladder': need(cg, 'k_ladder'),
    'contact_cutoff_A': need(cg, 'contact_cutoff_A'),
    'n_targets': len(need(cg, 'targets')),
}
L15 = logtext('15_coupling_guided_rescue_log.txt')
N['coupling_guided']['rankable_positions'] = grab(
    L15, r'positions rankable by coupling (\d+)', 'rankable positions', int)
N['coupling_guided']['unrankable_positions'] = grab(
    L15, r'positions rankable by coupling \d+, unrankable (\d+)',
    'unrankable positions', int)
N['coupling_guided']['matched_k'] = grab(
    L15, r'the published window of \d+ admits (\d+) positions', 'matched K', int)

dr = load('double_rescue_report.json')
N['pairs_407'] = {
    'target': need(dr, 'target_chicken_pos'),
    'human_pos': need(dr, 'target_human_pos'),
    'n_doubles': need(dr, 'n_doubles'),
    'n_singles': need(dr, 'n_singles'),
    'n_candidate_positions': len(need(dr, 'candidate_positions')),
    'baseline': need(dr, 'p_baseline'),
    'worst_relative_difference_vs_step2': need(dr, 'worst_relative_difference_vs_step2'),
    'n_crossing': need(dr, 'n_doubles_crossing'),
    'best_single_p': need(dr, 'best_single', 'p'),
    'best_double': need(dr, 'best_double'),
}
cr = load('construct_report.json')
pa = need(cr, 'part_a')
N['pairs_144'] = {
    'target': need(pa, 'target'),
    'n_doubles': need(pa, 'n_doubles'),
    'n_crossing': need(pa, 'n_crossing'),
    'baseline': need(pa, 'baseline'),
    'best_double': need(pa, 'best_double'),
    'worst_relative_difference_vs_step2': need(pa, 'worst_relative_difference_vs_step2'),
}

surv = need(cr, 'survival')
N['set_effect'] = {
    'threshold': need(cr, 'threshold'),
    'targets': need(surv, 'targets'),
    'n_holding_final': sum(1 for t in need(surv, 'targets') if t['holds_final']),
    'n_targets': len(need(surv, 'targets')),
    'overstatement_min': min(t['overstatement_final'] for t in need(surv, 'targets')
                             if t['overstatement_final'] > 1),
    'overstatement_max': max(t['overstatement_final'] for t in need(surv, 'targets')),
}

crr = load('construct_rescue_report.json')
N['construct_background'] = {
    'n_scan_rows': need(crr, 'n_scan_rows'),
    'n_targets': len(need(crr, 'targets')),
    'background_effect': need(crr, 'background_effect'),
    'best_policy': need(crr, 'best_policy'),
}

ac = load('admissible_construct.json')
N['construct'] = {
    'threshold': need(ac, 'threshold'),
    'n_off_limits_strict': len(need(ac, 'off_limits_strict')),
    'n_primaries': len(need(ac, 'primaries')),
    'n_rescuer_positions': len(need(ac, 'rescuer_positions')),
    'n_positions_changed': need(ac, 'n_positions_changed'),
    'rescuer_positions': need(ac, 'rescuer_positions'),
    'unrescued': need(ac, 'unrescued'),
    'chosen': need(ac, 'chosen'),
    'choice_rule': need(ac, 'choice_rule'),
    'overruled_by_background': need(ac, 'overruled_by_background'),
    'sequence_length': len(need(ac, 'sequence')),
}

s7 = load('step7_report.json')
pb = need(s7, 'part_b')
N['pairs_322'] = {
    'target': need(pb, 'target'),
    'n_doubles': need(pb, 'n_doubles'),
    'n_singles': need(pb, 'n_singles'),
    'n_candidates': len(need(pb, 'candidates')),
    'baseline': need(pb, 'baseline'),
    'n_crossing': need(pb, 'n_crossing'),
    'n_crossing_sparing_the_set': need(pb, 'n_crossing_sparing_the_set'),
}

s8 = load('step8_report.json')
N['step8'] = {
    'target': need(s8, 'target'),
    'human_pos': need(s8, 'human_pos'),
    'n_pairs': need(s8, 'n_pairs'),
    'baseline': need(s8, 'baseline'),
    'n_candidates': len(need(s8, 'candidates')),
    'n_candidates_in_set': len(need(s8, 'candidates_in_set')),
    'threshold': need(s8, 'threshold'),
    'background': need(s8, 'background'),
}

_sing = [r for r in tsv(os.path.join('esm2_step8', 'c407_singles.tsv'))
         if r['p_wt_step2'] not in ('', 'nan', None)
         and float(r['p_wt_step2']) > 0]
if len(_sing) < 10:
    die('only %d of the candidate singles at the rescue site have a wild-type '
        'counterpart, which is too few to quote a median' % len(_sing))
_f = sorted(float(r['p_construct']) / float(r['p_wt_step2']) for r in _sing)


def _q(v, q):
    i = (len(v) - 1) * q
    lo, hi = int(i), min(int(i) + 1, len(v) - 1)
    return v[lo] + (v[hi] - v[lo]) * (i - lo)


N['step8']['singles_vs_wt'] = {
    'n': len(_f), 'n_scored': len(tsv(os.path.join('esm2_step8',
                                                   'c407_singles.tsv'))),
    'median': _q(_f, 0.5), 'p05': _q(_f, 0.05), 'p95': _q(_f, 0.95),
}

s9 = load('step9_report.json')
N['step9'] = {
    'threshold': need(s9, 'threshold'),
    'checks': need(s9, 'checks'),
    'n_profile': len(need(s9, 'profile')),
    'below_threshold': need(s9, 'below_threshold'),
    'below_threshold_step7': need(s9, 'below_threshold_step7'),
    'fold_over_step7': need(s9, 'fold_over_step7'),
    'rescues_holding': need(s9, 'rescues_holding'),
    'rescues_lost': need(s9, 'rescues_lost'),
    'rescuer_positions_swapped': need(s9, 'rescuer_positions_swapped'),
    'displaced_rescues': need(s9, 'displaced_rescues'),
    'neighbourhood_of_the_swap': need(s9, 'neighbourhood_of_the_swap'),
}

af3 = load('af3_models_report.json')
N['af3'] = {
    'seed': need(af3, 'seed'),
    'core_positions': need(af3, 'core_positions'),
    'references': need(af3, 'references'),
    'reference_gap': need(af3, 'reference_gap'),
    'arms': need(af3, 'arms'),
    'arm_pairs': need(af3, 'arm_pairs'),
    'furthest_arm_pair': need(af3, 'furthest_arm_pair'),
    'furthest_arm_to_reference': need(af3, 'furthest_arm_to_reference'),
    'arms_are_one_set': need(af3, 'arms_are_one_set'),
    'plddt_deltas_vs_baseline': need(af3, 'plddt_deltas_vs_baseline'),
    'watched_positions': need(af3, 'watched_positions'),
    'models_per_arm': {a: len(v) for a, v in need(af3, 'all_models').items()},
}
watch = need(af3, 'watched_positions')
N['af3']['plddt_range_at_watch'] = {
    a: {str(q): [min(m['watched'][str(q)] for m in models),
                 max(m['watched'][str(q)] for m in models)]
        for q in watch}
    for a, models in need(af3, 'all_models').items()
}

N['af3']['per_arm'] = {
    a: {k: v[k] for k in ('plddt_chain', 'plddt_set', 'plddt_substituted',
                          'plddt_rescuers', 'plddt_outside_core', 'ptm',
                          'ranking_score', 'fraction_disordered', 'has_clash',
                          'n_outside_core')}
    for a, v in need(af3, 'arms').items()
}
N['af3']['rmsd_to_reference'] = {
    a: {e: need(af3, 'arms', a, 'rmsd_to', e, 'rmsd') for e in need(af3, 'references')}
    for a in need(af3, 'arms')
}

plan = load('md_plan.json', MD)
spec = load(os.path.join('inputs', 'build_spec.json'), MD)
N['md'] = {
    'budget': need(plan, 'run_budget'),
    'build': need(plan, 'build'),
    'arms': need(spec, 'builds'),
    'formal_charge': need(spec, 'formal_charge'),
    'disulfides': need(spec, 'disulfides_expected'),
    'residue_range': need(spec, 'residue_range'),
    'n_residues': need(spec, 'n_residues'),
    'entry': need(spec, 'entry'),
    'state': need(spec, 'state'),
    'substitutions': {a: len(v) for a, v in need(spec, 'substitutions').items()},
}

MDR = os.path.join(D, 'md')
runs = {}
for a in sorted(need(spec, 'builds')):
    found = sorted(glob.glob(os.path.join(MDR, '%s_r*_report.json' % a)))
    if not found:
        raise SystemExit('STOP: no report under %s for the %s arm, so that arm has '
                         'no simulation values to report' % (MDR, a))
    runs[a] = [json.load(open(f)) for f in found]
    runs[a].sort(key=lambda r: r['replicate'])
ARMO = ['wt', 'primaries', 'construct']
if sorted(runs) != sorted(ARMO):
    raise SystemExit('STOP: the reports cover %s and the arms are %s'
                     % (sorted(runs), sorted(ARMO)))
nrep = {a: len(v) for a, v in runs.items()}
if len(set(nrep.values())) != 1:
    raise SystemExit('STOP: the arms have unequal replicate counts, %s, so a '
                     'comparison across them is not one' % nrep)
lengths = {r['ns'] for v in runs.values() for r in v}
if len(lengths) != 1:
    raise SystemExit('STOP: the replicates reached different lengths, %s, so a '
                     'comparison across them is not one' % sorted(lengths))
shifts = {r['numbering_shift'] for v in runs.values() for r in v}
if len(shifts) != 1:
    raise SystemExit('STOP: the replicates resolved different numbering shifts, %s, '
                     'so their per-position values are not on one numbering' % shifts)
for a, v in runs.items():
    seeds = [r['seed'] for r in v]
    if len(set(seeds)) != len(seeds):
        raise SystemExit('STOP: the %s arm has replicates sharing a seed, %s. They '
                         'are the same trajectory twice.' % (a, seeds))
WATCH_R = sorted(int(k) for k in runs['wt'][0]['rmsf_A_at_watched'])
PAIRS_R = sorted(runs['wt'][0]['pair_distances_A'])


def _dig(d, path):
    for k in path:
        d = d[k]
    return float(d)


def across(arm, *path):
    return [_dig(r, path) for r in runs[arm]]


def band(*path):
    out = {}
    for a in runs:
        v = across(a, *path)
        out[a] = {'mean': sum(v) / len(v), 'lo': min(v), 'hi': max(v),
                  'spread': max(v) - min(v), 'values': v}
    return out


def order_of(*path):
    means = {a: sum(across(a, *path)) / len(across(a, *path)) for a in runs}
    rank = sorted(means, key=lambda a: -means[a])
    holds = all([across(a, *path)[i] for a in rank] ==
                sorted([across(a, *path)[i] for a in rank], reverse=True)
                for i in range(nrep['wt']))
    return rank, holds


RMSF_W = {str(w): band('rmsf_A_at_watched', str(w)) for w in WATCH_R}
ORDER_W = {str(w): order_of('rmsf_A_at_watched', str(w)) for w in WATCH_R}
N['md']['run'] = {
    'ns': sorted(lengths)[0],
    'replicates': sorted(nrep.values())[0],
    'arms': ARMO,
    'seeds': {a: [r['seed'] for r in runs[a]] for a in runs},
    'frames': {a: [r['frames'] for r in runs[a]] for a in runs},
    'offload': sorted({r.get('offload_path') for v in runs.values() for r in v
                       if r.get('offload_path')}),
    'numbering_shift': sorted(shifts)[0],
    'watched': WATCH_R,
    'rmsd_A': band('rmsd_A', 'mean'),
    'rg_A': band('radius_of_gyration_A', 'mean'),
    'rmsf_chain_A': band('rmsf_A_whole_chain', 'mean'),
    'area_per_lipid_A2': band('box', 'area_per_lipid_A2', 'mean'),
    'temperature_K': band('box', 'temperature_K', 'mean'),
    'rmsf_watched_A': RMSF_W,
    'rmsf_watched_order': {w: ORDER_W[w][0] for w in ORDER_W},
    'rmsf_watched_order_holds': {w: ORDER_W[w][1] for w in ORDER_W},
    'pair_A': {nm: band('pair_distances_A', nm, 'mean') for nm in PAIRS_R},
    'position_set_size': runs['wt'][0]['position_set_size'],
    'energy_columns': runs['wt'][0]['box'].get('columns_read'),
    'rmsd_spread_max': max(band('rmsd_A', 'mean')[a]['spread'] for a in runs),
    'total_ns': sorted(lengths)[0] * sorted(nrep.values())[0] * len(ARMO),
}

ENVP = os.path.join(MDR, 'environment.json')
if not os.path.isfile(ENVP):
    raise SystemExit('STOP: %s is absent, so the compute environment cannot '
                     'be stated' % ENVP)
ENV = json.load(open(ENVP))
if ENV.get('n_runs') != N['md']['run']['replicates'] * len(ARMO):
    raise SystemExit('STOP: the environment record covers %s runs and the '
                     'reports cover %s, so the two do not describe the same '
                     'set' % (ENV.get('n_runs'),
                              N['md']['run']['replicates'] * len(ARMO)))
if not ENV.get('uniform'):
    raise SystemExit('STOP: the runs did not share one card, one GROMACS '
                     'build and one thread count, so one sentence cannot '
                     'state the environment: %s' % ENV.get('per_run'))
N['md']['environment'] = ENV

DAP = os.path.join(D, 'dynamics_aware_summary.json')
if not os.path.isfile(DAP):
    raise SystemExit('STOP: %s is absent, so the dynamics-aware scores cannot '
                     'be quoted. Run codes/35_dynamics_aware.py first' % DAP)
DA = json.load(open(DAP))
if sorted(DA['arms']) != sorted(ARMO):
    raise SystemExit('STOP: the dynamics-aware summary covers %s and the '
                     'simulation reports cover %s'
                     % (', '.join(sorted(DA['arms'])), ', '.join(sorted(ARMO))))
if DA['method']['replicas'] != N['md']['run']['replicates']:
    raise SystemExit('STOP: the dynamics-aware summary averaged %d replicates '
                     'and the simulation reports cover %d'
                     % (DA['method']['replicas'], N['md']['run']['replicates']))
DW = sorted(set(DA['arms'][ARMO[0]]['watched'])
            & set(DA['arms'][ARMO[-1]]['watched']), key=int)
if not DW:
    raise SystemExit('STOP: no position carries a dynamics-aware score in both '
                     'the %s and the %s arm' % (ARMO[0], ARMO[-1]))
DBET = DA['between_arms']['%s_against_%s' % (ARMO[-1], ARMO[0])]
DDIR = DA['single_replicate_direction']['positions']
N['dynamics'] = {
    'replicas': DA['method']['replicas'],
    'replicas_in_published_method': DA['method']['replicas_in_published_method'],
    'top_fraction': DA['method']['top_fraction'],
    'alignment_sequences': DA['method']['alignment_sequences'],
    'minimum_column_coverage': DA['method']['minimum_column_coverage'],
    'covered_columns': DA['method']['covered_columns'],
    'published_dcs_cut': DA['method']['published_cuts']['dcs'],
    'published_rcs_cut': DA['method']['published_cuts']['rcs'],
    'residues': DA['arms'][ARMO[0]]['residues'],
    'scored': DA['arms'][ARMO[0]]['scored'],
    'with_enrichment': DA['arms'][ARMO[0]]['with_enrichment'],
    'top_n': DA['arms'][ARMO[0]]['top_n'],
    'watched': DW,
    'set_positions_scored':
        DA['arms'][ARMO[0]]['groups']['position_set']['positions'],
    'set_conservation_mean':
        DA['arms'][ARMO[0]]['groups']['position_set']['conservation_mean'],
    'not_scored': {k: [{'chicken': r['chicken'], 'human': r['human'],
                        'gap_percent': round(100.0 * r['gap_frequency'], 1),
                        'why': r['why']}
                       for r in v if r.get('gap_frequency') is not None]
                   for k, v in DA['not_scored'].items()},
    'rank_correlation': DBET['dcs_rank_correlation'],
    'rmsf_correlation': DBET['rmsf_correlation'],
    'arm': {a: {
        'above_published_dcs_cut': DA['arms'][a]['above_published_dcs_cut'],
        'above_published_rcs_cut': DA['arms'][a]['above_published_rcs_cut'],
        'chain_share_in_top': DA['arms'][a]['share_of_chain_in_top_dcs'],
        'set_share_top_dcs':
            DA['arms'][a]['groups']['position_set']['share_in_top_dcs'],
        'set_share_top_rcs':
            DA['arms'][a]['groups']['position_set']['share_in_top_rcs'],
        'rescuers': {k: DA['arms'][a]['groups']['rescuers'][k]
                     for k in ('positions', 'dcs_mean', 'rcs_mean',
                               'rmsf_A_mean', 'share_in_top_dcs',
                               'share_in_top_rcs')},
        'targets': {k: DA['arms'][a]['groups']['rescue_targets'][k]
                    for k in ('positions', 'dcs_mean', 'rcs_mean',
                              'rmsf_A_mean')},
        'watched': {w: DA['arms'][a]['watched'][w] for w in DW}} for a in ARMO},
    'direction': {w: {'dcs_holds': DDIR[w]['dcs']['every_replicate_agrees'],
                      'rcs_holds': DDIR[w]['rcs']['every_replicate_agrees'],
                      'dcs_change': DDIR[w]['dcs']['pooled_change'],
                      'rcs_change': DDIR[w]['rcs']['pooled_change']}
                  for w in DW if w in DDIR},
}

REPO_OWNER = 'karagol-taner'
GH = os.path.join(ROOT, 'github')
_repos = [d for d in sorted(os.listdir(GH))
          if os.path.isdir(os.path.join(GH, d)) and not d.startswith('.')]
if len(_repos) != 1:
    raise SystemExit('STOP: %d repository folders under %s, and the availability '
                     'statement names one: %s' % (len(_repos), GH, _repos))
C407 = os.path.join(D, 'esm2_step8', 'c407_double_rescue.tsv')
if not os.path.isfile(C407):
    raise SystemExit('STOP: %s is absent, so the number of scored pairs at the '
                     'target cannot be quoted' % C407)
with open(C407) as fh:
    _scored = sum(1 for _ in fh) - 1
_expect = N['step8']['n_pairs'] * 19 * 19
if _scored != _expect:
    raise SystemExit('STOP: the step 8 scan holds %d scored pairs and %d '
                     'position pairs, which give %d, so the two do not agree'
                     % (_scored, N['step8']['n_pairs'], _expect))
N['step8']['n_scored_pairs'] = _scored

RGP = os.path.join(D, 'rescue_geometry.json')
if not os.path.isfile(RGP):
    raise SystemExit('STOP: %s is absent, so the rescue geometry cannot be '
                     'quoted. Run codes/36_rescue_geometry.py first' % RGP)
RG = json.load(open(RGP))
if sorted(RG['entries']) != sorted(N['af3']['references']):
    raise SystemExit('STOP: the rescue geometry measures %s and the predicted '
                     'models were compared against %s'
                     % (', '.join(RG['entries']), ', '.join(N['af3']['references'])))
if RG['n_targets'] != len(ac['chosen']):
    raise SystemExit('STOP: the rescue geometry covers %d targets and the '
                     'construct carries %d'
                     % (RG['n_targets'], len(ac['chosen'])))
if RG['window'] != N['coupled']['window']:
    raise SystemExit('STOP: the rescue geometry uses a window of %d and the '
                     'coupling analysis uses %d'
                     % (RG['window'], N['coupled']['window']))
if RG['cutoff_A'] != N['coupled']['reported_cutoff_A']:
    raise SystemExit('STOP: the rescue geometry uses a cutoff of %.1f and the '
                     'coupling analysis reports %.1f'
                     % (RG['cutoff_A'], N['coupled']['reported_cutoff_A']))
N['rescue'] = {
    'pairs': RG['n_pairs'],
    'targets': RG['n_targets'],
    'touching': RG['n_touching'],
    'beyond_window': RG['n_beyond_window'],
    'window': RG['window'],
    'cutoff': RG['cutoff_A'],
    'separation_min': RG['separation_min'],
    'separation_max': RG['separation_max'],
    'distance_min': RG['distance_min'],
    'distance_max': RG['distance_max'],
}

JFP = os.path.join(D, 'junction_figure.json')
if not os.path.isfile(JFP):
    raise SystemExit('STOP: %s is absent, so the junction cannot be described. '
                     'Run codes/33_structure_figures.py first' % JFP)
JF = json.load(open(JFP))
if sorted(JF['watched']) != sorted(N['af3']['watched_positions']):
    raise SystemExit('STOP: the junction figure watches %s and the predicted '
                     'models watch %s' % (JF['watched'],
                                          N['af3']['watched_positions']))
if JF['crossing_six'][1] + 1 != JF['connection'][0]:
    raise SystemExit('STOP: the sixth crossing ends at %d and the connection '
                     'starts at %d, so they do not meet'
                     % (JF['crossing_six'][1], JF['connection'][0]))
if JF['fluctuation_replicates'] != N['md']['run']['replicates']:
    raise SystemExit('STOP: the fluctuation panel averaged %d replicates and '
                     'the simulation reports %d'
                     % (JF['fluctuation_replicates'],
                        N['md']['run']['replicates']))
N['junction'] = {
    'crossing_six': JF['crossing_six'],
    'crossing_six_end': JF['crossing_six'][1],
    'connection': JF['connection'],
    'connection_start': JF['connection'][0],
    'connection_length': JF['connection_length'],
    'replicates': JF['fluctuation_replicates'],
    'fluctuation_low': JF['fluctuation_low'],
    'fluctuation_high': JF['fluctuation_high'],
}

SFP = os.path.join(D, 'structure_figure.json')
if not os.path.isfile(SFP):
    raise SystemExit('STOP: %s is absent, so the structural figure cannot be '
                     'described. Run codes/33_structure_figures.py first' % SFP)
SF = json.load(open(SFP))
if SF['residue_range'] != N['md']['residue_range']:
    raise SystemExit('STOP: the structural figure covers %s and the simulation '
                     'covers %s' % (SF['residue_range'], N['md']['residue_range']))
if SF['reference'] not in N['af3']['references']:
    raise SystemExit('STOP: the structural figure draws %s and the predicted '
                     'models were compared against %s'
                     % (SF['reference'], ', '.join(N['af3']['references'])))
if SF['membrane_crossings'] != 2 * SF['crossings_amino_half']:
    raise SystemExit('STOP: the structural figure finds %d membrane crossings '
                     'and %d of them in the amino terminal half, so the two '
                     'halves do not carry the same number'
                     % (SF['membrane_crossings'], SF['crossings_amino_half']))
N['structure'] = {
    'crossings': SF['membrane_crossings'],
    'crossings_per_half': SF['crossings_amino_half'],
    'split_residue': SF['split_residue'],
    'half_thickness': SF['half_thickness_A'],
    'plddt_floor': SF['plddt_floor'],
    'plddt_below_floor': SF['plddt_below_floor'],
    'site_shell': SF['site_shell_A'],
    'site_positions': SF['site_positions_shown'],
    'site_humanising': SF['site_positions_humanising'],
    'site_rescuers': SF['site_positions_rescuer'],
    'reference': SF['reference'],
}

N['repo'] = {'owner': REPO_OWNER, 'name': _repos[0],
             'url': 'github.com/%s/%s' % (REPO_OWNER, _repos[0])}

out = os.path.join(D, 'derived_values.json')
with open(out, 'w') as fh:
    json.dump(N, fh, indent=1, sort_keys=True)

missing = []


def walk(d, path=()):
    if isinstance(d, dict):
        for k, v in d.items():
            walk(v, path + (str(k),))
    elif d is None:
        missing.append('.'.join(path))


walk(N)
print('wrote %s (%d bytes)' % (out, os.path.getsize(out)))
print('sections: %s' % ', '.join(sorted(N)))
if missing:
    print('\nVALUES NOT FOUND:')
    for m in missing:
        print('  ' + m)
