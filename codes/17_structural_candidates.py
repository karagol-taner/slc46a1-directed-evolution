import json
import os
import sys

import gemmi

STRUCTS = {'7BC6': 'seq/7BC6.cif', '7BC7': 'seq/7BC7.cif'}
SWEEP = 'data/rescue_sweep_report.json'
POSSET = 'data/pharmacophore_position_set.tsv'
TOPN = 10
OUT = 'data/structural_candidates.json'
LOG = 'data/17_structural_candidates_log.txt'

lines = []


def log(m):
    print(m, flush=True)
    lines.append(str(m))


def die(m):
    with open(LOG, 'w') as f:
        f.write('\n'.join(lines + ['FATAL ' + m]) + '\n')
    sys.exit('FATAL ' + m)


for p in [SWEEP, POSSET] + list(STRUCTS.values()):
    if not os.path.exists(p):
        die(p + ' not found')

ent = {}
for rec in open('seq/7BC6.fasta').read().split('>'):
    if rec.strip():
        parts = rec.strip().split('\n')
        ent[parts[0].split('|')[0]] = ''.join(parts[1:])
chicken = ent['7BC6_1'][:473]
if len(chicken) != 473:
    die('chicken native length %d, expected 473' % len(chicken))

atoms = {}
for name, path in STRUCTS.items():
    st = gemmi.read_structure(path)
    st.setup_entities()
    st.remove_ligands_and_waters()
    poly = st[0]['A'].get_polymer()
    per = {}
    bad = []
    for r in poly:
        n = r.seqid.num
        one = gemmi.find_tabulated_residue(r.name)
        letter = one.one_letter_code.upper() if one else 'X'
        if not 1 <= n <= len(chicken):
            continue
        if letter != chicken[n - 1]:
            bad.append((n, letter, chicken[n - 1]))
        per[n] = [a.pos for a in r if a.element != gemmi.Element('H')]
    if bad:
        die('%s author numbering disagrees with the chicken sequence at %d residues, '
            'the first being %d where the structure has %s and the sequence has %s'
            % (name, len(bad), bad[0][0], bad[0][1], bad[0][2]))
    if not per:
        die(name + ' yielded no chain A residues')
    atoms[name] = per
    ks = sorted(per)
    log('%s: %d modelled chain A residues, author %d to %d, every residue agrees with '
        'the chicken sequence' % (name, len(per), ks[0], ks[-1]))

common = sorted(set(atoms['7BC6']) & set(atoms['7BC7']))
log('modelled in both structures: %d positions' % len(common))

sweep = json.load(open(SWEEP))
targets = sorted(int(t) for t in sweep['per_target'])
rows = [l.rstrip('\n').split('\t') for l in open(POSSET)]
col = {k: i for i, k in enumerate(rows[0])}
by_chicken = {int(r[col['chicken_pos']]): r for r in rows[1:]}

out = {}
for t in targets:
    if t not in common:
        log('')
        log('chicken %d is not modelled in both structures, so no distance ranking is '
            'possible for it' % t)
        out[str(t)] = {'modelled': False, 'candidates': []}
        continue
    scored = []
    for q in common:
        if q == t:
            continue
        d6 = min(a.dist(b) for a in atoms['7BC6'][t] for b in atoms['7BC6'][q])
        d7 = min(a.dist(b) for a in atoms['7BC7'][t] for b in atoms['7BC7'][q])
        scored.append({'chicken_pos': q, 'separation': abs(q - t),
                       'min_dist_7BC6': d6, 'min_dist_7BC7': d7,
                       'min_dist': min(d6, d7),
                       'in_position_set': q in by_chicken,
                       'is_primary': bool(q in by_chicken
                                          and by_chicken[q][col['differs']] == '1')})
    scored.sort(key=lambda r: r['min_dist'])
    if len(scored) < TOPN:
        die('chicken %d has only %d ranked positions, %d requested' % (t, len(scored), TOPN))
    top = scored[:TOPN]
    rec = by_chicken.get(t)
    log('')
    log('chicken %d %s to human %s: %d positions ranked by distance'
        % (t, chicken[t - 1],
           rec[col['human_aa']] if rec else '?', len(scored)))
    log('  pos  sep  min dist  7BC6   7BC7  in the set  primary')
    for r in top:
        log('  %4d %4d %9.2f %6.2f %6.2f %11s %8s'
            % (r['chicken_pos'], r['separation'], r['min_dist'],
               r['min_dist_7BC6'], r['min_dist_7BC7'],
               'yes' if r['in_position_set'] else 'no',
               'yes' if r['is_primary'] else 'no'))
    out[str(t)] = {'modelled': True, 'n_ranked': len(scored),
                   'top_n': TOPN, 'candidates': top}

json.dump({'structures': sorted(STRUCTS), 'top_n': TOPN,
           'modelled_in_both': len(common), 'targets': out},
          open(OUT, 'w'), indent=1)
log('')
log('wrote ' + OUT)
with open(LOG, 'w') as f:
    f.write('\n'.join(lines) + '\n')
