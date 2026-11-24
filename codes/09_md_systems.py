import json
import os
import sys

import gemmi

OUTDIR = 'md/systems'
MUTATIONS = {407: ('GLU', 'ASN'), 166: ('ASN', 'GLY')}
SYSTEMS = {
    'wt': [],
    'E407N': [407],
    'N166G': [166],
    'E407N_N166G': [407, 166],
}
STATES = {'apo_pH7.5': ('7BC6', False), 'pmx_pH6.0': ('7BC7', True)}
BACKBONE = ['N', 'CA', 'C', 'O']
KEEP = {('GLU', 'ASN'): BACKBONE + ['CB', 'CG'],
        ('ASN', 'GLY'): BACKBONE}


def die(msg):
    sys.exit('FATAL ' + msg)


pos = [l.rstrip('\n').split('\t') for l in open('data/pharmacophore_position_set.tsv')]
pc = {k: i for i, k in enumerate(pos[0])}
expect = {}
for r in pos[1:]:
    if r[pc['chicken_pos']]:
        expect[int(r[pc['chicken_pos']])] = (r[pc['chicken_aa']], r[pc['human_aa']])
one = {'GLU': 'E', 'ASN': 'N', 'GLY': 'G'}
for cp, (wt3, mut3) in MUTATIONS.items():
    if cp not in expect:
        die('chicken %d is not in the position set' % cp)
    if expect[cp][0] != one[wt3] or expect[cp][1] != one[mut3]:
        die('position set has chicken %s / human %s at %d, this script assumes %s to %s'
            % (expect[cp][0], expect[cp][1], cp, one[wt3], one[mut3]))

os.makedirs(OUTDIR, exist_ok=True)
written = []

for state, (entry, keep_lig) in STATES.items():
    for name, sites in SYSTEMS.items():
        st = gemmi.read_structure(f'seq/{entry}.cif')
        st.setup_entities()
        st.remove_waters()
        st.remove_alternative_conformations()
        st.remove_hydrogens()
        model = st[0]

        chain = None
        for ch in model:
            poly = ch.get_polymer()
            if len(poly) > 300:
                chain = ch
                break
        if chain is None:
            die('%s has no chain with more than 300 residues' % entry)

        res = {r.seqid.num: r for r in chain}
        for cp in sites:
            wt3, mut3 = MUTATIONS[cp]
            if cp not in res:
                die('%s chain %s has no residue %d' % (entry, chain.name, cp))
            r = res[cp]
            if r.name != wt3:
                die('%s %s%d is %s, expected %s' % (entry, chain.name, cp, r.name, wt3))
            keep = KEEP[(wt3, mut3)]
            before = {a.name for a in r}
            if not set(BACKBONE) <= before:
                die('%s %s%d is missing backbone atoms %s'
                    % (entry, chain.name, cp, sorted(set(BACKBONE) - before)))
            for nm_ in [a.name for a in r if a.name not in keep]:
                r.remove_atom(nm_, '*')
            r.name = mut3
            after = {a.name for a in r}
            if not set(BACKBONE) <= after:
                die('%s %s%d lost backbone atoms during substitution' % (entry, chain.name, cp))
            if after - set(keep):
                die('%s %s%d kept unexpected atoms %s'
                    % (entry, chain.name, cp, sorted(after - set(keep))))

        out = gemmi.Structure()
        out.spacegroup_hm = st.spacegroup_hm
        out.cell = st.cell
        m = gemmi.Model('1')
        m.add_chain(chain)
        if keep_lig:
            lig = gemmi.Chain('L')
            found = 0
            for ch in model:
                for r in ch:
                    if r.name == 'LYA':
                        lig.add_residue(r)
                        found += 1
            if found == 0:
                die('%s has no LYA to keep for the bound state' % entry)
            m.add_chain(lig)
        out.add_model(m)
        out.setup_entities()

        path = os.path.join(OUTDIR, f'{name}__{state}.pdb')
        out.write_pdb(path)
        n_atom = sum(1 for ch in out[0] for r in ch for a in r)
        written.append((path, name, state, len(sites), n_atom, os.path.getsize(path)))

print(f'{len(written)} structures written to {OUTDIR}/')
print('file                              system        state        muts  atoms   bytes')
for p, name, state, nm, na, sz in written:
    print(f'{os.path.basename(p):34s}{name:14s}{state:13s}{nm:4d}{na:8d}{sz:8d}')

ref = {(n, s): a for _, n, s, _, a, _ in written}
for state in STATES:
    base = ref[('wt', state)]
    for name in SYSTEMS:
        d = base - ref[(name, state)]
        exp = sum(len(set(KEEP_ALL := {'GLU': ['CD', 'OE1', 'OE2'],
                                       'ASN': ['CB', 'CG', 'OD1', 'ND2']}[MUTATIONS[cp][0]]))
                  for cp in SYSTEMS[name])
        if d != exp:
            die('%s %s lost %d atoms against the wild type, expected %d'
                % (name, state, d, exp))
print('atom-count differences against each wild type match the truncations exactly')

plan = {
    'purpose': 'test the Glu407 with Asn166 coupling that Parker and colleagues '
               'quantified by alchemical transformation',
    'systems': {k: {'mutations': [f'{MUTATIONS[c][0]}{c}{MUTATIONS[c][1]}' for c in v]}
                for k, v in SYSTEMS.items()},
    'states': {k: {'entry': v[0], 'ligand_kept': v[1]} for k, v in STATES.items()},
    'measurements_chicken_numbering': [
        ['E407N_to_N166', 407, 'ND2', 166, 'OD1'],
        ['E407_to_N166_wildtype', 407, 'OE1', 166, 'ND2'],
        ['R156_to_D164_control', 156, 'NH1', 164, 'OD1'],
        ['E193_to_R384_control', 193, 'OE1', 384, 'NH1'],
    ],
    'measurement_note': 'residue indices are chicken author numbering as deposited. '
                        'CHARMM-GUI renumbers, so resolve these against step3_input.pdb '
                        'before putting them in MEASUREMENTS in the analysis template.',
    'files': [os.path.basename(p) for p, *_ in written],
}
json.dump(plan, open('md/md_plan.json', 'w'), indent=1)
print('wrote md/md_plan.json')
