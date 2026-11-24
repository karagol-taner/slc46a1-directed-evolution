import json
import os
import sys

import gemmi

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = os.path.join(ROOT, 'md', 'systems')
SCRATCH = os.path.join(ROOT, 'md', 'systems', '.wt_check.pdb')
ADM = os.path.join(ROOT, 'data', 'admissible_construct.json')
POSSET = os.path.join(ROOT, 'data', 'pharmacophore_position_set.tsv')
SEQF = os.path.join(ROOT, 'seq', '7BC6.fasta')
CIF = os.path.join(ROOT, 'seq', '7BC6.cif')
MAPOUT = os.path.join(ROOT, 'md', 'construct_systems.json')
LOG = os.path.join(ROOT, 'data', '24_md_construct_systems_log.txt')
STATE = 'apo_pH7.5'
ENTRY = '7BC6'
BACKBONE = ['N', 'CA', 'C', 'O']
THREE = {'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE', 'G': 'GLY',
         'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
         'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER', 'T': 'THR', 'V': 'VAL',
         'W': 'TRP', 'Y': 'TYR'}
ONE = {v: k for k, v in THREE.items()}

lines = []


def log(m):
    print(m, flush=True)
    lines.append(str(m))


def die(m):
    with open(LOG, 'w') as f:
        f.write('\n'.join(lines + ['FATAL ' + m]) + '\n')
    sys.exit('FATAL ' + m)


for p in (ADM, POSSET, SEQF, CIF):
    if not os.path.exists(p):
        die('missing input %s' % p)

ent = {}
for rec in open(SEQF).read().split('>'):
    if rec.strip():
        parts = rec.strip().split('\n')
        ent[parts[0].split('|')[0]] = ''.join(parts[1:])
if '7BC6_1' not in ent:
    die('%s holds no 7BC6_1 record' % SEQF)
chicken = ent['7BC6_1'][:473]

rows = [l.rstrip('\n').split('\t') for l in open(POSSET)]
pc = {k: i for i, k in enumerate(rows[0])}
PSET = {int(r[pc['chicken_pos']]): {'chicken_aa': r[pc['chicken_aa']],
                                    'human_aa': r[pc['human_aa']],
                                    'differs': r[pc['differs']] == '1'}
        for r in rows[1:]}
PRIM = sorted(p for p, v in PSET.items() if v['differs'])

adm = json.load(open(ADM))
if sorted(adm['primaries']) != PRIM:
    die('the construct file lists different primaries from the position set')
RESC = {int(k): v for k, v in adm['rescuer_positions'].items()}
for q in RESC:
    if q in PSET:
        die('rescuer position %d is inside the pharmacophore set' % q)

SUBS = {p: PSET[p]['human_aa'] for p in PRIM}
SUBS.update(RESC)
if len(SUBS) != adm['n_positions_changed']:
    die('%d substitution sites assembled, the construct file records %d changed positions'
        % (len(SUBS), adm['n_positions_changed']))
target = list(chicken)
for p, y in SUBS.items():
    target[p - 1] = y
if ''.join(target) != adm['sequence']:
    n = sum(1 for a, b in zip(''.join(target), adm['sequence']) if a != b)
    die('the substitutions assembled here do not reproduce the construct sequence, %d '
        'positions differ' % n)
log('substitutions assembled from %s and checked against its own sequence: %d sites, '
    '%d primaries and %d rescuers' % (os.path.basename(ADM), len(SUBS), len(PRIM),
                                      len(RESC)))

SYSTEMS = {'wt': {}, 'primaries': {p: PSET[p]['human_aa'] for p in PRIM}, 'construct': SUBS}
EXPECT_N = {'wt': 0, 'primaries': len(PRIM), 'construct': len(SUBS)}


def keep_for(target_one):
    if THREE[target_one] in ('GLY', 'PRO'):
        return list(BACKBONE)
    return BACKBONE + ['CB']


def build(sites, path):
    st = gemmi.read_structure(CIF)
    st.setup_entities()
    st.remove_waters()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    model = st[0]
    chain = None
    for ch in model:
        if len(ch.get_polymer()) > 300:
            chain = ch
            break
    if chain is None:
        die('%s has no chain with more than 300 residues' % ENTRY)
    res = {r.seqid.num: r for r in chain}
    done = 0
    for cp, to_one in sorted(sites.items()):
        if cp not in res:
            die('%s chain %s has no residue %d' % (ENTRY, chain.name, cp))
        r = res[cp]
        if r.name not in ONE:
            die('%s %s%d is %s, not a standard residue' % (ENTRY, chain.name, cp, r.name))
        if ONE[r.name] != chicken[cp - 1]:
            die('%s %s%d is %s, the chicken sequence has %s at that index, so the author '
                'numbering and the sequence index do not agree there'
                % (ENTRY, chain.name, cp, ONE[r.name], chicken[cp - 1]))
        if to_one == chicken[cp - 1]:
            die('the substitution at chicken %d writes the residue already there' % cp)
        keep = keep_for(to_one)
        before = {a.name for a in r}
        if not set(BACKBONE) <= before:
            die('%s %s%d is missing backbone atoms %s'
                % (ENTRY, chain.name, cp, sorted(set(BACKBONE) - before)))
        for nm in [a.name for a in r if a.name not in keep]:
            r.remove_atom(nm, '*')
        r.name = THREE[to_one]
        after = {a.name for a in r}
        if not set(BACKBONE) <= after:
            die('%s %s%d lost backbone atoms during substitution' % (ENTRY, chain.name, cp))
        if after - set(keep):
            die('%s %s%d kept unexpected atoms %s'
                % (ENTRY, chain.name, cp, sorted(after - set(keep))))
        done += 1
    out = gemmi.Structure()
    out.spacegroup_hm = st.spacegroup_hm
    out.cell = st.cell
    m = gemmi.Model('1')
    m.add_chain(chain)
    out.add_model(m)
    out.setup_entities()
    out.write_pdb(path)
    return done, chain.name, len(res)


os.makedirs(OUTDIR, exist_ok=True)
ref = os.path.join(OUTDIR, 'wt__%s.pdb' % STATE)
if not os.path.exists(ref):
    die('no wild-type structure on disk at %s to check this script against' % ref)
n, cname, nres = build({}, SCRATCH)
with open(SCRATCH, 'rb') as a, open(ref, 'rb') as b:
    same = a.read() == b.read()
os.remove(SCRATCH)
if not same:
    die('the wild-type arm written here differs from %s, so this script does not reproduce '
        'the route codes/09_md_systems.py used. Nothing was overwritten.'
        % os.path.basename(ref))
log('route check: the wild-type arm reproduces %s byte for byte, chain %s, %d residues'
    % (os.path.basename(ref), cname, nres))

written = {}
for name in ('primaries', 'construct'):
    path = os.path.join(OUTDIR, '%s__%s.pdb' % (name, STATE))
    n, cname, nres = build(SYSTEMS[name], path)
    if n != EXPECT_N[name]:
        die('%s applied %d substitutions, %d expected' % (name, n, EXPECT_N[name]))
    st = gemmi.read_structure(path)
    got = {r.seqid.num: ONE.get(r.name) for r in st[0][0]}
    for cp, to_one in SYSTEMS[name].items():
        if got.get(cp) != to_one:
            die('%s carries %s at chicken %d, %s intended'
                % (name, got.get(cp), cp, to_one))
    for cp in set(PSET) - set(SYSTEMS[name]):
        if cp in got and got[cp] != PSET[cp]['chicken_aa']:
            die('%s changed chicken %d, which it should not touch' % (name, cp))
    written[name] = {'file': os.path.basename(path),
                     'substitutions': {str(k): v for k, v in sorted(SYSTEMS[name].items())},
                     'n_substitutions': n,
                     'bytes': os.path.getsize(path)}
    log('wrote %s, %d substitutions, %d bytes'
        % (os.path.basename(path), n, os.path.getsize(path)))

st = gemmi.read_structure(CIF)
st.setup_entities()
st.remove_waters()
chain = next((c for c in st[0] if len(c.get_polymer()) > 300), None)
if chain is None:
    die('%s has no chain with more than 300 residues' % ENTRY)
obs = {r.seqid.num for r in chain}
N = len(chicken)
gaps, run = [], []
for i in range(1, N + 1):
    if i not in obs:
        run.append(i)
    elif run:
        gaps.append([run[0], run[-1]])
        run = []
if run:
    gaps.append([run[0], run[-1]])
internal = [g for g in gaps if g[0] > 1 and g[1] < N]
coverage = {'length': N, 'observed': len(obs), 'unmodelled': N - len(obs),
            'gaps': gaps, 'internal_gaps': internal,
            'set_positions': len(PSET),
            'set_observed': len(set(PSET) & obs),
            'set_missing': sorted(set(PSET) - obs),
            'substituted': len(SUBS),
            'substituted_observed': len(set(SUBS) & obs),
            'substituted_missing': sorted(set(SUBS) - obs)}
if coverage['substituted_missing']:
    die('substituted positions %s are not observed in %s, so they cannot be threaded'
        % (coverage['substituted_missing'], ENTRY))
log('')
log('coverage of %s chain %s: %d of %d residues observed, %d unmodelled in %d stretches '
    '%s, of which %d internal'
    % (ENTRY, chain.name, coverage['observed'], N, coverage['unmodelled'], len(gaps),
       ', '.join('%d to %d' % tuple(g) if g[0] != g[1] else str(g[0]) for g in gaps),
       len(internal)))
log('  pharmacophore set: %d of %d observed, missing %s'
    % (coverage['set_observed'], coverage['set_positions'],
       coverage['set_missing'] or 'none'))
log('  substituted sites: %d of %d observed'
    % (coverage['substituted_observed'], coverage['substituted']))

written['wt'] = {'file': os.path.basename(ref), 'substitutions': {},
                 'n_substitutions': 0, 'bytes': os.path.getsize(ref),
                 'note': 'written by codes/09_md_systems.py and reproduced by this script'}
json.dump({'state': STATE, 'entry': ENTRY, 'ligand_kept': False,
           'coverage': coverage,
           'truncation': 'backbone, plus CB where the target residue has one',
           'arms': written,
           'construct_sequence': adm['sequence'],
           'rescuer_positions': {str(k): v for k, v in sorted(RESC.items())},
           'primaries': PRIM,
           'next': 'upload each file to CHARMM-GUI Membrane Builder with the settings in '
                   'md/md_plan.json, then run the GROMACS notebook on the step3 output'},
          open(MAPOUT, 'w'), indent=1)
log('')
log('three arms for the %s state: %s' % (STATE, ', '.join(sorted(written))))
log('wrote ' + os.path.relpath(MAPOUT, ROOT))
with open(LOG, 'w') as f:
    f.write('\n'.join(lines) + '\n')
