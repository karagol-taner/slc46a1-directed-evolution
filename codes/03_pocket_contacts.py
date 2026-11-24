import gemmi

ent = {}
for rec in open('seq/7BC6.fasta').read().split('>'):
    if not rec.strip():
        continue
    head, *rest = rec.strip().split('\n')
    ent[head.split('|')[0]] = ''.join(rest)
chainA_entity = ent['7BC6_1']
print(f'deposited entity 7BC6_1 (chain A) length: {len(chainA_entity)}')

st = gemmi.read_structure('seq/7BC6.cif')
st.setup_entities()
poly = st[0]['A'].get_polymer()
obs = gemmi.one_letter_code(poly.extract_sequence()).upper()
nums = [r.seqid.num for r in poly]
print(f'observed chain A: {len(obs)} residues, author {nums[0]}..{nums[-1]}')

idx = chainA_entity.find(obs)
print(f'observed substring found in entity at 0-based index {idx} '
      f'(1-based entity position {idx + 1})')
print(f'=> author number {nums[0]} corresponds to entity position {idx + 1}; '
      f'offset author - entity = {nums[0] - (idx + 1)}')

with open('data/7BC6_chainA_authornum.tsv', 'w') as fh:
    fh.write('author_num\taa\n')
    for n, a in zip(nums, obs):
        fh.write(f'{n}\t{a}\n')
print('wrote 7BC6_chainA_authornum.tsv')

st7 = gemmi.read_structure('seq/7BC7.cif')
st7.setup_entities()
st7.remove_waters()
ns = gemmi.NeighborSearch(st7, 6.0, 0).populate()
model = st7[0]

lig_atoms = []
for ch in model:
    for r in ch:
        if r.name == 'LYA':
            lig_atoms.append((ch.name, r.seqid.num, r))
print(f'\nLYA copies in 7BC7: {[(c, n) for c, n, _ in lig_atoms]}')

for cutoff in (4.0, 4.5, 5.0):
    hits = {}
    for cname, lnum, lres in lig_atoms:
        for atom in lres:
            marks = ns.find_atoms(atom.pos, '\0', radius=cutoff)
            for m in marks:
                cra = m.to_cra(model)
                if cra.residue.name in ('LYA', 'HOH'):
                    continue
                if cra.chain.name != 'A':
                    continue
                d = cra.atom.pos.dist(atom.pos)
                if d <= cutoff:
                    key = (cra.residue.seqid.num,
                           gemmi.find_tabulated_residue(cra.residue.name).one_letter_code.upper())
                    hits[key] = min(hits.get(key, 99), round(d, 2))
    ordered = sorted(hits.items())
    print(f'\n--- chain A residues within {cutoff} A of pemetrexed: {len(ordered)} ---')
    print(', '.join(f'{aa}{num}({d})' for (num, aa), d in ordered))
