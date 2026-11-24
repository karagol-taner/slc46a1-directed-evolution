import json
from Bio import Align
from Bio.Align import substitution_matrices

SEQ_DIR = 'seq/'


def read_fasta_single(path):
    return ''.join(l.strip() for l in open(path) if not l.startswith('>'))


def read_entity(path, key):
    ent = {}
    for rec in open(path).read().split('>'):
        if not rec.strip():
            continue
        head, *rest = rec.strip().split('\n')
        ent[head.split('|')[0]] = ''.join(rest)
    return ent[key]


human = read_fasta_single(SEQ_DIR + 'Q96NT5.fasta')
construct = read_entity(SEQ_DIR + '7BC6.fasta', '7BC6_1')

TAG = 'GGENLYFQ'
assert construct.endswith('ENLYFQ'), construct[-12:]
chicken = construct[:473]
assert len(chicken) == 473

print(f'human  Q96NT5        : {len(human)} aa')
print(f'chicken 7BC6 chain A : {len(chicken)} aa native part '
      f'(construct entity {len(construct)} aa, C-terminal tag removed)')

aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
aligner.open_gap_score = -11
aligner.extend_gap_score = -1
aligner.mode = 'global'
aln = aligner.align(chicken, human)[0]
print(f'\nalignment score {aln.score}, {len(aligner.align(chicken, human))} optimal path(s)')

ca, ha = aln[0], aln[1]
ck2hu, hu2ck = {}, {}
ci = hi = 0
ident = aligned = 0
cols = []
for a, b in zip(ca, ha):
    cn = hn = None
    if a != '-':
        ci += 1
        cn = ci
    if b != '-':
        hi += 1
        hn = hi
    if cn and hn:
        aligned += 1
        if a == b:
            ident += 1
        ck2hu[cn] = hn
        hu2ck[hn] = cn
    cols.append((cn, a, hn, b))

print(f'aligned columns {aligned}, identical {ident} '
      f'({100.0 * ident / aligned:.1f}% over aligned positions)')
print(f'identity over full chicken length: {100.0 * ident / len(chicken):.1f}%')

PARKER_PAIRS = {
    156: 148,
    164: 156,
    289: 281,
    193: 185,
    384: 376,
    98:  90,
    411: 403,
    290: 282,
    323: 315,
    407: 399,
    350: 342,
    43:  39,
    135: 127,
}
print('\n--- validation: Parker et al. 2021 stated pairs vs this alignment ---')
print('chick  stated_human  aligned_human  chick_aa  human_aa  offset  agree')
ok = 0
for cn in sorted(PARKER_PAIRS):
    stated = PARKER_PAIRS[cn]
    got = ck2hu.get(cn)
    caa = chicken[cn - 1]
    haa = human[got - 1] if got else '-'
    agree = (got == stated)
    ok += agree
    print(f'{cn:5d}  {stated:12d}  {str(got):13s}  {caa:8s}  {haa:8s}'
          f'  {("" if got is None else cn - got):6}  {"yes" if agree else "NO"}')
print(f'{ok}/{len(PARKER_PAIRS)} stated pairs reproduced by the alignment')

offs = sorted({cn - PARKER_PAIRS[cn] for cn in PARKER_PAIRS})
print(f'distinct chicken-minus-human offsets among the stated pairs: {offs}')

runs = []
prev = None
for cn in sorted(ck2hu):
    o = cn - ck2hu[cn]
    if o != prev:
        runs.append((cn, o))
        prev = o
print('\noffset changes along the chain (chicken position where offset becomes value):')
print(', '.join(f'{p}->{o}' for p, o in runs))

json.dump({'ck2hu': ck2hu, 'hu2ck': hu2ck,
           'identity_aligned_pct': round(100.0 * ident / aligned, 2),
           'aligned_columns': aligned, 'identical': ident,
           'parker_pairs_reproduced': f'{ok}/{len(PARKER_PAIRS)}'},
          open('structures/7BC6_to_human.json', 'w'), indent=0)
print('\nwrote pcft_residue_map.json')

with open('data/pcft_alignment.txt', 'w') as fh:
    fh.write(str(aln))
print('wrote pcft_alignment.txt')
