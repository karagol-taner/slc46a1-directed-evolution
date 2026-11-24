import itertools
from collections import defaultdict
from Bio import Align
from Bio.Align import substitution_matrices


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


human = read_fasta_single('seq/Q96NT5.fasta')
chicken = read_entity('seq/7BC6.fasta', '7BC6_1')[:473]

PARKER = [43, 98, 135, 156, 164, 193, 289, 290, 323, 350, 384, 407, 411]
CONTACT_5A = [37, 40, 44, 98, 135, 156, 164, 165, 169, 193, 196, 286, 290, 293,
              320, 323, 380, 384, 407, 411, 415]
PROBE = sorted(set(PARKER + CONTACT_5A + [419]))

matrices = ['BLOSUM62', 'BLOSUM80', 'PAM250']
opens = [-7, -9, -11, -13, -15]
extends = [-0.5, -1, -2]

seen = defaultdict(set)
n_settings = 0
n_paths = 0
for mat, go, ge in itertools.product(matrices, opens, extends):
    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load(mat)
    al.open_gap_score = go
    al.extend_gap_score = ge
    al.mode = 'global'
    alns = al.align(chicken, human)
    n_settings += 1
    npaths = min(len(alns), 20)
    for k in range(npaths):
        a = alns[k]
        n_paths += 1
        ci = hi = 0
        m = {}
        for x, y in zip(a[0], a[1]):
            cn = hn = None
            if x != '-':
                ci += 1
                cn = ci
            if y != '-':
                hi += 1
                hn = hi
            if cn and hn:
                m[cn] = hn
        for p in PROBE:
            seen[p].add(m.get(p))

print(f'settings swept: {n_settings} '
      f'({len(matrices)} matrices x {len(opens)} gap-open x {len(extends)} gap-extend)')
print(f'alignment paths examined: {n_paths}')

stable, unstable = [], []
for p in PROBE:
    vals = seen[p]
    if len(vals) == 1:
        stable.append((p, next(iter(vals))))
    else:
        unstable.append((p, sorted(v for v in vals if v is not None),
                         None in vals))

print(f'\nunambiguous positions: {len(stable)}/{len(PROBE)}')
print(f'positions whose human partner changes with the parameters: {len(unstable)}')
for p, vals, hasgap in unstable:
    print(f'  chicken {p}: human partners seen {vals}'
          + (' (and gapped in some alignments)' if hasgap else ''))

print('\nchicken -> human, every setting agreeing:')
for p, h in stable:
    print(f'  {chicken[p-1]}{p} -> {human[h-1]}{h}'
          f'{"   DIFFERENT" if chicken[p-1] != human[h-1] else ""}')
