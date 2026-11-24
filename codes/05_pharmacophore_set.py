import json
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict

import gemmi

EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
TOOL = 'tool=pcft-position-set&email=karagolt@mit.edu'

m = json.load(open('structures/7BC6_to_human.json'))
ck2hu = {int(k): v for k, v in m['ck2hu'].items()}
hu2ck = {int(k): v for k, v in m['hu2ck'].items()}

human = ''.join(l.strip() for l in open('seq/Q96NT5.fasta') if not l.startswith('>'))
ent = {}
for rec in open('seq/7BC6.fasta').read().split('>'):
    if rec.strip():
        h, *r = rec.strip().split('\n')
        ent[h.split('|')[0]] = ''.join(r)
chicken = ent['7BC6_1'][:473]

_pr = json.load(open('data/parker_residues.json'))['residues']
PARKER = {int(k): v['human'] for k, v in _pr.items()}
PARKER_ROLE = {int(k): v['role'] for k, v in _pr.items()}

st = gemmi.read_structure('seq/7BC7.cif')
st.setup_entities()
st.remove_waters()
ns = gemmi.NeighborSearch(st, 6.0, 0).populate()
model = st[0]
lig = [r for ch in model for r in ch if r.name == 'LYA']


def shell(cutoff):
    found = set()
    for lres in lig:
        for atom in lres:
            for mark in ns.find_atoms(atom.pos, '\0', radius=cutoff):
                cra = mark.to_cra(model)
                if cra.residue.name in ('LYA', 'HOH') or cra.chain.name != 'A':
                    continue
                if cra.atom.pos.dist(atom.pos) <= cutoff:
                    found.add(cra.residue.seqid.num)
    return sorted(found)


rows = [l.rstrip('\n').split('\t') for l in open('data/hfm_positions.tsv')][1:]
pos_subs = {int(r[0]): (r[1], r[2].split(','), r[3].split(',')) for r in rows}

all_pmids = sorted({p for _, (_, _, pm) in pos_subs.items() for p in pm}, key=int)


def fetch_titles(pmids):
    titles = {}
    for i in range(0, len(pmids), 100):
        url = (EUTILS + 'efetch.fcgi?db=pubmed&retmode=xml&' + TOOL
               + '&id=' + ','.join(pmids[i:i + 100]))
        for attempt in range(4):
            try:
                with urllib.request.urlopen(url, timeout=60) as r:
                    blob = r.read()
                break
            except Exception:
                if attempt == 3:
                    raise
                time.sleep(2 + 2 * attempt)
        root = ET.fromstring(blob)
        for art in root.iter('PubmedArticle'):
            pm = art.find('.//PMID').text
            t = art.find('.//ArticleTitle')
            titles[pm] = ''.join(t.itertext()) if t is not None else ''
        time.sleep(0.4)
    return titles


titles = fetch_titles(all_pmids)

MUTAGENESIS = re.compile(
    r'cysteine|cys-|substituted[- ]cysteine|scanning|accessibility|cross-link|'
    r'crosslink|alanine|site-directed|mutagenesis|topolog|oligomer|'
    r'structure-?function|residues? (?:in|of|that|required)|'
    r'transmembrane domain|proton coupling|substrate binding|'
    r'glutamate|histidine|tryptophan|serine|role of',
    re.I)
CLINICAL = re.compile(
    r'patient|subject|case|famil|infant|child|novel mutation|'
    r'hereditary folate malabsorption in|clinical|diagnos|report of|'
    r'identification of (?:a|two|three|novel)|syndrome|phenotype|cohort',
    re.I)

kind = {}
for pm, t in titles.items():
    cl, mu = bool(CLINICAL.search(t)), bool(MUTAGENESIS.search(t))
    kind[pm] = 'clinical' if (cl and not mu) else ('mutagenesis' if mu and not cl
                                                   else ('both' if cl and mu else 'unclear'))

print('=== source papers behind the extracted positions ===')
print('pmid       class         title')
for pm in all_pmids:
    print(f'{pm:10s} {kind[pm]:13s} {titles[pm][:110]}')

H, L = {}, {}
for pos, (wt, subs, pms) in sorted(pos_subs.items()):
    classes = {kind[p] for p in pms}
    patient_subs = [s for s in subs if s not in ('C', 'A')]
    if ('clinical' in classes or 'both' in classes) and patient_subs:
        H[pos] = (wt, patient_subs, [p for p in pms])
    else:
        L[pos] = (wt, subs, [p for p in pms])

print(f'\n=== component H, hereditary folate malabsorption patient missense '
      f'positions: {len(H)} ===')
for pos, (wt, subs, pms) in sorted(H.items()):
    print(f'  {wt}{pos} -> {",".join(sorted(subs))}   pmids {",".join(sorted(pms, key=int))}')

print(f'\n=== component L, laboratory substitution positions: {len(L)} ===')
print('  ' + ', '.join(f'{v[0]}{p}' for p, v in sorted(L.items())))


def c_count(human_positions, label):
    rows, miss = [], []
    for hp in sorted(human_positions):
        cp = hu2ck.get(hp)
        if cp is None:
            miss.append(hp)
            continue
        caa, haa = chicken[cp - 1], human[hp - 1]
        rows.append((hp, haa, cp, caa, caa != haa))
    diff = [r for r in rows if r[4]]
    print(f'\n--- {label} ---')
    print(f'positions in set: {len(human_positions)}, mapped: {len(rows)}, '
          f'unmapped: {len(miss)}'
          + (f' {miss}' if miss else ''))
    print(f'C = {len(diff)} of {len(rows)} mapped positions differ '
          f'({100.0 * len(diff) / len(rows):.0f}%)' if rows else 'C = 0')
    if diff:
        print('  differing: ' + ', '.join(
            f'chick {caa}{cp} / human {haa}{hp}' for hp, haa, cp, caa, _ in diff))
    return len(diff), len(rows)


print('\n' + '=' * 72)
print('C COUNTS')
print('=' * 72)

P_hu = sorted(PARKER.values())
c_count(P_hu, 'set P: Parker et al. 2021 named functional residues (13)')

for cutoff in (4.0, 4.5, 5.0):
    sh = shell(cutoff)
    sh_hu = [ck2hu[c] for c in sh if c in ck2hu]
    c_count(sh_hu, f'set S: pemetrexed contact shell at {cutoff} A '
                   f'({len(sh)} chicken positions)')

c_count(sorted(H), f'set H: HFM patient missense positions ({len(H)})')

for cutoff in (4.0, 4.5, 5.0):
    sh = shell(cutoff)
    union = sorted(set(P_hu) | {ck2hu[c] for c in sh if c in ck2hu} | set(H))
    c_count(union, f'UNION P + S({cutoff} A) + H  = the pharmacophore position set')
