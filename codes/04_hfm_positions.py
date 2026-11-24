import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict

EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
TOOL = 'tool=pcft-position-set&email=karagolt@mit.edu'

THREE = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q',
    'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
    'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W',
    'Tyr': 'Y', 'Val': 'V',
}
AA1 = set(THREE.values())

QUERIES = [
    'hereditary folate malabsorption AND (PCFT OR SLC46A1)',
    'SLC46A1 AND mutation AND folate',
    'proton-coupled folate transporter AND hereditary folate malabsorption',
    'PCFT AND folate malabsorption AND (missense OR variant OR mutation)',
]


def get(url):
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read()
        except Exception as e:
            if attempt == 3:
                raise
            time.sleep(2 + 2 * attempt)


def esearch(q):
    url = (EUTILS + 'esearch.fcgi?db=pubmed&retmax=300&retmode=json&'
           + TOOL + '&term=' + urllib.parse.quote(q))
    import json
    return json.loads(get(url))['esearchresult'].get('idlist', [])


def efetch(pmids):
    out = []
    for i in range(0, len(pmids), 100):
        chunk = pmids[i:i + 100]
        url = (EUTILS + 'efetch.fcgi?db=pubmed&retmode=xml&' + TOOL
               + '&id=' + ','.join(chunk))
        out.append(get(url))
        time.sleep(0.4)
    return out


human = ''.join(l.strip() for l in open('seq/Q96NT5.fasta')
                if not l.startswith('>'))
assert len(human) == 459

pmids = []
for q in QUERIES:
    ids = esearch(q)
    print(f'query {q!r} -> {len(ids)} records')
    pmids.extend(ids)
    time.sleep(0.4)
pmids = sorted(set(pmids), key=int)
print(f'\nunion of records searched: {len(pmids)}')

pat1 = re.compile(r'\b([' + ''.join(sorted(AA1)) + r'])(\d{1,3})([' + ''.join(sorted(AA1)) + r'])\b')
pat3 = re.compile(r'\b(?:p\.)?(' + '|'.join(THREE) + r')(\d{1,3})(' + '|'.join(THREE) + r')\b')

kept = defaultdict(lambda: {'wt': None, 'subs': set(), 'pmids': set()})
rejected = defaultdict(set)
scanned = 0

for blob in efetch(pmids):
    root = ET.fromstring(blob)
    for art in root.iter('PubmedArticle'):
        pmid_el = art.find('.//PMID')
        pmid = pmid_el.text if pmid_el is not None else '?'
        parts = []
        for tag in ('ArticleTitle', 'AbstractText'):
            for el in art.iter(tag):
                parts.append(''.join(el.itertext()))
        text = ' '.join(parts)
        if not text:
            continue
        scanned += 1
        hits = [(m.group(1), int(m.group(2)), m.group(3)) for m in pat1.finditer(text)]
        hits += [(THREE[m.group(1)], int(m.group(2)), THREE[m.group(3)])
                 for m in pat3.finditer(text)]
        for wt, pos, mut in hits:
            if not (1 <= pos <= len(human)):
                rejected[(wt, pos, mut)].add(pmid)
                continue
            if human[pos - 1] != wt:
                rejected[(wt, pos, mut)].add(pmid)
                continue
            rec = kept[pos]
            rec['wt'] = wt
            rec['subs'].add(mut)
            rec['pmids'].add(pmid)

print(f'records with text scanned: {scanned}')
print(f'\n=== positions whose wild-type residue matches Q96NT5: {len(kept)} ===')
print('pos  wt  substitutions          n_pmids  pmids')
for pos in sorted(kept):
    r = kept[pos]
    print(f'{pos:4d}  {r["wt"]}   {",".join(sorted(r["subs"])):22s} '
          f'{len(r["pmids"]):7d}  {",".join(sorted(r["pmids"], key=int))}')

print(f'\n=== tokens rejected because the wild-type residue does not match '
      f'Q96NT5 or the position is out of range: {len(rejected)} ===')
for (wt, pos, mut), pm in sorted(rejected.items(), key=lambda kv: kv[0][1]):
    at = human[pos - 1] if 1 <= pos <= len(human) else 'out-of-range'
    print(f'  {wt}{pos}{mut}: Q96NT5 has {at} at {pos}; pmids {",".join(sorted(pm, key=int))}')

with open('data/hfm_positions.tsv', 'w') as fh:
    fh.write('human_pos\twt_aa\tsubstitutions\tpmids\n')
    for pos in sorted(kept):
        r = kept[pos]
        fh.write(f'{pos}\t{r["wt"]}\t{",".join(sorted(r["subs"]))}\t'
                 f'{",".join(sorted(r["pmids"], key=int))}\n')
print('\nwrote hfm_positions.tsv')
