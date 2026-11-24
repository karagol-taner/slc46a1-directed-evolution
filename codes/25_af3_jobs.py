import json
import os
import sys

SEED = 20261113

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ADM = os.path.join(ROOT, 'data', 'admissible_construct.json')
POSSET = os.path.join(ROOT, 'data', 'pharmacophore_position_set.tsv')
SEQF = os.path.join(ROOT, 'seq', '7BC6.fasta')
FASTA = os.path.join(ROOT, 'data', 'construct_sequences.fasta')
AF3 = os.path.join(ROOT, 'data', 'af3')
MANIFEST = os.path.join(ROOT, 'data', 'af3_jobs_manifest.json')
LOG = os.path.join(ROOT, 'data', '25_af3_jobs_log.txt')
AF3_DIALECT, AF3_VERSION = 'alphafold3', 4
SERVER_DIALECT, SERVER_VERSION = 'alphafoldserver', 1
AMINO_ACIDS = set('ACDEFGHIKLMNPQRSTVWY')

lines = []


def log(m):
    print(m, flush=True)
    lines.append(str(m))


def die(m):
    with open(LOG, 'w') as f:
        f.write('\n'.join(lines + ['FATAL ' + m]) + '\n')
    sys.exit('FATAL ' + m)


for p in (ADM, POSSET, SEQF):
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
for p in PRIM:
    if chicken[p - 1] != PSET[p]['chicken_aa']:
        die('the chicken sequence has %s at chicken %d, the position set says %s'
            % (chicken[p - 1], p, PSET[p]['chicken_aa']))

adm = json.load(open(ADM))
if sorted(adm['primaries']) != PRIM:
    die('the construct file lists different primaries from the position set')
RESC = {int(k): v for k, v in adm['rescuer_positions'].items()}

primaries = list(chicken)
for p in PRIM:
    primaries[p - 1] = PSET[p]['human_aa']
primaries = ''.join(primaries)
construct = list(primaries)
for q, y in RESC.items():
    construct[q - 1] = y
construct = ''.join(construct)

if construct != adm['sequence']:
    n = sum(1 for a, b in zip(construct, adm['sequence']) if a != b)
    die('the construct assembled here differs from the stored sequence at %d positions' % n)
d_prim = [i + 1 for i, (a, b) in enumerate(zip(chicken, primaries)) if a != b]
if d_prim != PRIM:
    die('the primaries arm differs from chicken at %s, the %d primary positions expected'
        % (d_prim, len(PRIM)))
d_con = [i + 1 for i, (a, b) in enumerate(zip(primaries, construct)) if a != b]
if d_con != sorted(RESC):
    die('the construct differs from the primaries arm at %s, the rescuer positions %s '
        'expected' % (d_con, sorted(RESC)))

ARMS = [('wt', chicken, 'wild-type chicken transporter, chain A of 7BC6 as deposited'),
        ('primaries', primaries,
         'the %d humanising substitutions alone, no rescuers' % len(PRIM)),
        ('construct', construct,
         'the %d humanising substitutions with the %d rescuers the scans chose'
         % (len(PRIM), len(RESC)))]
for name, seq, _ in ARMS:
    stray = set(seq) - AMINO_ACIDS
    if stray:
        die('%s carries characters that are not standard amino acids: %s'
            % (name, sorted(stray)))
    if len(seq) != len(chicken):
        die('%s is %d residues, the chicken sequence is %d' % (name, len(seq), len(chicken)))
if len({seq for _, seq, _ in ARMS}) != len(ARMS):
    die('two arms carry the same sequence')

with open(FASTA, 'w') as fh:
    for name, seq, note in ARMS:
        fh.write('>%s %s\n' % (name, note))
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + '\n')

for sub in ('server', 'local'):
    os.makedirs(os.path.join(AF3, sub), exist_ok=True)


def server_job(name, seq):
    return {'name': name, 'modelSeeds': [SEED],
            'sequences': [{'proteinChain': {'sequence': seq, 'count': 1}}],
            'dialect': SERVER_DIALECT, 'version': SERVER_VERSION}


def local_job(name, seq):
    return {'name': name, 'modelSeeds': [SEED],
            'sequences': [{'protein': {'id': 'A', 'sequence': seq}}],
            'dialect': AF3_DIALECT, 'version': AF3_VERSION}


written = []
server_path = os.path.join(AF3, 'server', 'slc46a1_af3_jobs.json')
json.dump([server_job(n, s) for n, s, _ in ARMS], open(server_path, 'w'), indent=1)
written.append(server_path)
for name, seq, _ in ARMS:
    path = os.path.join(AF3, 'local', name + '.json')
    json.dump(local_job(name, seq), open(path, 'w'), indent=1)
    written.append(path)

for path in written:
    reloaded = json.load(open(path))
    for job in (reloaded if isinstance(reloaded, list) else [reloaded]):
        if not job.get('sequences'):
            die('%s reloaded with no sequences' % path)
        got = [c.get('protein', c.get('proteinChain', {})).get('sequence')
               for c in job['sequences']]
        want = dict((n, s) for n, s, _ in ARMS)[job['name']]
        if got != [want]:
            die('%s reloaded with a sequence that is not the one intended for %s'
                % (path, job['name']))

readme = os.path.join(AF3, 'README.txt')
with open(readme, 'w') as fh:
    fh.write(
        'AlphaFold3 inputs for the SLC46A1 engraftment construct set.\n\n'
        'Sequences are taken from ../construct_sequences.fasta, which '
        'codes/25_af3_jobs.py derives\nfrom seq/7BC6.fasta and '
        'data/admissible_construct.json. The files here carry the project\nseed %d, '
        'so a run can be repeated.\n\n'
        'server/slc46a1_af3_jobs.json\n'
        '    All %d jobs in one upload for the AlphaFold Server, protein only, in the\n'
        '    proteinChain form, one chain of count 1 for each arm. If the server '
        'declines the\n    seed, replace each modelSeeds value with an empty list and it '
        'will assign one,\n    which the downloaded job request then records.\n\n'
        'local/<job>.json\n'
        '    One file per job for AlphaFold3 run locally, dialect %s, protein only, '
        'chain A.\n\n'
        'The three arms are the wild-type chicken transporter, the %d humanising\n'
        'substitutions alone, and the construct, which is those %d plus the %d '
        'rescuers. The\nmiddle arm is there because without it nothing separates the '
        'rescuers helping from\nthe humanisation being tolerated without them.\n\n'
        'The predictions are not the deposited frame, so codes/26_af3_models_read.py '
        'reports\neach model against 7BC6 and 7BC7 before anything downstream uses it: '
        'a transporter\nhas more than one state and the arms have to be in the same one '
        'to be comparable.\n'
        % (SEED, len(ARMS), AF3_DIALECT, len(PRIM), len(PRIM), len(RESC)))
written.append(readme)
written.append(FASTA)

json.dump({'seed': SEED, 'arms': [n for n, _s, _t in ARMS], 'length': len(chicken),
           'references': ['7BC6', '7BC7'],
           'server_file': os.path.relpath(server_path, ROOT),
           'local_files': [os.path.relpath(os.path.join(AF3, 'local', n + '.json'), ROOT)
                           for n, _s, _t in ARMS],
           'fasta': os.path.relpath(FASTA, ROOT),
           'readme': os.path.relpath(readme, ROOT),
           'reader': 'codes/26_af3_models_read.py',
           'dialects': {'server': [SERVER_DIALECT, SERVER_VERSION],
                        'local': [AF3_DIALECT, AF3_VERSION]},
           'substitutions': {'primaries': PRIM,
                             'rescuers': {str(k): v for k, v in sorted(RESC.items())}}},
          open(MANIFEST, 'w'), indent=1)
written.append(MANIFEST)

log('project seed %d, defined in codes/25_af3_jobs.py' % SEED)
log('sequences derived from %s and %s, checked against the stored construct'
    % (os.path.basename(SEQF), os.path.basename(ADM)))
log('')
log('%-12s %-8s %-28s %s' % ('arm', 'length', 'differs from chicken at', 'note'))
for name, seq, note in ARMS:
    d = [i + 1 for i, (a, b) in enumerate(zip(chicken, seq)) if a != b]
    log('%-12s %-8d %-28s %s' % (name, len(seq), '%d positions' % len(d), note))
log('')
for path in written:
    log('  %-58s %6.1f kB' % (os.path.relpath(path, ROOT), os.path.getsize(path) / 1024))
with open(LOG, 'w') as f:
    f.write('\n'.join(lines) + '\n')
