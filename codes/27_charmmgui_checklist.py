import collections
import json
import os
import sys

import gemmi

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CS = os.path.join(ROOT, 'md', 'construct_systems.json')
NB_CUT = 3.0
PLAN = os.path.join(ROOT, 'md', 'md_plan.json')
ADM = os.path.join(ROOT, 'data', 'admissible_construct.json')
POSSET = os.path.join(ROOT, 'data', 'pharmacophore_position_set.tsv')
OUTDIR = os.path.join(ROOT, 'md', 'inputs')
OUT = os.path.join(OUTDIR, 'CHARMM-GUI-checklist.txt')
MEASURED = os.path.join(OUTDIR, 'build_measurements.json')
LOG = os.path.join(ROOT, 'data', '27_charmmgui_checklist_log.txt')
SS_CUT = 2.5
POS, NEG = {'ARG', 'LYS'}, {'ASP', 'GLU'}
TITRATABLE = ('HIS', 'ASP', 'GLU', 'LYS', 'ARG', 'CYS', 'TYR')

lines = []


def log(m):
    print(m, flush=True)
    lines.append(str(m))


def die(m):
    with open(LOG, 'w') as f:
        f.write('\n'.join(lines + ['FATAL ' + m]) + '\n')
    sys.exit('FATAL ' + m)


for p in (CS, PLAN, ADM, POSSET):
    if not os.path.exists(p):
        die('missing input %s' % p)
cs = json.load(open(CS))
plan = json.load(open(PLAN))
adm = json.load(open(ADM))
BU = plan.get('build')
for key in ('route', 'force_field', 'lipid', 'salt', 'water_thickness_angstrom',
            'lipids_per_leaflet_target', 'lipid_ratio',
            'lipid_counts_per_leaflet', 'hydrogen_mass_repartitioning',
            'system_basename', 'production_mdp', 'drive_builds_root',
            'template_requires'):
    if key not in BU:
        die('md/md_plan.json records no %s, which the checklist has to state' % key)

rows = [l.rstrip('\n').split('\t') for l in open(POSSET)]
pc = {k: i for i, k in enumerate(rows[0])}
PSET = {int(r[pc['chicken_pos']]): r[pc['human_aa']] for r in rows[1:]}

arms = sorted(cs['arms'])
m = {}
for arm in arms:
    path = os.path.join(ROOT, 'md', 'systems', cs['arms'][arm]['file'])
    if not os.path.exists(path):
        die('%s is named in %s but not on disk' % (cs['arms'][arm]['file'], CS))
    st = gemmi.read_structure(path)
    if len(st[0]) != 1:
        die('%s holds %d chains, the builds are single chain' % (arm, len(st[0])))
    ch = st[0][0]
    nums = sorted(r.seqid.num for r in ch)
    counts = collections.Counter(r.name for r in ch)
    cys = [r for r in ch if r.name == 'CYS']
    ss = []
    for i in range(len(cys)):
        for j in range(i + 1, len(cys)):
            a, b = cys[i].find_atom('SG', '*'), cys[j].find_atom('SG', '*')
            if a is None or b is None:
                continue
            d = a.pos.dist(b.pos)
            if d <= SS_CUT:
                ss.append({'pair': [cys[i].seqid.num, cys[j].seqid.num], 'distance': d})
    m[arm] = {'file': cs['arms'][arm]['file'], 'chain': ch.name,
              'first': nums[0], 'last': nums[-1], 'residues': len(nums),
              'formal_charge': sum(counts[x] for x in POS) - sum(counts[x] for x in NEG),
              'n_cysteine': len(cys),
              'cysteines': sorted(r.seqid.num for r in cys),
              'disulfides': sorted(ss, key=lambda x: x['pair']),
              'titratable': {t: counts[t] for t in TITRATABLE}}

ranges = {(m[a]['first'], m[a]['last']) for a in arms}
if len(ranges) != 1:
    die('the arms span different residue ranges, %s, so one checklist cannot cover them'
        % sorted(ranges))
FIRST, LAST = ranges.pop()
ss_sets = {tuple(tuple(d['pair']) for d in m[a]['disulfides']) for a in arms}
if len(ss_sets) != 1:
    die('the arms carry different disulfides, %s' % sorted(ss_sets))
SS = m[arms[0]]['disulfides']
charges = {a: m[a]['formal_charge'] for a in arms}
diff_tit = {t for t in TITRATABLE if len({m[a]['titratable'][t] for a in arms}) > 1}
SUBS = {p: PSET[p] for p in adm['primaries']}
SUBS.update({int(k): v for k, v in adm['rescuer_positions'].items()})
tit_subs = sorted(p for p, y in SUBS.items()
                  if y in {'H', 'D', 'E', 'K', 'R', 'C', 'Y'})

BB = {'N', 'C', 'CA', 'O'}
SSP = {tuple(d['pair']) for d in SS} | {tuple(reversed(d['pair'])) for d in SS}
NB_MIN, NB_PAIR = 999.0, None
for a in arms:
    st = gemmi.read_structure(os.path.join(ROOT, 'md', 'systems', cs['arms'][a]['file']))
    st.setup_entities()
    ns = gemmi.NeighborSearch(st, 5.0).populate()
    for r in st[0][0]:
        for at in r:
            for mark in ns.find_atoms(at.pos, '\0', radius=NB_CUT):
                cra = mark.to_cra(st[0])
                n2 = cra.residue.seqid.num
                if n2 == r.seqid.num:
                    continue
                if abs(n2 - r.seqid.num) == 1 and {at.name, cra.atom.name} <= BB:
                    continue
                if (r.seqid.num, n2) in SSP and {at.name, cra.atom.name} == {'SG'}:
                    continue
                d = at.pos.dist(cra.atom.pos)
                if d < NB_MIN:
                    NB_MIN, NB_PAIR = d, ('%s%d %s to %s%d %s'
                                          % (r.name, r.seqid.num, at.name,
                                             cra.residue.name, n2, cra.atom.name))
if NB_PAIR is None:
    die('no non-bonded heavy-atom contact within %.1f angstrom was found in any arm, '
        'which cannot be right for a folded protein' % NB_CUT)

os.makedirs(OUTDIR, exist_ok=True)
w = []
w.append('CHARMM-GUI settings for the SLC46A1 engraftment construct set')
w.append('=' * 61)
w.append('')
w.append('Input Generator > Membrane Builder > Bilayer Builder > Upload PDB.')
w.append('Run it %d times, once per file, from md/systems:' % len(arms))
w.append('')
for a in arms:
    w.append('  %-28s %s' % (m[a]['file'],
                             'no substitutions' if not cs['arms'][a]['n_substitutions']
                             else '%d substitutions'
                                  % cs['arms'][a]['n_substitutions']))
w.append('')
w.append('Chains            keep chain %s, the only chain in each file. No ligand is'
         % m[arms[0]]['chain'])
w.append('                  present, since this is the substrate-free state.')
w.append('Disulfides        %s' % (
    'tick the %d pair%s CHARMM-GUI offers: %s. Measured SG to SG at %s'
    % (len(SS), '' if len(SS) == 1 else 's',
       ' and '.join('%d with %d' % tuple(d['pair']) for d in SS),
       ' and '.join('%.2f angstrom' % d['distance'] for d in SS))
    if SS else 'none within %.1f angstrom SG to SG, so nothing is to be ticked'
               % SS_CUT))
if SS:
    w.append('                  in all %d files, so the same pair%s in every build. The '
             % (len(arms), ' appears' if len(SS) == 1 else 's appear'))
    w.append('                  %d cysteines are at %s and no other pair is close.'
             % (m[arms[0]]['n_cysteine'],
                ', '.join(str(x) for x in m[arms[0]]['cysteines'])))
w.append('Missing residues  none inside the chain. The %d residues span %d to %d with no'
         % (m[arms[0]]['residues'], FIRST, LAST))
w.append('                  internal gap, so there is nothing for CHARMM-GUI to model.')
w.append('Terminal groups   neutral caps, ACE and CT3 or their equivalents. Residue %d is'
         % FIRST)
w.append('                  not the start of the protein and residue %d is not its end:'
         % LAST)
w.append('                  both are where the deposited model stops, so charged NTER and')
w.append('                  CTER would put formal charges in the middle of the sequence.')
w.append('Protonation       defaults at the pH the state carries. The arms are not')
w.append('                  identical in titratable content: %s'
         % ('%s differ across the arms'
            % ', '.join(t.title() for t in sorted(diff_tit))
            if diff_tit else 'they match on every titratable type'))
if 'HIS' in diff_tit:
    w.append('                  and one humanising substitution installs a histidine, so')
    w.append('                  the wild type carries %d and the others %d.'
             % (m['wt']['titratable']['HIS'], m['construct']['titratable']['HIS']))
w.append('Lipid             %s, %d per leaflet as the target, entered as a ratio.'
         % (BU['lipid'], BU['lipids_per_leaflet_target']))
if BU.get('lipid_counts_per_leaflet'):
    w.append('                  The ratio resolves to %s per leaflet.'
             % ', '.join('%d %s' % (v, k)
                         for k, v in BU['lipid_counts_per_leaflet'].items()))
w.append('Water             %.1f angstrom thickness.' % BU['water_thickness_angstrom'])
w.append('Ions              %s, neutralising. The arms do not carry the same formal'
         % BU['salt'])
w.append('                  charge, %s, so the counter-ion counts will differ by %d'
         % (', '.join('%s %+d' % (a, charges[a]) for a in arms),
            max(charges.values()) - min(charges.values())))
w.append('                  between the extremes. A build whose charge disagrees with the')
w.append('                  figure above means the wrong file was uploaded.')
w.append('Force field       %s.' % BU['force_field'])
w.append('Dynamics          GROMACS.')
w.append('')
w.append('')
w.append('Step by step')
w.append('-' * 12)
w.append('')
w.append('1.  Open charmm-gui.org, log in, and choose Input Generator >')
w.append('    Membrane Builder > Bilayer Builder. The protein is a transporter in a')
w.append('    bilayer, so this is the builder, not Solution Builder.')
w.append('2.  Choose "Upload PDB" and upload the first file from md/systems. The three')
w.append('    files, in the order they are listed above, are:')
for a in arms:
    w.append('      %s' % m[a]['file'])
w.append('    Do one at a time and keep the three downloads apart. The run is the same')
w.append('    for each, and only the uploaded file changes.')
w.append('3.  PDB Info page. Keep chain %s and nothing else. Leave the het records'
         % m[arms[0]]['chain'])
w.append('    alone, there are none to keep.')
w.append('4.  On the same page, the disulfide and protonation panel. %s'
         % ('Tick the pair %s. Leave the protonation at its defaults.'
            % ' and '.join('%d with %d' % tuple(d['pair']) for d in SS)
            if SS else 'There is no pair to tick. Leave the protonation at its defaults.'))
w.append('    Read off which histidine form it assigns at the %d titratable substituted'
         % len(tit_subs))
w.append('    positions, %s, before going on.' % ', '.join(str(x) for x in tit_subs))
w.append('5.  Terminal group patching. Set both ends to neutral caps rather than')
w.append('    charged NTER and CTER, because residue %d is not the start of the'
         % FIRST)
w.append('    protein and residue %d is not its end.' % LAST)
w.append('6.  Orientation. Choose "Run PPM 2.0" and leave every Positioning Option')
w.append('    unticked. PPM sets the membrane normal along z itself, so a rotation or a')
w.append('    translation entered by hand afterwards undoes what it just did. The')
w.append('    computed energy printed above those options is a CHARMM single point on')
w.append('    the rebuilt coordinates: a very large van der Waals term there comes from')
w.append('    the side chains and hydrogens CHARMM-GUI has just added, since in the')
w.append('    uploaded files the closest non-bonded heavy-atom contact is %.2f angstrom'
         % NB_MIN)
w.append('    (%s), and the only pair nearer than that is the %s disulfide at %.2f.'
         % (NB_PAIR, ' with '.join(str(x) for x in SS[0]['pair']) if SS else 'none',
            SS[0]['distance'] if SS else 0.0))
w.append('    It falls away during minimisation and equilibration, so read it after')
w.append('    those rather than here.')

w.append('7.  Lipid page. %s, %d per leaflet as the target, symmetric leaflets.'
         % (BU['lipid'], BU['lipids_per_leaflet_target']))
if BU.get('lipid_ratio'):
    w.append('    Enter the ratio, %s, and leave the counts to the page. It'
             % ' to '.join(str(v) for v in BU['lipid_ratio'].values()))
    w.append('    resolves to %s per leaflet. The XY box length the page asks'
             % ', '.join('%d %s' % (v, k)
                         for k, v in BU['lipid_counts_per_leaflet'].items()))
    w.append('    for in ratio mode is an initial guess and not a value to leave')
    w.append('    blank, since a blank there selects no lipid and the build stops.')
w.append('8.  Water and ions. %.1f angstrom water thickness, %s, neutralising. Check'
         % (BU['water_thickness_angstrom'], BU['salt']))
w.append('    the net charge it reports against the figure above for that arm.')
w.append('9.  Force field and equilibration input. %s, and tick GROMACS as the'
         % BU['force_field'])
_T = plan['run_budget'].get('temperature_k')
w.append('    dynamics engine. Leave hydrogen mass repartitioning unticked.')
w.append('    Temperature: %s'
         % ('%g K.' % _T if _T is not None else
            '>>> NOT ON RECORD IN md/md_plan.json, CHECK BEFORE BUILDING <<<'))
w.append('10. Download the charmm-gui-*.tgz when the build finishes, then repeat from')
w.append('    step 2 for the next file. Three builds in all.')
w.append('11. Unpack each archive. Inside is a gromacs/ folder holding')
import textwrap as _tw
for _l in _tw.wrap(', '.join(BU['template_requires']) + '.', 70):
    w.append('    ' + _l)
w.append('    The contents of that folder are what go on Drive, without the')
w.append('    gromacs/ level around them.')
w.append('')
w.append('Unpack each charmm-gui-*.tgz and put the contents of its gromacs/ folder')
w.append('on Drive as')
w.append('')
for a in arms:
    w.append('  %s/%s' % (BU['drive_builds_root'], a))
w.append('')
w.append('A build folder is recognised by holding %s.pdb, so the folder names'
         % BU['system_basename'])
w.append('above matter and the names CHARMM-GUI gives the archives do not.')
w.append('')
w.append('The %d titratable positions among the %d substituted ones are %s, and the'
         % (len(tit_subs), len(SUBS), ', '.join(str(x) for x in tit_subs)))
w.append('protonation CHARMM-GUI assigns them is worth reading off step1_pdbreader.out')
w.append('before the runs start, since a changed protonation state at one of those would')
w.append('be a difference between arms that no substitution put there.')
with open(OUT, 'w') as fh:
    fh.write('\n'.join(w) + '\n')
json.dump({'arms': m, 'residue_range': [FIRST, LAST], 'disulfides': SS,
           'formal_charges': charges, 'charge_spread': max(charges.values())
           - min(charges.values()),
           'titratable_differing': sorted(diff_tit),
           'titratable_substituted_positions': tit_subs,
           'closest_nonbonded_contact': {'distance': NB_MIN, 'pair': NB_PAIR,
                                         'excludes_disulfide': True,
                                         'search_radius': NB_CUT},
           'drive_paths': ['%s/%s' % (BU['drive_builds_root'], a)
                           for a in arms],
           'settings_from': os.path.relpath(PLAN, ROOT)},
          open(MEASURED, 'w'), indent=1)

log('%d arms, chain %s, residues %d to %d, %d residues'
    % (len(arms), m[arms[0]]['chain'], FIRST, LAST, m[arms[0]]['residues']))
log('disulfides within %.1f angstrom: %s'
    % (SS_CUT, ', '.join('%d with %d at %.2f' % (d['pair'][0], d['pair'][1],
                                                 d['distance']) for d in SS) or 'none'))
log('formal charges: %s, spread %d'
    % (', '.join('%s %+d' % (a, charges[a]) for a in arms),
       max(charges.values()) - min(charges.values())))
log('titratable types differing across arms: %s' % (sorted(diff_tit) or 'none'))
log('titratable substituted positions: %s' % tit_subs)
log('')
log('wrote %s (%d bytes)' % (os.path.relpath(OUT, ROOT), os.path.getsize(OUT)))
log('wrote ' + os.path.relpath(MEASURED, ROOT))
with open(LOG, 'w') as f:
    f.write('\n'.join(lines) + '\n')
