AlphaFold3 inputs for the SLC46A1 engraftment construct set.

Sequences are taken from ../construct_sequences.fasta, which codes/25_af3_jobs.py derives
from seq/7BC6.fasta and data/admissible_construct.json. The files here carry the project
seed 20261113, so a run can be repeated.

server/slc46a1_af3_jobs.json
    All 3 jobs in one upload for the AlphaFold Server, protein only, in the
    proteinChain form, one chain of count 1 for each arm. If the server declines the
    seed, replace each modelSeeds value with an empty list and it will assign one,
    which the downloaded job request then records.

local/<job>.json
    One file per job for AlphaFold3 run locally, dialect alphafold3, protein only, chain A.

The three arms are the wild-type chicken transporter, the 20 humanising
substitutions alone, and the construct, which is those 20 plus the 6 rescuers. The
middle arm is there because without it nothing separates the rescuers helping from
the humanisation being tolerated without them.

The predictions are not the deposited frame, so codes/26_af3_models_read.py reports
each model against 7BC6 and 7BC7 before anything downstream uses it: a transporter
has more than one state and the arms have to be in the same one to be comparable.
