# Humanising the proton-coupled folate transporter on a chicken scaffold by in silico directed evolution

Taner Karagöl<sup>1,\*,‡</sup>, Alper Karagöl<sup>2,3,\*,‡</sup>

<sup>1</sup> McGovern Institute for Brain Research, Massachusetts Institute of Technology, Cambridge, MA, USA  
<sup>2</sup> Broad Institute of MIT and Harvard, Cambridge, MA, USA  
<sup>3</sup> Department of Mechanical Engineering, Massachusetts Institute of Technology, Cambridge, MA, USA  
<sup>\*</sup> These authors contributed equally.  
<sup>‡</sup> To whom correspondence should be addressed.

## Abstract

The proton-coupled folate transporter is the route by which dietary folates cross the intestinal brush border and by which pemetrexed and related antifolates enter cells in the acidic interior of a tumour. No experimental structure of the human protein exists. The coordinates available are of the chicken orthologue, which is 58.4 per cent identical over 459 aligned positions, so a human pharmacophore can only be studied on a scaffold that lacks it. We substituted the human residue at every position of a 67-position pharmacophore set that we assembled from the structure paper's own assignments, the antifolate contact shell and the hereditary folate malabsorption literature, and we rescued the substitutions the scaffold did not tolerate by computational epistatic rescue. The two species differ at 20 of these positions. At 6.0 Angstrom, 17 of them contacted another, and 13 of the 22 contacting pairs lay beyond the published search window of plus or minus 3 residues, the largest separation being 242 residues. Widening that window had no detectable effect on how many targets were rescued, and ranking candidates by coupling showed no detectable advantage at a matched count. These results suggest that the rest of the set, and not the search, governs the outcome. Scored one at a time on the wild-type background, five high-risk substitutions each admitted a rescue above a likelihood threshold of 0.001; scored with the other 19 present, 2 of 5 survived and the overstatement reached 9.3-fold. We rescored 35,872 single substitutions in the construct background, which placed the median at 0.92 of its wild-type value. Dropping one rescuer raised the 4 humanising substitutions within 20 residues of it while the rest held, which suggests that a rescuer lowers likelihood at the substitutions nearest it. The finished construct contains 20 humanising substitutions and 6 rescuer positions, 26 changes in all, with 19 of them above the threshold and no pharmacophore position overwritten. We then simulated all three arms in a lipid bilayer, 2 replicates of 100 nanoseconds each. All stayed folded, and at chicken 239 the residue fluctuation fell in the same order as the predicted confidence in every replicate, which may indicate that the lower confidence there reflects something other than a loss of order. The numbers are changes in sequence likelihood and not calculated free energies, a ratio between two such likelihoods is not a ratio of free energies, and the humanised surrogate requires experimental validation of folding, binding and transport before use in pharmacological screening.

## Keywords

Proton-coupled folate transporter, SLC46A1, Protein humanisation, Epistatic rescue, Protein language models, Antifolate pharmacology

## Repository layout

```
slc46a1-directed-evolution/
├── codes/
│   ├── 01_residue_map.py ... 36_rescue_geometry.py   analysis pipeline, run in order
│   ├── 30_derived_values.py                      collects every reported value into one file
│   ├── 31_build_figures.py                       figures 1 to 4, 7 and S4
│   ├── 33_structure_figures.py                   figures 5 and 6, drawn from the predicted models
│   ├── 34_supplementary.py                       figures S1 to S3 and their tables
│   ├── 35_dynamics_aware.py                      dynamic and rigid conserved scores
│   ├── SLC46A1_step1 ... step9*.ipynb            ESM-2 scans, run on a GPU runtime
│   ├── SLC46A1_RUNALL.ipynb                      the nine scans end to end
│   ├── colab/                                    molecular dynamics notebooks
│   └── fetch_sources.sh                          fetches and checks the sequences and structures
├── data/
│   ├── esm2_step1/ ... esm2_step9/               likelihood tables from every scan
│   ├── dca_step3/                                coupling diagnostics
│   ├── evcouplings/                              coupling scores for Q96NT5
│   ├── af3/                                      AlphaFold Server models and job requests
│   ├── md/                                       per-residue fluctuation and run reports
│   ├── dynamics_aware/                           dynamic and rigid conserved scores
│   └── *.json *.tsv *_log.txt                    derived tables and run logs
├── structures/                                   the chicken to human residue map
├── figures/                                      figures 1 to 7 and S1 to S4
├── md/                                           protocol, build specification, starting coordinates
├── Supplementary_Information.pdf
├── LICENSE
└── README.md
```

## Dataset summary

| Path | Contents | Files | Size |
| --- | --- | ---: | ---: |
| `data/esm2_step1 to esm2_step9` | Zero-shot likelihood tables from every scan, one row per variant | 46 | 25.7 MB |
| `data/dca_step3` | Coupling diagnostics against the deposited coordinates | 2 | 25 kB |
| `data/evcouplings` | Coupling scores for every residue pair of the human sequence | 1 | 6.1 MB |
| `data/af3` | AlphaFold Server models of the three arms, five each, with job requests | 36 | 4.5 MB |
| `data/md` | Per-residue fluctuation and the report of every production run | 14 | 86 kB |
| `data/dynamics_aware` | Dynamic and rigid conserved scores at every scored position of an arm | 3 | 87 kB |
| `data (top level)` | Derived tables, reports and the run log of every script | 61 | 546 kB |
| `structures` | The chicken to human residue map every correspondence rests on | 1 | 11 kB |
| `figures` | Figures 1 to 7 and S1 to S4 at 400 dpi, with vector copies | 22 | 7.8 MB |
| `md` | Simulation protocol, build specification and starting coordinates | 15 | 2.7 MB |
| whole repository | | 253 | 49.7 MB |

## Reproducing

The pipeline runs in four passes. Sequences and structures are fetched rather
than redistributed, so the first step is required before any other.

1. `bash codes/fetch_sources.sh` writes `seq/` with the two sequences and the two
   deposited structures, checks the length of each and prints a checksum for the
   coordinate files.
2. The nine notebooks `codes/SLC46A1_step1 ... step9` run the zero-shot scans on a
   GPU runtime and write `data/esm2_step1/ ... esm2_step9/`. `SLC46A1_RUNALL.ipynb`
   runs all nine in one session. Every scan writes a partial table as it goes, so an
   interrupted run resumes from the checkpoint instead of starting again.
3. `python3 codes/01_residue_map.py` through `codes/36_rescue_geometry.py`, in numeric
   order, read those tables together with the structures and write the reports in
   `data/`. Each script writes its own log beside its output and stops on a missing
   input rather than carrying a default forward.
4. `python3 codes/30_derived_values.py` collects every value the analysis reports
   into `data/derived_values.json`, and `31_build_figures.py`,
   `33_structure_figures.py` and `34_supplementary.py` draw the figures from it.

Structure prediction and simulation sit outside that sequence. The three AlphaFold
jobs are in `data/af3/server/` for the AlphaFold Server and in `data/af3/local/`
for a local AlphaFold3 install, both carrying the project seed. The bilayer systems
are built with CHARMM-GUI following `md/inputs/CHARMM-GUI-checklist.txt`, and the
notebooks in `codes/colab/` run the simulations and the comparison.

## Software dependencies

| Software | Version | Used for |
| --- | --- | --- |
| Python | 3.10.12 | Every script in codes/ |
| PyTorch | 2.14.1+cu130 | ESM-2 inference in the step notebooks, on a GPU |
| fair-esm | not recorded in the run logs | ESM-2 t33 650M weights and tokeniser |
| NumPy | 2.2.6 | Array arithmetic throughout |
| gemmi | 0.7.5 | Reading mmCIF coordinates and computing contacts |
| matplotlib | 3.10.9 | Figures 1 to 7 and S1 to S4 |
| Biopython | not recorded in the run logs | Pairwise alignment in 01 and 02 |
| EVcouplings | public run Q96NT5_b0.3 | Direct coupling analysis |
| AlphaFold Server | AlphaFold3 | Structure prediction of the three arms |
| CHARMM-GUI | Membrane Builder | Bilayer systems for the simulations |
| GROMACS | as provided by the Colab runtime | Simulation, through codes/colab |

## Data availability

The human sequence is UniProt Q96NT5 and the chicken sequence is UniProt
F1NJ67. The coordinates are PDB 7BC6, the substrate-free chicken transporter, and
7BC7, its pemetrexed complex. None of these are redistributed here;
`codes/fetch_sources.sh` downloads them from UniProt and the RCSB and checks what it
received before anything reads it.

The couplings come from a public EVcouplings run against UniProt, job `Q96NT5_b0.3`.
`data/evcouplings/` carries the coupling scores that run produced. The alignment
itself is not redistributed and is regenerated by rerunning EVcouplings on Q96NT5 at
a bitscore threshold of 0.3.

`data/af3/` contains AlphaFold Server Output. That Output and anything substantially
derived from it are provided under and subject to the AlphaFold Server Output Terms
of Use, reproduced in `data/af3/terms_of_use.md` and published at
https://alphafoldserver.com/output-terms. The models here are unmodified Server
output for the three sequences in `data/construct_sequences.fasta`.

## Citation

Karagöl T, Karagöl A. Humanising the proton-coupled folate transporter on a chicken scaffold by in silico directed evolution. Manuscript in preparation, 2026.

The Supplementary Information that accompanies the manuscript is in this repository
as `Supplementary_Information.pdf`.

## Contact

Taner Karagöl, karagolt@mit.edu  
Alper Karagöl, akaragol@mit.edu

## License

The code in this repository is released under the MIT License, in `LICENSE`. The upstream data keep their own licences: UniProt and the RCSB for the sequences and coordinates, EVcouplings for the couplings, and the AlphaFold Server Output Terms of Use for everything in `data/af3/`.
