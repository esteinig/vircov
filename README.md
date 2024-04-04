# vircov <a href='https://github.com/esteinig'><img src='docs/vircov.png' align="right" height="270"/></a>

[![build](https://github.com/esteinig/nanoq/actions/workflows/rust-ci.yaml/badge.svg?branch=master)](https://github.com/esteinig/nanoq/actions/workflows/rust-ci.yaml)
[![codecov](https://codecov.io/gh/esteinig/vircov/branch/main/graph/badge.svg?token=RG95F4C6FE)](https://codecov.io/gh/esteinig/vircov)
![](https://img.shields.io/badge/version-0.6.0-black.svg)

Minimal virus genome coverage assessment for metagenomic diagnostics

## Overview


**`v0.6.0`**

- [Purpose](#purpose)
- [Implementation](#implementation)
- [Installation](#installation)
- [Usage](#usage)
- [Tests](#tests)
- [Concept](#concept)
- [Clinical examples](#clinical-examples)
- [Performance](#performance)
- [Etymology](#concept)
- [Contributors](#contributors)

## Purpose

Viral metagenomic diagnostics from low-abundance clinical samples can be challenging in the absence of sufficient genome coverage. `Vircov` extracts distinct non-overlapping regions from a reference alignment and generates some helpful statistics. It can be used to flag potential hits without inspection of coverage plots in automated pipelines and reports.

## Implementation

`Vircov` is written in Rust and works with alignments in the standard `PAF` or `SAM/BAM/CRAM` (next release) formats. It is extremely fast and can process alignments against thousands of viral reference genomes in seconds. Basic input filters can be selected to remove spurious alignments and text-style coverage plots can be printed to the terminal for visual confirmation.

`Vircov` is written for implementation in metagenomics pipelines for human patients enroled in the `META-GP` network (Australia). As such, it attempts to be production-grade code with high test coverage, continuous integration, and versioned releases with precompiled binaries for Linux and MacOS.

Version `0.6.0` will be distributed on `BioConda` and `Cargo` and will have precompiled binaries available.

## Installation

```bash
git clone https://github.com/esteinig/vircov 
cd vircov && cargo build --release
```

## Tests

```bash
git clone https://github.com/esteinig/vircov 
cd vircov && cargo test && cargo tarpaulin 
```

## Concept

Definitive viral diagnosis from metagenomic clinical samples can be extremely challenging due to low sequencing depth, large amounts of host reads and low infectious titres, especially in blood or CSF. One way to distinguish a positive viral diagnosis is to look at alignment coverage against one or multiple reference sequences. When only few reads map to the reference, and genome coverage is low, positive infections often display multiple distinct alignment regions, as opposed to reads mapping to a single or few regions on the reference. This has been implemented for example in `SURPI+` used by Miller et al. (2019) for DNA and RNA viruses in CSF samples. [De Vries et al. (2021)](https://www.sciencedirect.com/science/article/pii/S1386653221000792) summarize this concept succinctly in this figure (adapted):

![devries](https://user-images.githubusercontent.com/12873366/158775480-447d847e-5b0d-487c-a39a-81bdf428e09d.png)

Positive calls in these cases can be made from coverage plots showing the distinct alignment regions and a threshold on the number of regions is chosen by the authors (> 3). However, coverage plots require visual assessment and may not be suitable for flagging potential hits in automated pipelines or summary reports. 

`Vircov` attempts to make visual inspection and automated flagging easier by counting the distinct (non-overlapping) coverage regions in an alignment and reports some helpful statistics to make an educated call based on coverage information. We have specifically implemented grouping options by for example `taxid` in the reference fasta headers and account for segmented viruses and those split into genes in some databases like `Virosaurus`.

In addition the most recent version implements single reference selection based on first grouping the reference sequences that have significant alignments by `taxid` or species name, and then selecting the reference equence with the highest number of unique mapped reads (for example). This allows for using `Vircov` with permissive (no filter) settings to select single reference genoems for re-mapping and and provides sensitive coverage information. 

## Usage examples

```bash
# Output coverage statistics from a SAM|BAM|PAF alignment
vircov --alignment test.paf > stats.tsv

# Output coverage statistics with headers of reference sequences
# used in the alignment, including non aligned references
vircov --alignment test.paf --fasta ref.fa --zero -v > stats.zero.tsv

# Output coverage statistics with threshold filters activated
vircov --alignment test.bam \
   --fasta ref.fa \
   --min-len 50 \                  # minimum query alignment length
   --min-cov 0 \                   # minimum query alignment coverage
   --min-mapq 50 \                 # minimum mapping quality
   --reads 3 \                     # minimum reads aligned against a reference
   --coverage 0.05 \               # minimum reference coverage fraction
   --regions 4 \                   # minimum distinct alignment regions
   --regions-coverage 0.3 \        # distinct regions filter applies < coverage of 30%
   --read-ids ids.txt \            # read identifiers of mapped reads of surviving alignments
   -v > stats.tsv

# Group alignments by field in reference description and select best coverage reference
# for each group - header field example: "name=virus; taxid=12345; segment=N/A". Use the
# "segment=" field to select a single best segment if aligned against and output selected
# segments in one sequence file (e.g. for Influenza genomes). Each selected reference
# or selected segments are output as sequence files into "outref".

vircov --alignment test.paf  \
   --fasta ref.fa \
   --min-len 50 \                      # minimum query alignment length
   --min-cov 0 \                       # minimum query alignment coverage
   --min-mapq 50 \                     # minimum mapping quality
   --reads 3 \                         # minimum reads aligned against a reference
   --coverage 0.05 \                   # minimum reference coverage fraction
   --regions 4 \                       # minimum distinct alignment regions
   --regions-coverage 0.3 \            # distinct regions filter applies < coverage of 30%
   --group-by "taxid=" \               # field to group alignments by
   --group-sep ";" \                   # field separator in reference description
   --group-select-by "cov" \           # select best coverage reference sequences
   --group-select-split outref \       # output selected references into this folder
   --group-select-order \              # order by highest coverage, include as index in filename 
   --segment-field "segment=" \        # group and select segmented genomes with field
   --segment-field-nan "N/A" \         # no segment value in segment field
   -v > stats.selected.tsv             # output selected references alignment stats
```

## Performance

Alignments conducted with `minimap2 -c -x sr` (PAF) and `minimap2 --sam-hits-only -ax sr` (SAM/BAM, only aligned sequences). Peak memory is mainly determined by aligned interval records that are stored for the overlap computations. It may vary depending on how many alignments remain after filtering, the number of aligned reads, output formats and size of the reference database.

* **Sample 1**: ~ 70 million Illumina PE reads against ~ 70k reference genomes, ~ 1.7 million alignments 
* **Sample 2**: ~ 80 million Illumina PE reads against ~ 70k reference genomes, ~ 63 million alignments 

`.paf`
    
  * **Sample 1** (212 MB): 0.56 seconds, peak memory: 12 MB 
  * **Sample 2** (8.9 GB): 55.4 seconds, peak memory: 2.9 GB


## Etymology

Not a very creative abbreviation of "virus coverage" but the little spectacles in the logo are a reference to [Rudolf Virchow](https://en.wikipedia.org/wiki/Rudolf_Virchow) who described such trivial concepts as cells, cancer and pathology. His surname is pronounced like `vircov` if you mumble the terminal `v`.

## Contributors

* Prof. Deborah Williamson and Prof. Lachlan Coin (principal investigators for `META-GP`)
* Dr. Leon Caly (samples and sequencing for testing on clinical data)

