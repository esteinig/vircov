# vircov <a href='https://github.com/esteinig'><img src='docs/vircov.png' align="right" height="200"/></a>

![](https://img.shields.io/badge/version-0.8.0-black.svg)

Viral genome coverage metrics from read alignments and genomic neighbor typing schemes for rapid subtyping, clade-assignment and genotype inference from assembled genomes.

## Overview


**`v0.8.0`**

- [Purpose](#purpose)
- [Implementation](#implementation)
- [Example applications](#example-applications)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Concepts](#concepts)
- [Usage examples](#usage-examples)
- [Etymology](#etymology)
- [Contributors](#contributors)

## Purpose

`Vircov` implements two modules that address common challenges associated with identification and recovery of viral genomes for metagenomic diagnostics applications:

* Viral metagenomic diagnostics from low-abundance clinical samples can be challenging in the absence of sufficient genome coverage. `Vircov` extracts distinct non-overlapping regions from a reference alignment and generates some coverage statistics. It can be used to flag potential hits without inspection of coverage plots in automated pipelines and reports. Coverage evaluations and automated selection of reference genomes for downstream consensus assembly form the initial step in detection of viral genomes from the "scan-remap" virus-focused pipeline in `Cerebro`.

* Viral assembly subtyping can be challenging even when full genomes can be recovered from enriched metagenomic data due to differences in genome-informed lineage subtyping schemes such as those inferred from phylogenies with Nextstrain and other traditional schemes e.g. based on serology and single gene typing for a range of highly divergent viruses. `Vircov` integrates custom, user-defined and reference-database derived subtyping schemes whose construction can be automated with `Cipher`. We make automatically updated [subtyping databases]() available for a range of common viruses where data-sharing arrangements make this possible (e.g. NCBI- but not GISAID-derived assemblies). Schemes are checked for inconsistiencies in automated annotation extraction by our bioinformatics team at the Victorian Infectious Diseases Reference Laboratory (VIDRL) in Melbourne.

* `Vircov` rapidly computes average amino acid and nucleotide identities (AAI and ANI) as well as mutual nearest neighbor population graphs ([Neuditschko et al. 2012](), [Steinig et al. 2016]()) based on genome similarity or phylogenetic distances to rank and assign subtypes, clades and other genome-associated meta-data from consensus assemblies - a form of genomic neighbor typing previously applied to bacterial genomes and AMR inference ([Brinda et al. 2020](https://www.nature.com/articles/s41564-019-0656-6), [Steinig et al. 2022](https://www.biorxiv.org/content/10.1101/2022.02.05.479210v1.full)).  



## Implementation

`Vircov` is written in Rust and primarily operates on alignments in the standard `PAF` or `SAM/BAM/CRAM` formats and consensus genome assemblies in `FASTA` format. It is fast and can process alignments against thousands of viral reference genomes and compute subtypes for genome assemblies in seconds. Basic filters can remove spurious alignments or subtype inferences. ASCII-style coverage plots can be printed to the terminal for visual checks. Mutual nearest neighbor graphs can optionally be visualized along with relevant metadata, for example using graph decorators syntax and `igraph` plots with `NetView`.

`Vircov` is written for implementation in metagenomic diagnostic pipelines for human patients enroled in the `META-GP` network (Australia).

## Example applications

For a walthrough and detailed usage examples please see the [documentation]() for `Vircov`. Some examples of its application to various metagenomic protocols and viral pathogens are published.

Respiratory pathogens from rapid antigen tests using pan-viral enrichment:

> Moso et al. (2024)

> Beutel-Simoes et al. (2024)

Hepatitis E from clinical samples using whole genome primer-scheme:

> O'Keefe et al. (2024)

Central nervous system infections from clinical samples using short-reads:

> Ramachandran et al. (2024)

Enterovirus D68 from clinical sample using pan-viral enrichment:

> Beutel-Simoes et al. (2024)

Dengue from clinical sample using whole genome primer-scheme:

> Beutel-Simoes et al. (2024)



## Dependencies




## Installation

From binaries:

```bash
git clone https://github.com/esteinig/vircov 
cd vircov && cargo build --release
```

From source:

```bash
git clone https://github.com/esteinig/vircov 
cd vircov && cargo build --release
```

ANI and AAI computation using `Vircov` from binary executable or compiled source code require additional dependencies on `$PATH`:

* `BLAST` (`blastn`) 
* `DIAMOND` (`blastx`)


## Usage examples

### Coverage assessment and selection

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


### Genomic neighbor subtyping


## Concepts

### Low-abundance infections and coverage assessment for detection from reads

Definitive viral diagnosis from metagenomic clinical samples can be extremely challenging due to low sequencing depth, large amounts of host reads and low infectious titres, especially in blood or CSF. One way to distinguish a positive viral diagnosis is to look at alignment coverage against one or multiple reference sequences. When only few reads map to the reference, and genome coverage is low, positive infections often display multiple distinct alignment regions, as opposed to reads mapping to a single or few regions on the reference. This has been implemented for example in `SURPI+` used by Miller et al. (2019) for DNA and RNA viruses in CSF samples. [De Vries et al. (2021)](https://www.sciencedirect.com/science/article/pii/S1386653221000792) summarize this concept succinctly in this figure (adapted):

![devries](https://user-images.githubusercontent.com/12873366/158775480-447d847e-5b0d-487c-a39a-81bdf428e09d.png)

Positive calls in these cases can be made from coverage plots showing the distinct alignment regions and a threshold on the number of regions is chosen by the authors (> 3). However, coverage plots require visual assessment and may not be suitable for flagging potential hits in automated pipelines or summary reports. 

`Vircov` attempts to make visual inspection and automated flagging easier by counting the distinct (non-overlapping) coverage regions in an alignment and reports some helpful statistics to make an educated call based on coverage information. We have specifically implemented grouping options by for example `taxid` in the reference fasta headers and account for segmented viruses and those split into genes in some databases like `Virosaurus`.

In addition the most recent version implements single reference selection based on first grouping the reference sequences that have significant alignments by `taxid` or species name, and then selecting the reference equence with the highest number of unique mapped reads (for example). This allows for using `Vircov` with permissive (no filter) settings to select single reference genoems for re-mapping and and provides sensitive coverage information. 

### Genomic neighbor subtyping schemes and genotype inference from assemblies

TBD

## Etymology

Not a very creative abbreviation of "virus coverage" but the little spectacles in the logo are a reference to [Rudolf Virchow](https://en.wikipedia.org/wiki/Rudolf_Virchow) who described such trivial concepts as cells, cancer and pathology. His surname is pronounced like `vircov` if you mumble the terminal `v`.

## Contributors

* Prof. Deborah Williamson and Prof. Lachlan Coin (principal investigators for `META-GP`)
* Dr. Leon Caly (samples and sequencing for testing on clinical data)
* Dr. Ammar Aziz (bioinformatics support and testing)

