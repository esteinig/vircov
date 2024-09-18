# vircov <a href='https://github.com/esteinig'><img src='docs/vircov.png' align="right" height="270"/></a>

![](https://img.shields.io/badge/version-0.8.0-black.svg)

Automated coverage statistics, genome recovery and subtyping for metagenomic diagnostics of viral infections from reads or alignments.

## Overview


**`v1.0.0`**

- [Purpose](#purpose)
- [Install](#install)
- [Examples](#examples)
- [Concepts](#concepts)
- [Etymology](#etymology)

## Purpose

`Vircov` implments two modules that address common challenges associated with identification and recovery of viral genomes for metagenomic diagnostics applications:

1. Viral metagenomic diagnostics from low-abundance clinical samples can be challenging in the absence of sufficient genome coverage. `Vircov` extracts distinct non-overlapping regions from a reference alignment and generates some helpful coverage statistics. It can be used to flag potential hits without inspection of coverage plots in automated pipelines and reports. Coverage evaluations and automated selection of reference genomes for subsequent consensus assembly form the initial step in detection of viral genomes from the scan-remap pipeline in `Cerebro`.

2. Viral subtyping can be challenging in the absence of consistent subtyping schemes or where detailed tracking of lineage emergence - such as for SARS-CoV-2 or other viruses covered by Nextstrain - is not desired. `Vircov` integrates NCBI Virus derived subtyping schemes and reference database construction automated with `Cipher`. It rapidly compute average amino acid and nucleotide identities (AAI and ANI) as well as mutual nearest neighbor population graphs based on genome similarity or phylogenetic trees to infer genotypes from consensus assemblies. 

## Install

Anaconda installation with dependencies:

```
mamba create -n vircov -c esteinig vircov
```

Dependencies for full pipeline from reads:

* `bowtie2` | `minimap2` | `strobealign` 
* `samtools`
* `ivar` 

Source installation without dependencies:

```
git clone https://github.com/esteinig/vircov
cd vircov && cargo build --release 
```

## Usage

### Virus detection and whole genome recovery


### Virus assembly subtyping


## Concepts

### Low-abundance infections and coverage assessment for detection from reads

Definitive viral diagnosis from metagenomic clinical samples can be extremely challenging due to low sequencing depth, large amounts of host reads and low infectious titres, especially in blood or CSF. One way to distinguish a positive viral diagnosis is to look at alignment coverage against one or multiple reference sequences. When only few reads map to the reference, and genome coverage is low, positive infections often display multiple distinct alignment regions, as opposed to reads mapping to a single or few regions on the reference. [De Vries et al. (2021)](https://www.sciencedirect.com/science/article/pii/S1386653221000792) summarize this concept succinctly in this figure (adapted):

![devries](https://user-images.githubusercontent.com/12873366/158775480-447d847e-5b0d-487c-a39a-81bdf428e09d.png)

Positive calls in these cases can be made from coverage plots showing the distinct alignment regions and a threshold on the number of regions is chosen by the authors (> 3). However, coverage plots require visual assessment and may not be suitable for flagging potential hits in automated pipelines or summary reports. 

### Genomic neighbor subtyping schemes and genotype inference from assemblies

TBD

## Etymology

Not a very creative abbreviation of "virus coverage" but the little spectacles in the logo are a reference to [Rudolf Virchow](https://en.wikipedia.org/wiki/Rudolf_Virchow) who described such trivial concepts as cells, cancer and pathology. His surname is pronounced like `vircov` if you mumble the terminal `v`.

