# vircov <a href='https://github.com/esteinig'><img src='docs/vircov.png' align="right" height="270"/></a>

[![build](https://github.com/esteinig/nanoq/actions/workflows/rust-ci.yaml/badge.svg?branch=master)](https://github.com/esteinig/nanoq/actions/workflows/rust-ci.yaml)
[![codecov](https://codecov.io/gh/esteinig/vircov/branch/main/graph/badge.svg?token=RG95F4C6FE)](https://codecov.io/gh/esteinig/vircov)
![](https://img.shields.io/badge/version-0.4.0-black.svg)

Minimal virus genome coverage assessment for metagenomic diagnostics

## Overview


**`v0.4.0`**

- [Purpose](#purpose)
- [Implementation](#implementation)
- [Install](#install)
- [Usage](#usage)
- [Concept](#concept)

## Purpose

Viral metagenomic diagnostics from low-abundance clinical samples can be challenging in the absence of sufficient genome coverage. `Vircov` extracts distinct non-overlapping regions from a reference alignment and generates some helpful statistics. It can be used to flag potential hits without inspection of coverage plots in automated pipelines and reports.

## Implementation

`Vircov` is written in Rust and works with alignments in the standard `PAF` or `SAM/BAM/CRAM` (next release) formats. It is extremely fast and can process alignments against thousands of viral reference genomes in seconds. Basic input filters can be selected to remove spurious alignments, and text-style coverage plots can be printed to the terminal for visual inspection.

`Vircov` is written for implementation in (accredited) metagenomics pipelines for human patients enroled in the `META-GP` network (Australia). As such, it attempts to be production-grade code with high test coverage, continuous integration, and versioned releases with precompiled binaries for Linux and MacOS.


## Install

```bash
git clone https://github.com/esteinig/vircov 
cd vircov && cargo build --release
```

## Usage

```
vircov tests/cases/test_full_ok.paf --fasta tests/cases/test_ok.fasta --table --cov-reg
```

## Tests

```bash
git clone https://github.com/esteinig/vircov 
cd vircov && cargo test && cargo tarpaulin 
```

## Concept

Definitive viral diagnosis from metagenomic clinical samples can be extremely challenging due to low sequencing depth, large amounts of host reads and low infectious titres, especially in blood or CSF. One way to distinguish a positive viral diagnosis is to look at alignment coverage against one or multiple reference sequences. When only few reads map to the reference and when genome coverage is therefore low, positive infections often display multiple distinct alignment regions, as opposed to reads mapping to a single or few regions on the reference.

[De Vries et al. (2021)](https://www.sciencedirect.com/science/article/pii/S1386653221000792) summarize this concept succinctly in this figure (adapted):

![devries](https://user-images.githubusercontent.com/12873366/158775480-447d847e-5b0d-487c-a39a-81bdf428e09d.png)

Positive calls in these cases can be made from coverage plots showing the distinct alignment regions and a threshold on the number of regions is chosen by the authors (> 3). However, coverage plots require visual assessment and may not be suitable for flagging potential hits in automated pipelines or summary reports. 

`Vircov` attempts to make visual inspection and automated flagging easier by counting the distinct (non-overlapping) coverage regions in an alignment and reports some helpful statistics to make an educated call without having to generate coverage plots. 


## Examples

```
```

## Etymology

Not a very creative abbreviation of virus coverage but the little spectacles in the logo are a reference to [Rudolf Virchow](https://en.wikipedia.org/wiki/Rudolf_Virchow) who described such trivial concepts as cells, cancer and pathology. `Vircov` is pronounced somewhat similarly to his surname if you mumble the sharp "v" right at the end.
