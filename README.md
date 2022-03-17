# vircov <a href='https://github.com/esteinig'><img src='docs/vircov.png' align="right" height="270"/></a>

[![build](https://github.com/esteinig/nanoq/actions/workflows/rust-ci.yaml/badge.svg?branch=master)](https://github.com/esteinig/nanoq/actions/workflows/rust-ci.yaml)
[![codecov](https://codecov.io/gh/esteinig/vircov/branch/main/graph/badge.svg?token=RG95F4C6FE)](https://codecov.io/gh/esteinig/vircov)
![](https://img.shields.io/badge/version-0.3.0-black.svg)

Minimal virus genome coverage assessment for metagenomic diagnostics

## Overview

Viral metagenomic diagnostics from low-abundance clinical samples can be challenging in the absence of sufficient genome coverage. `Vircov` extracts distinct non-overlapping regions from a reference alignment and generates some helpful statistics. It can be used to flag potential hits without inspection of coverage plots in automated pipelines and reports.

## Concept

Definitive viral diagnosis from metagenomic clinical samples can be extremely challenging due to low sequence depth, large amounts of host reads and low infectious titres, especially in samples like blood or CSF. One way to distinguish a positive viral diagnosis is to look at alignment coverage against one or multiple reference sequences. Even where genome coverage is low (e.g. < 20%) positive calls confirmed by orthogonal testing methods (e.g. PCR) often display multiple distinct alignment regions, as opposed to reads mapping to a single or few regions on the reference. 

De Vries et al. (2021) summarize this concept succinctly in this figure:




Positive calls in these cases can be made from coverage plots showing the distinct alignment regions, but these require visual assessment and may not be suitable for flagging potential hits in automated pipelines and reports. `Vircov` attempts to make inspection and automated flagging easier, by counting the distinct (non-overlapping) coverage regions in an alignment and reporting their number and coverage of the genome to make an educated call without having to generate coverage plots.




