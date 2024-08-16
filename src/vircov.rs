use crate::alignment::{Aligner, Preset, ReadAlignment};
use crate::error::VircovError;
use crate::subtype::{SubtypeDatabase, SubtypeSummary};
use crate::utils::init_logger;

use rayon::prelude::*;
use anyhow::Result;
use indexmap::IndexMap;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::PathBuf;
use std::io::{BufRead, BufReader};
use std::fmt;


/// The main Vircov application structure.
pub struct Vircov {

}

impl Vircov {

}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VircovConfig {
    pub aligner: AlignerConfig,
    pub reference: ReferenceConfig,
    pub filter: FilterConfig,
    pub coverage: CoverageConfig,
    pub consensus: ConsensusConfig,
    pub subtype: SubtypeConfig
} 


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignerConfig {
    pub aligner: Aligner,
    pub preset: Option<Preset>,
    pub args: Option<String>,
    pub input: Vec<PathBuf>,
    pub paired_end: bool,
    pub index: PathBuf,
    pub output: Option<PathBuf>,
    pub threads: usize
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReferenceConfig {
    pub fasta: Option<PathBuf>,
    pub exclude: Option<PathBuf>,
}
impl Default for ReferenceConfig {
    fn default() -> Self {
        Self {
            fasta: None,
            exclude: None
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FilterConfig {
    pub min_query_length: u64,
    pub min_query_coverage: f64,
    pub min_mapq: u8,
    pub min_alignments: u64,
    pub min_regions: u64,
    pub min_coverage: f64,
    pub min_unique_reads: u64,
    pub min_reference_length: u64,
    pub min_regions_coverage: Option<f64>,

    pub min_grouped_regions: u64,
    pub min_grouped_mean_coverage: f64,
    pub min_grouped_alignments: u64,
    pub min_grouped_reads: u64,
}
impl Default for FilterConfig {
    fn default() -> Self {
        Self {
            min_query_length: 0,
            min_query_coverage: 0.0,
            min_mapq: 0,
            min_alignments: 0,
            min_regions: 0,
            min_coverage: 0.0,
            min_unique_reads: 0,
            min_reference_length: 0,
            min_regions_coverage: None,
            min_grouped_regions: 0,
            min_grouped_mean_coverage: 0.0,
            min_grouped_alignments: 0,
            min_grouped_reads: 0,
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoverageConfig {
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConsensusConfig {
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SubtypeConfig {
}

/// Helper function to read lines from a file into a vector of strings.
fn read_lines_to_vec(filename: &PathBuf) -> Result<Vec<String>> {
    let file = std::fs::File::open(filename)?;
    let reader = BufReader::new(file);
    let mut lines = Vec::new();

    for line in reader.lines() {
        let line = line?;
        lines.push(line);
    }

    Ok(lines)
}

/// Helper function to get supported subtypes. This could be implemented based on actual data.
fn get_supported_subtypes() -> HashMap<&'static str, IndexMap<&'static str, Vec<&'static str>>> {
    HashMap::from([
        ("rsv", IndexMap::from([
            ("RSV-A", Vec::from(["Subgroup A", "virus A isolate", "RSV-A", "RSVA", "/A/", "A-TX", "syncytial virus A"])),
            ("RSV-B", Vec::from(["Subgroup B", "virus B isolate", "RSV-B", "RSVB", "/B/", "B-TX", "B-WaDC", "syncytial virus B"]))
        ])),
        ("rva", IndexMap::from([
            ("regex", Vec::from([r"[Rr]hinovirus A([A-Za-z0-9]+)\b"])),
        ])),
        ("hpiv", IndexMap::from([
            ("HPIV-1", Vec::from(["parainfluenza virus 1", "respirovirus 1", "HPIV1", "hPIV1", "PIV1"])),
            ("HPIV-2", Vec::from(["parainfluenza virus 2", "respirovirus 2", "HPIV2", "orthorubulavirus 2", "hPIV2", "PIV2"])),
            ("HPIV-3", Vec::from(["parainfluenza virus 3", "respirovirus 3", "HPIV3", "hPIV3", "PIV3"])),
            ("HPIV-4", Vec::from(["parainfluenza virus 4", "respirovirus 4", "HPIV4", "orthorubulavirus 4", "hPIV4", "PIV4"]))
        ])),
        ("hmpv", IndexMap::from([
            ("HMPV-A1", Vec::from(["/A1", "type A1", "A1/"])),
            ("HMPV-A2", Vec::from(["/A2", "type A2", "A2/"])),
            ("HMPV-B1", Vec::from(["/B1", "type B1", "B1/"])),
            ("HMPV-B2", Vec::from(["/B2", "type B2", "B2/"])),
            ("HMPV-A", Vec::from(["A/HMPV/Beijing", "type A", "/A,"])),
            ("HMPV-B", Vec::from(["B/HMPV/Beijing", "type B", "/B,"])),
        ])),
        ("hcov", IndexMap::from([
            ("229E", Vec::from(["229E", "Camel alphacoronavirus"])),
            ("HKU1", Vec::from(["HKU1"])),
            ("NL63", Vec::from(["NL63"])),
            ("OC43", Vec::from(["OC43"])),
        ])),
    ])
}
