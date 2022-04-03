use crate::covplot::CovPlot;
use anyhow::Result;
use crossterm::style::Color;
use itertools::Itertools;
use noodles::fasta::{Reader as FastaReader, Record as FastaRecord};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::{bam, bam::record::Cigar, bam::HeaderView, bam::Read};
use rust_lapper::{Interval, Lapper};
use std::collections::btree_map::Entry;
use std::collections::{BTreeMap, HashMap};
use std::convert::TryInto;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::str::from_utf8;
use tabled::{Column, MaxWidth, Modify, Style, Table, Tabled};
use thiserror::Error;

/*
========================
Custom error definitions
========================
*/

#[derive(Error, Debug)]
pub enum ReadAlignmentError {
    /// Indicates failure to read file from command line option
    #[error("failed to read file from option")]
    FileInputError(),
    /// Indicates failure to read file
    #[error("failed to read file")]
    FileIO(#[from] std::io::Error),
    /// Indicates a failure to get sequence length require for coverage plots
    #[error("failed to get sequence length for coverage plot")]
    CovPlotSeqLengthError(),
    /// Indicates failure with the coverage plot module
    #[error("failed to generate data for coverage plot")]
    CovPlot(#[from] crate::covplot::CovPlotError),
    /// Indicates failure to parse a BAM file
    #[error("failed to parse records from BAM")]
    HTSLIBError(#[from] rust_htslib::errors::Error),
    /// Indicates failure to parse a record name from BAM file
    #[error("failed to parse record name from BAM")]
    UTF8Error(#[from] std::str::Utf8Error),
    /// Indicates failure to parse a target name from BAM file
    #[error("failed to parse a valid record target name from BAM")]
    TIDError(#[from] std::num::TryFromIntError),
    /// Indicates failure to parse an u64 from PAF
    #[error("failed to parse a valid integer from PAF")]
    PafRecordIntError(#[from] std::num::ParseIntError),
    /// Indicates failure to conduct grouping because no reference sequences were parsed
    #[error("failed to group outputs due to missing reference sequences")]
    GroupSequenceError(),
    /// Indicates failure to plot coverage when data is grouped
    #[error("coverage plots are not enabled when grouping output")]
    GroupCovPlotError(),
}

/*
=====================
Table display helpers
=====================
*/

fn display_coverage(cov: &f64) -> String {
    format!("{:.4}", cov)
}

/// A struct for computed output fields
#[derive(Debug, Clone, PartialEq)]
pub struct CoverageFields<'a> {
    /// Name of the target sequence
    name: String,
    /// Number of non-overlapping alignment regions
    regions: u64,
    /// Number of unique reads aligned
    reads: u64,
    /// Number of alignments
    alignments: u64,
    /// Number of bases covered by alignments
    bases: u64,
    /// Length of the target sequence
    length: u64,
    /// Fractional coverage of the alignments
    coverage: f64,
    /// Descriptions of the reference sequence headers
    description: String,
    /// Tags for alignment regions in start:stop:aln format
    tags: String,
    /// Unique read identifiers for grouped operations
    unique_reads: Vec<&'a String>,
}

impl CoverageFields<'_> {
    fn as_tag(&self) -> String {
        format!(
            "{:}:{:}:{:}:{:}:{:}:{:}:{:}",
            self.name,
            self.regions,
            self.reads,
            self.alignments,
            self.bases,
            self.length,
            self.coverage
        )
    }
}

/// A struct for computed output fields for display
#[derive(Debug, Clone, PartialEq, Tabled)]
pub struct CoverageTableFields {
    #[header("Sequence")]
    /// Name of the target sequence
    name: String,
    #[header("Regions")]
    /// Number of non-overlapping alignment regions
    regions: u64,
    #[header("Aligned Reads")]
    /// Number of unique reads aligned
    reads: u64,
    #[header("Alignments")]
    /// Number of alignments
    alignments: u64,
    #[header("Covered Bases")]
    /// Number of bases covered by alignments
    bases: u64,
    #[header("Total Bases")]
    /// Length of the target sequence
    length: u64,
    #[header("Coverage")]
    #[field(display_with = "display_coverage")]
    /// Fractional coverage of the alignments
    coverage: f64,
    /// Descriptions of the reference sequence headers
    #[header("Description")]
    description: String,
    /// Tags for alignment regions in start:stop:aln format
    #[header("Tags")]
    tags: String,
}

impl CoverageTableFields {
    /// Takes a coverage field and subsets it to fields
    /// suitable for tabled output
    pub fn from(coverage_fields: &CoverageFields) -> Self {
        Self {
            name: coverage_fields.name.clone(),
            regions: coverage_fields.regions,
            reads: coverage_fields.reads,
            alignments: coverage_fields.alignments,
            bases: coverage_fields.bases,
            length: coverage_fields.length,
            coverage: coverage_fields.coverage,
            description: coverage_fields.description.clone(),
            tags: coverage_fields.tags.clone(),
        }
    }
}

/*
==========================================
Alignment parsing and interval extraction
==========================================
*/

type TargetIntervals = Vec<(String, Lapper<usize, String>)>;

// Parse an optional FASTA file into a optional HashMap
// with sequence name (key) and sequence record (value)
fn parse_fasta(
    fasta: Option<PathBuf>,
) -> Result<Option<HashMap<String, FastaRecord>>, ReadAlignmentError> {
    match fasta {
        Some(fasta_path) => {
            let mut reader = File::open(fasta_path)
                .map(BufReader::new)
                .map(FastaReader::new)?;

            let mut lengths: HashMap<String, FastaRecord> = HashMap::new();
            for result in reader.records() {
                let fasta_record = result?;
                lengths.insert(fasta_record.name().to_string(), fasta_record);
            }
            Ok(Some(lengths))
        }
        None => Ok(None),
    }
}

/// Paf struct
#[derive(Debug, Clone)]
pub struct ReadAlignment {
    /// PafRecords parsed from file (PAF)
    pub target_intervals: TargetIntervals,
    /// Reference sequence names and lengths from file (FASTA)
    pub target_sequences: Option<HashMap<String, FastaRecord>>,
}

impl ReadAlignment {
    // Parse alignments from file
    pub fn from_paf(
        // Path to alignment file [PAF or "-"]
        path: PathBuf,
        // Reference sequence fasta file
        fasta: Option<PathBuf>,
        // Minimum query alignment length
        min_qaln_len: u64,
        // Minimum query coverage
        min_qaln_cov: f64,
        // Minimum mapping quality
        min_mapq: u8,
    ) -> Result<Self, ReadAlignmentError> {
        let reader: Box<dyn BufRead> = match path.file_name() {
            Some(os_str) => match os_str.to_str() {
                Some("-") => Box::new(BufReader::new(std::io::stdin())),
                Some(_) => Box::new(BufReader::new(File::open(path)?)),
                None => return Err(ReadAlignmentError::FileInputError()),
            },
            None => return Err(ReadAlignmentError::FileInputError()),
        };

        let mut target_intervals: BTreeMap<String, Vec<Interval<usize, String>>> = BTreeMap::new();

        for result in reader.lines() {
            let record: PafRecord = PafRecord::from_str(result?)?;
            if record.query_aligned_length() >= min_qaln_len
                && record.query_coverage() >= min_qaln_cov
                && record.mapq >= min_mapq
            {
                match target_intervals.entry(record.tname.clone()) {
                    Entry::Occupied(mut entry) => {
                        entry.get_mut().push(Interval {
                            start: record.tstart,
                            stop: record.tend,
                            val: record.qname, // add the query read name here for reference, see if it affects performance
                        });
                    }
                    Entry::Vacant(entry) => {
                        entry.insert(vec![Interval {
                            start: record.tstart,
                            stop: record.tend,
                            val: record.qname,
                        }]);
                    }
                }
            }
        }

        let target_lappers = target_intervals
            .into_iter()
            .map(|entry| (entry.0, Lapper::new(entry.1)))
            .collect::<TargetIntervals>();

        match parse_fasta(fasta)? {
            None => Ok(Self {
                target_intervals: target_lappers,
                target_sequences: None,
            }),
            Some(fasta_records) => Ok(Self {
                target_intervals: target_lappers,
                target_sequences: Some(fasta_records),
            }),
        }
    }
    // Parse alignments from file
    pub fn from_bam(
        // Path to alignment file [SAM/BAM/CRAM or "-"]
        path: PathBuf,
        // Reference sequence fasta file
        fasta: Option<PathBuf>,
        // Minimum query alignment length
        min_qaln_len: u64,
        // Minimum query coverage
        min_qaln_cov: f64,
        // Minimum mapping quality
        min_mapq: u8,
    ) -> Result<Self, ReadAlignmentError> {
        let mut bam: bam::Reader = match path.file_name() {
            Some(os_str) => match os_str.to_str() {
                Some("-") => bam::Reader::from_stdin()?,
                Some(_) => bam::Reader::from_path(path)?,
                None => return Err(ReadAlignmentError::FileInputError()),
            },
            None => return Err(ReadAlignmentError::FileInputError()),
        };
        let header_view = bam.header().to_owned();

        let mut target_intervals: BTreeMap<String, Vec<Interval<usize, String>>> = BTreeMap::new();

        for result in bam.records() {
            let record = result?;
            if record.is_unmapped() {
                continue;
            }

            let bam_record = BamRecord::from(&record, &header_view)?;

            if bam_record.qalen >= min_qaln_len
                && bam_record.query_coverage() >= min_qaln_cov
                && bam_record.mapq >= min_mapq
            {
                match target_intervals.entry(bam_record.tname.clone()) {
                    Entry::Occupied(mut entry) => {
                        entry.get_mut().push(Interval {
                            start: bam_record.tstart,
                            stop: bam_record.tend,
                            val: bam_record.qname, // add the query read name here for reference, see if it affects performance
                        });
                    }
                    Entry::Vacant(entry) => {
                        entry.insert(vec![Interval {
                            start: bam_record.tstart,
                            stop: bam_record.tend,
                            val: bam_record.qname,
                        }]);
                    }
                }
            }
        }

        let target_lappers = target_intervals
            .into_iter()
            .map(|entry| (entry.0, Lapper::new(entry.1)))
            .collect::<TargetIntervals>();

        match parse_fasta(fasta)? {
            None => Ok(Self {
                target_intervals: target_lappers,
                target_sequences: None,
            }),
            Some(fasta_records) => Ok(Self {
                target_intervals: target_lappers,
                target_sequences: Some(fasta_records),
            }),
        }
    }

    /// Compute coverage distribution by target sequence
    pub fn coverage_statistics(
        &self,
        cov_reg: u64,
        seq_len: u64,
        coverage: f64,
        reads: u64,
        group_by: &Option<String>,
        verbosity: u64,
    ) -> Result<Vec<CoverageFields>, ReadAlignmentError> {
        let mut coverage_fields: Vec<CoverageFields> = Vec::new();
        for (target_name, targets) in &self.target_intervals {
            // Bases of target sequence covered
            let target_cov_bp = targets.cov() as u64;

            let target_record = match &self.target_sequences {
                None => None,
                Some(seqs) => seqs.get(target_name),
            };

            let (target_cov, target_len) = match target_record {
                None => (0_f64, 0_u64),
                Some(fasta_record) => match fasta_record.sequence().len() {
                    0 => (0_f64, 0_u64),
                    _ => (
                        target_cov_bp as f64 / fasta_record.sequence().len() as f64,
                        fasta_record.sequence().len() as u64,
                    ),
                },
            };

            let mut merged_targets = targets.clone();
            merged_targets.merge_overlaps();

            let target_cov_n = merged_targets.intervals.len() as u64;

            let target_tags = match verbosity {
                2 => {
                    // Get the merged target interval data [start:stop:count]
                    let merged_target_intervals = merged_targets
                        .depth()
                        .collect::<Vec<Interval<usize, usize>>>();
                    let merged_target_interval_strings = merged_target_intervals
                        .iter()
                        .map(|i| {
                            let _start = i.start;
                            let _stop = i.stop;
                            let _count = targets.count(i.start, i.stop);
                            format!("{_start}:{_stop}:{_count}") // format tag strings here
                        })
                        .collect::<Vec<String>>();
                    merged_target_interval_strings.join(" ")
                }
                _ => "-".to_string(),
            };

            let target_description = match (verbosity, group_by) {
                (0, Some(_)) | (1, _) | (2, _) | (_, Some(_)) => match target_record {
                    None => "-".to_string(),
                    Some(fasta_record) => match fasta_record.description() {
                        None => "-".to_string(),
                        Some(descr) => descr.to_owned(),
                    },
                },

                (_, None) => "-".to_string(),
            };

            match group_by {
                // If we are not grouping, ignore read identifiers
                None => {
                    // Get the number of unique reads in the alignments
                    let unique_reads = targets
                        .iter()
                        .map(|interval| interval.val.to_owned())
                        .unique()
                        .count() as u64;

                    if target_len >= seq_len
                        && target_cov_n >= cov_reg
                        && target_cov >= coverage
                        && unique_reads >= reads
                    {
                        coverage_fields.push(CoverageFields {
                            name: target_name.to_owned(),
                            regions: target_cov_n,
                            reads: unique_reads,
                            alignments: targets.len() as u64,
                            bases: target_cov_bp,
                            length: target_len,
                            coverage: target_cov,
                            description: target_description,
                            unique_reads: Vec::new(),
                            tags: target_tags,
                        });
                    }
                }
                // If we are grouping, store the read identifiers to summarize for grouped outputs
                Some(_) => {
                    // Get the unique read identifiers in the alignments
                    let unique_read_ids = targets
                        .iter()
                        .map(|interval| &interval.val)
                        .unique()
                        .collect::<Vec<&String>>();
                    let unique_reads = unique_read_ids.len() as u64;

                    if target_len >= seq_len && unique_reads >= reads {
                        coverage_fields.push(CoverageFields {
                            name: target_name.to_owned(),
                            regions: target_cov_n,
                            reads: unique_reads,
                            alignments: targets.len() as u64,
                            bases: target_cov_bp,
                            length: target_len,
                            coverage: target_cov,
                            description: target_description,
                            unique_reads: unique_read_ids,
                            tags: target_tags,
                        });
                    }
                }
            }
        }

        Ok(coverage_fields)
    }

    #[cfg(not(tarpaulin_include))]
    /// Print the coverage statistics to console, alternatively in pretty table format
    pub fn to_console(
        &self,
        coverage_fields: &[CoverageFields],
        table: bool,
    ) -> Result<(), ReadAlignmentError> {
        match table {
            true => {
                let table_fields = coverage_fields
                    .iter()
                    .map(|x| CoverageTableFields::from(x))
                    .collect::<Vec<CoverageTableFields>>();
                let _table = Table::new(table_fields)
                    .with(Modify::new(Column(7..)).with(MaxWidth::wrapping(32)))
                    .with(Style::modern());
                println!("{}", _table)
            }
            false => {
                for cov_fields in coverage_fields {
                    println!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}",
                        cov_fields.name,
                        cov_fields.regions,
                        cov_fields.reads,
                        cov_fields.alignments,
                        cov_fields.bases,
                        cov_fields.length,
                        cov_fields.coverage,
                        cov_fields.description,
                        cov_fields.tags
                    )
                }
            }
        }

        Ok(())
    }

    #[cfg(not(tarpaulin_include))]
    /// Print the coverage distributions to console as a text-based coverage plot
    pub fn coverage_plots(
        &self,
        coverage_fields: &[CoverageFields],
        max_width: u64,
    ) -> Result<(), ReadAlignmentError> {
        println!();
        for (target_name, targets) in &self.target_intervals {
            let target_record = match &self.target_sequences {
                None => None,
                Some(seqs) => seqs.get(target_name),
            };
            let target_len = match target_record {
                None => return Err(ReadAlignmentError::CovPlotSeqLengthError()),
                Some(fasta_record) => fasta_record.sequence().len() as u64,
            };

            if coverage_fields
                .iter()
                .map(|x| &x.name)
                .any(|x| x == target_name)
            {
                let covplot = CovPlot::new(targets, target_len, max_width)?;
                covplot.to_console(target_name, target_len, Color::Red)?;
            }
        }

        Ok(())
    }

    pub fn group_output(
        &self,
        coverage_fields: &[CoverageFields],
        cov_reg: u64,
        coverage: f64,
        group_by: String,
        group_sep: String,
    ) -> Result<Vec<CoverageFields>, ReadAlignmentError> {
        let mut grouped_fields: BTreeMap<String, Vec<&CoverageFields>> = BTreeMap::new();

        for cov_field in coverage_fields {
            let groups = cov_field.description.split(&group_sep);

            for group in groups {
                if group.contains(&group_by) {
                    match grouped_fields.entry(group.to_string()) {
                        Entry::Occupied(mut entry) => {
                            entry.get_mut().push(cov_field);
                        }
                        Entry::Vacant(entry) => {
                            entry.insert(vec![&cov_field]);
                        }
                    }
                }
            }
        }

        let mut grouped_coverage_fields: Vec<CoverageFields> = Vec::new();
        for (group, fields) in grouped_fields {
            let mut grouped_fields = CoverageFields {
                name: format!("{} ({})", group.trim().replace(&group_by, ""), fields.len()),
                regions: 0,
                reads: 0,
                alignments: 0,
                bases: 0,
                length: 0,
                coverage: 0.,
                description: String::new(),
                tags: String::new(),
                unique_reads: Vec::new(),
            };

            let mut ureads = Vec::new();

            for field in &fields {
                grouped_fields.regions += field.regions;
                grouped_fields.alignments += field.alignments;
                grouped_fields.coverage += field.coverage;

                if !grouped_fields.description.contains(&field.description) {
                    grouped_fields.description.push_str(&field.description);
                }

                match field.tags.as_str() {
                    "-" => {
                        if !grouped_fields.tags.contains("-") {
                            grouped_fields.tags.push_str("-");
                        }
                    }
                    _ => {
                        grouped_fields.tags.push_str(&field.as_tag());
                        grouped_fields.tags.push_str(" ");
                    }
                }
                for read in &field.unique_reads {
                    ureads.push(read)
                }
            }

            grouped_fields.coverage = grouped_fields.coverage / fields.len() as f64;

            let reads = ureads.iter().unique().count();

            grouped_fields.reads = reads as u64;

            if grouped_fields.regions >= cov_reg && grouped_fields.coverage >= coverage {
                grouped_coverage_fields.push(grouped_fields);
            }
        }

        Ok(grouped_coverage_fields)
    }
}

/*
=================
Alignment records
=================
*/

/// Return the query alignment length from a CIGAR string
/// as the sum of all matches (M) and insertions (I).
///
/// PAF considers insertaions but not deletions when computing
/// query start and end positions from which the query alignment
/// length is then calculated in PafRecord (qend - qstart).
fn qalen_from_cigar<'a>(cigar: impl Iterator<Item = &'a Cigar>) -> u32 {
    cigar
        .map(|x| match x {
            Cigar::Match(_) => x.len(),
            Cigar::Ins(_) => x.len(),
            _ => 0,
        })
        .sum()
}

#[derive(Debug, Clone)]
pub struct BamRecord {
    /// Query sequence name.
    pub qname: String,
    /// Query sequence length.
    pub qlen: u64,
    /// Query alignment length.
    pub qalen: u64,
    /// Target sequence name.
    pub tname: String,
    /// Target start on original strand (0-based).
    pub tstart: usize,
    /// Target end on original strand (0-based).
    pub tend: usize,
    /// Mapping quality (0-255; 255 for missing).
    pub mapq: u8,
}

impl BamRecord {
    /// Create a new (reduced) BamRecord from a BAM HTS LIB record
    pub fn from(record: &bam::Record, header: &HeaderView) -> Result<Self, ReadAlignmentError> {
        let tid: u32 = record.tid().try_into()?; // from i32 [-1 allowed, but should not occur]
        let tname = from_utf8(header.tid2name(tid))?.to_string();
        let qname = from_utf8(record.qname())?.to_string();

        let tstart: usize = record.reference_start().try_into()?; // from i32 [-1 allowed, but should not occur]
        let tend: usize = record.reference_end().try_into()?; // from i32 [-1 allowed, but should not occur]

        let qlen = record.seq_len() as u64;
        let mapq = record.mapq();

        let qalen: u32 = qalen_from_cigar(record.cigar().iter());

        Ok(Self {
            qname,
            qlen,
            qalen: qalen as u64,
            tname,
            tstart,
            tend,
            mapq,
        })
    }
    /// Coverage of the aligned query sequence.
    pub fn query_coverage(&self) -> f64 {
        match self.qlen == 0 {
            true => 0f64,
            false => self.qalen as f64 / self.qlen as f64,
        }
    }
}

/// PAF record without tags
#[derive(Debug, Clone)]
pub struct PafRecord {
    /// Query sequence name.
    pub qname: String,
    /// Query sequence length.
    pub qlen: u64,
    /// Query start (0-based; BED-like; closed).
    pub qstart: usize,
    /// Query end (0-based; BED-like; open).
    pub qend: usize,
    /// Relative strand: "+" or "-".
    pub strand: String,
    /// Target sequence name.
    pub tname: String,
    /// Target sequence length.
    pub tlen: u64,
    /// Target start on original strand (0-based).
    pub tstart: usize,
    /// Target end on original strand (0-based).
    pub tend: usize,
    /// Number of matching bases in the mapping.
    pub mlen: u64,
    /// Alignment block length. Number of bases, including gaps, in the mapping.
    pub blen: u64,
    /// Mapping quality (0-255; 255 for missing).
    pub mapq: u8,
}

impl PafRecord {
    // Create a record from a parsed line
    pub fn from_str(paf: String) -> Result<Self, ReadAlignmentError> {
        let fields: Vec<&str> = paf.split('\t').collect();

        let record = Self {
            qname: fields[0].to_string(),
            qlen: fields[1].parse::<u64>()?,
            qstart: fields[2].parse::<usize>()?,
            qend: fields[3].parse::<usize>()?,
            strand: fields[4].to_string(),
            tname: fields[5].to_string(),
            tlen: fields[6].parse::<u64>()?,
            tstart: fields[7].parse::<usize>()?,
            tend: fields[8].parse::<usize>()?,
            mlen: fields[9].parse::<u64>()?,
            blen: fields[10].parse::<u64>()?,
            mapq: fields[11].parse::<u8>()?,
        };

        Ok(record)
    }
    /// Length of the aligned query sequence.
    pub fn query_aligned_length(&self) -> u64 {
        (self.qend - self.qstart) as u64
    }
    /// Coverage of the aligned query sequence.
    pub fn query_coverage(&self) -> f64 {
        match self.qlen == 0 {
            true => 0f64,
            false => self.query_aligned_length() as f64 / self.qlen as f64,
        }
    }
}

#[cfg(test)]
#[cfg(not(tarpaulin_include))]
mod tests {

    use super::*;
    use std::path::PathBuf;

    /*
    ===============
      Test cases
    ===============
    */

    struct TestCases {
        // Valid PafRecord
        paf_record_ok: PafRecord,
        // BAM record ok
        bam_record_ok: BamRecord,
        // Valid PAF string, sufficient fields to parse
        paf_test_str_ok: String,
        // Invalid PAF string, too few fields to parse
        paf_test_str_size_fail: String,
        // Valid PAF file
        paf_test_file_ok: PathBuf,
        // Valid FASTA file
        paf_test_fasta_ok: PathBuf,
        // Invalid PAF format, record has too few fields
        paf_test_file_record_size_fail: PathBuf,
        // PAF test alignment intervals for L segment sequence
        paf_test_intervals_l_segment: Vec<Interval<usize, String>>,
        // PAF test alignment intervals for S segment sequence
        paf_test_intervals_s_segment: Vec<Interval<usize, String>>,
        // PAF test alignment extracted coverage statistics, without reference sequences
        paf_test_coverage_statistics_no_ref: Vec<CoverageFields>,
        // PAF test alignment extracted coverage statistics, without reference sequences, without tags (verbosity: 0)
        paf_test_coverage_statistics_no_ref_no_tags: Vec<CoverageFields>,
        // PAF test alignment extracted coverage statistics, with reference sequences
        paf_test_coverage_statistics_ref: Vec<CoverageFields>,

        // Valid SAM file
        sam_test_file_ok: PathBuf,
        // Valid BAM file
        bam_test_file_ok: PathBuf,
        // Valid FASTA file
        bam_test_fasta_ok: PathBuf,
        // Valid FASTA but without sequence edge case
        bam_test_fasta_zero_ok: PathBuf,

        // General failure to parse input file name
        input_file_name_fail: PathBuf,
    }

    impl TestCases {
        fn new() -> Self {
            Self {
                paf_record_ok: PafRecord {
                    qname: "query".to_string(),
                    qlen: 4,
                    qstart: 400,
                    qend: 404,
                    strand: "+".to_string(),
                    tname: "target".to_string(),
                    tlen: 5,
                    tstart: 500,
                    tend: 504,
                    mlen: 4,
                    blen: 4,
                    mapq: 60,
                },
                bam_record_ok: BamRecord {
                    qname: "query".to_string(),
                    qlen: 4,
                    qalen: 4,
                    tname: "target".to_string(),
                    tstart: 500,
                    tend: 504,
                    mapq: 60,
                },
                paf_test_str_ok: String::from(
                    "query\t4\t400\t404\t+\ttarget\t5\t500\t504\t4\t4\t60",
                ),
                paf_test_str_size_fail: String::from(
                    "query\t4\t400\t404\t+\ttarget\t5\t500\t504\t4\t4",
                ),
                paf_test_file_ok: PathBuf::from("tests/cases/test_ok.paf"),
                paf_test_fasta_ok: PathBuf::from("tests/cases/test_paf_ok.fasta"),
                paf_test_file_record_size_fail: PathBuf::from(
                    "tests/cases/test_record_size_fail.paf",
                ),
                paf_test_intervals_l_segment: vec![
                    Interval {
                        start: 1786,
                        stop: 1834,
                        val: "FS10001392:17:BPL20314-1135:1:1113:8980:1660".to_string(),
                    },
                    Interval {
                        start: 4538,
                        stop: 4665,
                        val: "FS10001392:17:BPL20314-1135:1:1101:5600:2170".to_string(),
                    },
                    Interval {
                        start: 4758,
                        stop: 4904,
                        val: "FS10001392:17:BPL20314-1135:1:1101:5600:2170".to_string(),
                    },
                ],
                paf_test_intervals_s_segment: vec![
                    Interval {
                        start: 1574,
                        stop: 1671,
                        val: "FS10001392:17:BPL20314-1135:1:1108:6180:1130".to_string(),
                    },
                    Interval {
                        start: 2188,
                        stop: 2228,
                        val: "FS10001392:17:BPL20314-1135:1:1116:1700:4010".to_string(),
                    },
                ],
                paf_test_coverage_statistics_no_ref: vec![
                    CoverageFields {
                        name: "21172389_LCMV_L-segment_final".to_string(),
                        regions: 3,
                        reads: 2,
                        alignments: 3,
                        bases: 321,
                        length: 0,
                        coverage: 0.,
                        description: "-".to_string(),
                        tags: "1786:1834:1 4538:4665:1 4758:4904:1".to_string(),
                    },
                    CoverageFields {
                        name: "21172389_LCMV_S-segment_final".to_string(),
                        regions: 2,
                        reads: 2,
                        alignments: 2,
                        bases: 137,
                        length: 0,
                        coverage: 0.,
                        description: "-".to_string(),
                        tags: "1574:1671:1 2188:2228:1".to_string(),
                    },
                ],
                paf_test_coverage_statistics_no_ref_no_tags: vec![
                    CoverageFields {
                        name: "21172389_LCMV_L-segment_final".to_string(),
                        regions: 3,
                        reads: 2,
                        alignments: 3,
                        bases: 321,
                        length: 0,
                        coverage: 0.,
                        description: "-".to_string(),
                        tags: "-".to_string(),
                    },
                    CoverageFields {
                        name: "21172389_LCMV_S-segment_final".to_string(),
                        regions: 2,
                        reads: 2,
                        alignments: 2,
                        bases: 137,
                        length: 0,
                        coverage: 0.,
                        description: "-".to_string(),
                        tags: "-".to_string(),
                    },
                ],
                paf_test_coverage_statistics_ref: vec![
                    CoverageFields {
                        name: "21172389_LCMV_L-segment_final".to_string(),
                        regions: 3,
                        reads: 2,
                        alignments: 3,
                        bases: 321,
                        length: 7194,
                        coverage: 0.044620517097581316,
                        description: "-".to_string(),
                        tags: "1786:1834:1 4538:4665:1 4758:4904:1".to_string(),
                    },
                    CoverageFields {
                        name: "21172389_LCMV_S-segment_final".to_string(),
                        regions: 2,
                        reads: 2,
                        alignments: 2,
                        bases: 137,
                        length: 3407,
                        coverage: 0.040211329615497504,
                        description: "-".to_string(),
                        tags: "1574:1671:1 2188:2228:1".to_string(),
                    },
                ],

                sam_test_file_ok: PathBuf::from("tests/cases/test_ok.sam"),
                bam_test_file_ok: PathBuf::from("tests/cases/test_ok.bam"),
                bam_test_fasta_ok: PathBuf::from("tests/cases/test_bam_ok.fasta"),
                bam_test_fasta_zero_ok: PathBuf::from("tests/cases/test_bam_zero_ok.fasta"),

                input_file_name_fail: PathBuf::from("tests/cases/.."),
            }
        }
    }

    #[test]
    fn paf_record_from_str_ok() {
        let test_cases = TestCases::new();
        let record = PafRecord::from_str(test_cases.paf_test_str_ok).unwrap();

        assert_eq!(record.qname, test_cases.paf_record_ok.qname);
        assert_eq!(record.qlen, test_cases.paf_record_ok.qlen);
        assert_eq!(record.qstart, test_cases.paf_record_ok.qstart);
        assert_eq!(record.qend, test_cases.paf_record_ok.qend);
        assert_eq!(record.strand, test_cases.paf_record_ok.strand);
        assert_eq!(record.tname, test_cases.paf_record_ok.tname);
        assert_eq!(record.tlen, test_cases.paf_record_ok.tlen);
        assert_eq!(record.tstart, test_cases.paf_record_ok.tstart);
        assert_eq!(record.tend, test_cases.paf_record_ok.tend);
        assert_eq!(record.mlen, test_cases.paf_record_ok.blen);
        assert_eq!(record.mapq, test_cases.paf_record_ok.mapq);
    }
    #[test]
    #[should_panic]
    fn paf_record_from_str_size_fail() {
        let test_cases = TestCases::new();
        PafRecord::from_str(test_cases.paf_test_str_size_fail).unwrap();
    }
    #[test]
    fn paf_record_query_aligned_length_ok() {
        let test_cases = TestCases::new();
        assert_eq!(test_cases.paf_record_ok.query_aligned_length(), 4)
    }
    #[test]
    fn paf_record_query_coverage_ok() {
        let test_cases = TestCases::new();
        assert_eq!(test_cases.paf_record_ok.query_coverage(), 1.0_f64)
    }
    #[test]
    fn paf_record_query_coverage_zero_len_ok() {
        let test_cases = TestCases::new();
        let mut record = test_cases.paf_record_ok;
        record.qlen = 0;
        assert_eq!(record.query_coverage(), 0_f64)
    }
    #[test]
    #[should_panic]
    fn paf_parser_input_file_name_fail() {
        let test_cases = TestCases::new();
        ReadAlignment::from_paf(
            test_cases.input_file_name_fail,
            Some(test_cases.paf_test_fasta_ok),
            0_u64,
            0_f64,
            0_u8,
        )
        .unwrap();
    }
    #[test]
    #[should_panic]
    fn paf_parser_create_new_record_size_fail() {
        let test_cases = TestCases::new();
        // PafAlignmentError does not implement PartialEq to assure standard Error type for parsing, let test fail instead
        let _ = ReadAlignment::from_paf(
            test_cases.paf_test_file_record_size_fail,
            None,
            0_u64,
            0_f64,
            0_u8,
        )
        .unwrap();
    }
    #[test]
    fn paf_parser_create_new_filter_mapq_ok() {
        let test_cases = TestCases::new();
        let paf_aln =
            ReadAlignment::from_paf(test_cases.paf_test_file_ok, None, 0_u64, 0_f64, 30_u8)
                .unwrap();
        assert_eq!(paf_aln.target_intervals[0].1.len(), 2);
    }
    #[test]
    fn paf_parser_create_new_filter_min_len_ok() {
        let test_cases = TestCases::new();
        let paf_aln =
            ReadAlignment::from_paf(test_cases.paf_test_file_ok, None, 50_u64, 0_f64, 0_u8)
                .unwrap();
        assert_eq!(paf_aln.target_intervals[0].1.len(), 2);
        assert_eq!(paf_aln.target_intervals[1].1.len(), 1);
    }
    #[test]
    fn paf_parser_create_new_filter_min_cov_ok() {
        let test_cases = TestCases::new();
        let paf_aln =
            ReadAlignment::from_paf(test_cases.paf_test_file_ok, None, 0_u64, 0.5, 0_u8).unwrap();
        assert_eq!(paf_aln.target_intervals[0].1.len(), 2);
        assert_eq!(paf_aln.target_intervals[1].1.len(), 1);
    }

    #[test]
    fn paf_parser_create_new_fasta_input_not_provided_ok() {
        let test_cases = TestCases::new();
        ReadAlignment::from_paf(test_cases.paf_test_file_ok, None, 0_u64, 0_f64, 0_u8).unwrap();
    }

    #[test]
    #[should_panic]
    fn paf_parser_create_new_fasta_input_not_exists_fail() {
        let test_cases = TestCases::new();
        let _ = ReadAlignment::from_paf(
            test_cases.paf_test_file_ok,
            Some(PathBuf::from("tests/cases/not_exist.fasta")),
            0_u64,
            0_f64,
            0_u8,
        )
        .unwrap();
    }

    #[test]
    fn paf_parser_compute_target_intervals_ok() {
        let test_cases = TestCases::new();
        let paf_aln =
            ReadAlignment::from_paf(test_cases.paf_test_file_ok, None, 0_u64, 0_f64, 0_u8).unwrap();
        let actual_intervals_l = paf_aln.target_intervals[0].1.intervals.clone();
        let actual_intervals_s = paf_aln.target_intervals[1].1.intervals.clone();
        assert_eq!(actual_intervals_l, test_cases.paf_test_intervals_l_segment);
        assert_eq!(actual_intervals_s, test_cases.paf_test_intervals_s_segment);
    }

    #[test]
    fn paf_parser_compute_statistics_no_ref_ok() {
        let test_cases = TestCases::new();
        let paf_aln =
            ReadAlignment::from_paf(test_cases.paf_test_file_ok, None, 0_u64, 0_f64, 0_u8).unwrap();
        let actual_statistics = paf_aln.coverage_statistics(0, 0, 2).unwrap();
        assert_eq!(
            actual_statistics,
            test_cases.paf_test_coverage_statistics_no_ref
        );
    }

    #[test]
    fn paf_parser_compute_statistics_no_ref_no_tags_ok() {
        let test_cases = TestCases::new();
        let paf_aln =
            ReadAlignment::from_paf(test_cases.paf_test_file_ok, None, 0_u64, 0_f64, 0_u8).unwrap();
        let actual_statistics = paf_aln.coverage_statistics(0, 0, 0).unwrap();
        assert_eq!(
            actual_statistics,
            test_cases.paf_test_coverage_statistics_no_ref_no_tags
        );
    }

    #[test]
    fn paf_parser_compute_statistics_ref_ok() {
        let test_cases = TestCases::new();
        let paf_aln = ReadAlignment::from_paf(
            test_cases.paf_test_file_ok,
            Some(test_cases.paf_test_fasta_ok),
            0_u64,
            0_f64,
            0_u8,
        )
        .unwrap();
        let actual_statistics = paf_aln.coverage_statistics(0, 0, 2).unwrap();
        assert_eq!(
            actual_statistics,
            test_cases.paf_test_coverage_statistics_ref
        );
    }

    #[test]
    fn tabled_helper_format_coverage_ok() {
        let expected = display_coverage(&0.333333333);
        assert_eq!("0.3333", expected);
    }

    #[test]
    fn query_alignment_length_from_cigar_match() {
        // Derived from a real example comparison between PAF and BAM of the same alignment
        let cigars = vec![Cigar::Match(149)];
        let qalen: u32 = qalen_from_cigar(cigars.iter());
        assert_eq!(qalen, 149)
    }

    #[test]
    fn query_alignment_length_from_cigar_with_del() {
        // Derived from a real example comparison between PAF and BAM of the same alignment
        // Deletion is not considered
        let cigars = vec![Cigar::Match(8), Cigar::Del(1), Cigar::Match(141)];
        let qalen: u32 = qalen_from_cigar(cigars.iter());
        assert_eq!(qalen, 149)
    }

    #[test]
    fn query_alignment_length_from_cigar_with_ins() {
        // Insertion is considered
        let cigars = vec![Cigar::Match(7), Cigar::Ins(1), Cigar::Match(141)];
        let qalen: u32 = qalen_from_cigar(cigars.iter());
        assert_eq!(qalen, 149)
    }

    #[test]
    fn query_alignment_length_from_cigar_with_indel() {
        // Derived from a real example comparison between PAF and BAM of the same alignment
        // Insertion is considered, Deletion is not considered
        let cigars = vec![
            Cigar::Match(35),
            Cigar::Del(3),
            Cigar::Match(15),
            Cigar::Ins(1),
            Cigar::Match(64),
            Cigar::SoftClip(34),
        ];
        let qalen: u32 = qalen_from_cigar(cigars.iter());
        assert_eq!(qalen, 115)
    }
    #[test]
    fn bam_parser_sam_file_ok() {
        let test_cases = TestCases::new();
        ReadAlignment::from_bam(
            test_cases.sam_test_file_ok,
            Some(test_cases.bam_test_fasta_ok),
            0_u64,
            0_f64,
            0_u8,
        )
        .unwrap();
    }
    #[test]
    fn bam_parser_bam_file_ok() {
        let test_cases = TestCases::new();
        ReadAlignment::from_bam(
            test_cases.bam_test_file_ok,
            Some(test_cases.bam_test_fasta_ok),
            0_u64,
            0_f64,
            0_u8,
        )
        .unwrap();
    }
    #[test]
    fn bam_parser_no_fasta_file_ok() {
        let test_cases = TestCases::new();
        ReadAlignment::from_bam(test_cases.bam_test_file_ok, None, 0_u64, 0_f64, 0_u8).unwrap();
    }
    #[test]
    #[should_panic]
    fn bam_parser_input_file_name_fail() {
        let test_cases = TestCases::new();
        ReadAlignment::from_bam(
            test_cases.input_file_name_fail,
            Some(test_cases.bam_test_fasta_ok),
            0_u64,
            0_f64,
            0_u8,
        )
        .unwrap();
    }
    #[test]
    fn bam_record_query_coverage_zero_len_ok() {
        let test_cases = TestCases::new();
        let mut record = test_cases.bam_record_ok;
        record.qlen = 0;
        assert_eq!(record.query_coverage(), 0_f64)
    }
    #[test]
    fn coverage_statistics_zero_len_ok() {
        let test_cases = TestCases::new();
        let bam_aln = ReadAlignment::from_bam(
            test_cases.bam_test_file_ok,
            Some(test_cases.bam_test_fasta_zero_ok),
            0_u64,
            0_f64,
            0_u8,
        )
        .unwrap();
        bam_aln.coverage_statistics(0, 0, 0).unwrap();
    }
}
