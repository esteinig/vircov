use crate::covplot::CovPlot;
use anyhow::Result;
use crossterm::style::Color;
use itertools::Itertools;
use noodles::fasta;
use rust_lapper::{Interval, Lapper};
use std::collections::btree_map::Entry;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use tabled::{Column, MaxWidth, Modify, Style, Table, Tabled};
use thiserror::Error;

type TargetIntervals = Vec<(String, Lapper<usize, String>)>;

/*
========================
Custom error definitions
========================
*/

#[derive(Error, Debug)]
pub enum PafAlignmentError {
    /// Indicates failure to read file
    #[error("failed to read file")]
    FileIO(#[from] std::io::Error),
    /// Indicates failure to read a line into record
    #[error("failed to parse record")]
    Record(#[from] PafRecordError),
    /// Indicates a failure to get sequence length require for coverage plots
    #[error("failed to get sequence length for coverage plot")]
    CovPlotSeqLength(),
    /// Indicates failure with the coverage plot module
    #[error("failed to generate data for coverage plot")]
    CovPlot(#[from] crate::covplot::CovPlotError),
}

#[derive(Error, Debug, PartialEq)]
pub enum PafRecordError {
    /// Indicates failure to parse a record due to insufficient fields
    #[error("record has too few fields")]
    RecordSize(),
    /// Indicates failure to parse an integer from a record string
    #[error("failed to parse an integer field")]
    RecordInteger(#[from] std::num::ParseIntError),
}

/*
=====================
Table display helpers
=====================
*/

fn display_coverage(cov: &f64) -> String {
    format!("{:.4}", cov)
}

fn display_tags(tags: &str) -> String {
    match tags {
        "" => "-".to_string(),
        _ => tags.to_owned(),
    }
}

/// A struct for computed output fields
#[derive(Debug, Clone, PartialEq, Tabled)]
pub struct CoverageFields {
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
    /// Tags for alignment regions
    #[header("Tags")]
    #[field(display_with = "display_tags")]
    tags: String,
}

/*
=====================================
PAF alignment parssing and extraction
=====================================
*/

/// Paf struct
#[derive(Debug, Clone)]
pub struct PafAlignment {
    /// PafRecords parsed from file (PAF)
    pub records: Vec<PafRecord>,
    /// Reference sequence names and lengths from file (FASTA)
    pub seq_lengths: Option<HashMap<String, u64>>,
}

impl PafAlignment {
    // Parse alignments from file without filtering
    pub fn from(
        path: PathBuf,
        // Reference sequence fasta
        fasta: Option<PathBuf>,
        // Minimum query alignment length
        min_qaln_len: u64,
        // Minimum query coverage
        min_qaln_cov: f64,
        // Minimum mapping quality
        min_mapq: u8,
    ) -> Result<Self, PafAlignmentError> {
        let paf_file = File::open(path).map(BufReader::new)?;
        let mut records = Vec::new();
        for line in paf_file.lines() {
            let record = PafRecord::from_str(line?)?;
            if record.query_aligned_length()? >= min_qaln_len
                && record.query_coverage()? >= min_qaln_cov
                && record.mapq >= min_mapq
            {
                records.push(record);
            }
        }

        match fasta {
            Some(fasta_path) => {
                let mut reader = File::open(fasta_path)
                    .map(BufReader::new)
                    .map(fasta::Reader::new)?;

                let mut lengths: HashMap<String, u64> = HashMap::new();
                for result in reader.records() {
                    let record = result?;
                    lengths.insert(record.name().to_string(), record.sequence().len() as u64);
                }

                Ok(Self {
                    records,
                    seq_lengths: Some(lengths),
                })
            }
            None => Ok(Self {
                records,
                seq_lengths: None,
            }),
        }
    }
    /// Compute coverage distribution by target sequence
    pub fn coverage_statistics(
        &self,
        verbosity: u64,
    ) -> Result<Vec<CoverageFields>, PafAlignmentError> {
        let mut coverage_fields: Vec<CoverageFields> = Vec::new();
        for (target_name, targets) in self.target_intervals()? {
            // Bases of target sequence covered
            let target_cov_bp = targets.cov();

            let target_seq_len = match &self.seq_lengths {
                None => None,
                Some(seqs) => seqs.get(&target_name),
            };

            let target_seq_cov = match target_seq_len {
                None => 0_f64,
                Some(0) => 0_f64,
                Some(seq_len) => target_cov_bp as f64 / *seq_len as f64,
            };

            let target_seq_len_display = match target_seq_len {
                None => 0,
                Some(seq_len) => *seq_len,
            };

            let mut merged_targets = targets.clone();
            merged_targets.merge_overlaps();

            let target_cov_n = merged_targets.intervals.len();

            let tags = match verbosity {
                1 => {
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
                            format!("{_start}:{_stop}:{_count}")
                        })
                        .collect::<Vec<String>>();
                    merged_target_interval_strings.join(" ")
                }
                _ => "".to_string(),
            };

            // Get the number of unique reads in the alignments
            let unique_reads: u64 = targets
                .iter()
                .map(|interval| interval.val.to_owned())
                .unique()
                .count() as u64;

            coverage_fields.push(CoverageFields {
                name: target_name,
                regions: target_cov_n as u64,
                reads: unique_reads,
                alignments: targets.len() as u64,
                bases: target_cov_bp as u64,
                length: target_seq_len_display,
                coverage: target_seq_cov,
                tags,
            });
        }

        Ok(coverage_fields)
    }

    #[cfg(not(tarpaulin_include))]
    /// Print the coverage statistics to console, alternatively in pretty table format
    pub fn to_console(
        &self,
        coverage_fields: Vec<CoverageFields>,
        pretty: bool,
    ) -> Result<(), PafAlignmentError> {
        match pretty {
            true => {
                let table = Table::new(coverage_fields)
                    .with(Modify::new(Column(7..)).with(MaxWidth::wrapping(32)))
                    .with(Style::modern());
                println!("{}", table)
            }
            false => {
                for cov_fields in coverage_fields {
                    println!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}",
                        cov_fields.name,
                        cov_fields.regions,
                        cov_fields.reads,
                        cov_fields.alignments,
                        cov_fields.bases,
                        cov_fields.length,
                        cov_fields.coverage,
                        cov_fields.tags
                    )
                }
            }
        }

        Ok(())
    }

    #[cfg(not(tarpaulin_include))]
    /// Print the coverage distributions to console as a text-based coverage plot
    pub fn coverage_plots(&self, max_width: u64) -> Result<(), PafAlignmentError> {
        println!();
        for (target_name, targets) in self.target_intervals()? {
            let target_seq_len = match &self.seq_lengths {
                None => None,
                Some(seqs) => seqs.get(&target_name),
            };
            let seq_length = match target_seq_len {
                None => return Err(PafAlignmentError::CovPlotSeqLength()),
                Some(value) => *value,
            };

            let covplot = CovPlot::new(targets, seq_length, max_width)?;

            covplot.to_console(target_name, seq_length, Color::Red)?;
        }

        Ok(())
    }

    /// Get target alignment interval lappers by target sequence
    fn target_intervals(&self) -> Result<TargetIntervals, PafAlignmentError> {
        let mut target_intervals: BTreeMap<String, Vec<Interval<usize, String>>> = BTreeMap::new();
        for record in &self.records {
            match target_intervals.entry(record.tname.clone()) {
                Entry::Occupied(mut entry) => {
                    entry.get_mut().push(Interval {
                        start: record.tstart,
                        stop: record.tend,
                        val: record.qname.clone(), // add the query read name here for reference, see if it affects performance
                    });
                }
                Entry::Vacant(entry) => {
                    entry.insert(vec![Interval {
                        start: record.tstart,
                        stop: record.tend,
                        val: record.qname.clone(),
                    }]);
                }
            }
        }

        let target_lappers = target_intervals
            .into_iter()
            .map(|entry| (entry.0, Lapper::new(entry.1)))
            .collect::<Vec<(String, Lapper<usize, String>)>>();

        Ok(target_lappers)
    }
}

/*
==========
PAF record
==========
*/

/// PAF record without tags
#[derive(Debug, Clone)]
#[allow(dead_code)]
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
    /// Populate a record from a record string, includes some error checks
    pub fn from_str(paf_str: String) -> Result<Self, PafRecordError> {
        let fields: Vec<&str> = paf_str.split('\t').collect();

        if fields.len() < 12 {
            return Err(PafRecordError::RecordSize());
        }

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
    pub fn query_aligned_length(&self) -> Result<u64, PafRecordError> {
        let aln_length = (self.qend - self.qstart) as u64;
        Ok(aln_length)
    }
    /// Coverage of the aligned query sequence.
    pub fn query_coverage(&self) -> Result<f64, PafRecordError> {
        match self.qlen == 0 {
            true => Ok(0f64),
            false => Ok(self.query_aligned_length()? as f64 / self.qlen as f64),
        }
    }
}

#[cfg(test)]
#[cfg(not(tarpaulin_include))]
mod tests {

    use float_eq::float_eq;
    use std::path::PathBuf;

    use super::*;

    /*
    ===============
      Test cases
    ===============
    */

    struct TestCases {
        // Valid PAF record struct instance
        paf_test_record_ok: PafRecord,
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
    }

    impl TestCases {
        fn new() -> Self {
            Self {
                paf_test_record_ok: PafRecord {
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
                paf_test_str_ok: String::from(
                    "query\t4\t400\t404\t+\ttarget\t5\t500\t504\t4\t4\t60",
                ),
                paf_test_str_size_fail: String::from(
                    "query\t4\t400\t404\t+\ttarget\t5\t500\t504\t4\t4",
                ),
                paf_test_file_ok: PathBuf::from("tests/cases/test_ok.paf"),
                paf_test_fasta_ok: PathBuf::from("tests/cases/test_ok.fasta"),
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
                        tags: "".to_string(),
                    },
                    CoverageFields {
                        name: "21172389_LCMV_S-segment_final".to_string(),
                        regions: 2,
                        reads: 2,
                        alignments: 2,
                        bases: 137,
                        length: 0,
                        coverage: 0.,
                        tags: "".to_string(),
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
                        tags: "1574:1671:1 2188:2228:1".to_string(),
                    },
                ],
            }
        }
    }

    /*
    ===============
       PafRecord
    ===============
    */

    #[test]
    fn paf_record_from_str_ok() {
        let test_cases = TestCases::new();
        let record = PafRecord::from_str(test_cases.paf_test_str_ok).unwrap();

        assert_eq!(record.qname, test_cases.paf_test_record_ok.qname);
        assert_eq!(record.qlen, test_cases.paf_test_record_ok.qlen);
        assert_eq!(record.qstart, test_cases.paf_test_record_ok.qstart);
        assert_eq!(record.qend, test_cases.paf_test_record_ok.qend);
        assert_eq!(record.strand, test_cases.paf_test_record_ok.strand);
        assert_eq!(record.tname, test_cases.paf_test_record_ok.tname);
        assert_eq!(record.tlen, test_cases.paf_test_record_ok.tlen);
        assert_eq!(record.tstart, test_cases.paf_test_record_ok.tstart);
        assert_eq!(record.tend, test_cases.paf_test_record_ok.tend);
        assert_eq!(record.mlen, test_cases.paf_test_record_ok.blen);
        assert_eq!(record.mapq, test_cases.paf_test_record_ok.mapq);
    }
    #[test]
    fn paf_record_from_str_size_fail() {
        let test_cases = TestCases::new();
        let actual_error = PafRecord::from_str(test_cases.paf_test_str_size_fail).unwrap_err();
        let expected_error = PafRecordError::RecordSize();
        assert_eq!(actual_error, expected_error);
    }
    #[test]
    fn paf_record_query_aligned_length_ok() {
        let test_cases = TestCases::new();
        let record = PafRecord::from_str(test_cases.paf_test_str_ok).unwrap();
        let test_length = record.query_aligned_length().unwrap();
        assert_eq!(test_length, 4_u64);
    }
    #[test]
    fn paf_record_query_coverage_ok() {
        let test_cases = TestCases::new();
        let record = PafRecord::from_str(test_cases.paf_test_str_ok).unwrap();
        let query_coverage = record.query_coverage().unwrap();
        float_eq!(query_coverage, 1_f64, abs <= f64::EPSILON);
    }

    #[test]
    fn paf_record_query_coverage_zero_division_ok() {
        let test_cases = TestCases::new();
        let mut record = PafRecord::from_str(test_cases.paf_test_str_ok).unwrap();
        record.qlen = 0;
        let zero_length = record.query_coverage().unwrap();
        float_eq!(zero_length, 0_f64, abs <= f64::EPSILON);
    }

    /*
    ====================
          PafFile
    ====================
    */

    #[test]
    fn paf_parser_create_new_no_filter_ok() {
        let test_cases = TestCases::new();
        let paf_aln =
            PafAlignment::from(test_cases.paf_test_file_ok, None, 0_u64, 0_f64, 0_u8).unwrap();
        assert_eq!(paf_aln.records.len(), 5);
    }
    #[test]
    #[should_panic]
    fn paf_parser_create_new_record_size_fail() {
        let test_cases = TestCases::new();
        // PafAlignmentError does not implement PartialEq to assure standard Error type for parsing, let test fail instead
        let _ = PafAlignment::from(
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
            PafAlignment::from(test_cases.paf_test_file_ok, None, 0_u64, 0_f64, 30_u8).unwrap();
        assert_eq!(paf_aln.records.len(), 2);
    }
    #[test]
    fn paf_parser_create_new_filter_min_len_ok() {
        let test_cases = TestCases::new();
        let paf_aln =
            PafAlignment::from(test_cases.paf_test_file_ok, None, 50_u64, 0_f64, 0_u8).unwrap();
        assert_eq!(paf_aln.records.len(), 3);
    }
    #[test]
    fn paf_parser_create_new_filter_min_cov_ok() {
        let test_cases = TestCases::new();
        let paf_aln =
            PafAlignment::from(test_cases.paf_test_file_ok, None, 0_u64, 0.5, 0_u8).unwrap();
        assert_eq!(paf_aln.records.len(), 3);
    }

    #[test]
    fn paf_parser_create_new_fasta_input_provided_ok() {
        let test_cases = TestCases::new();
        let paf_aln = PafAlignment::from(
            test_cases.paf_test_file_ok,
            Some(test_cases.paf_test_fasta_ok),
            0_u64,
            0_f64,
            0_u8,
        )
        .unwrap();
        let mut expected: HashMap<String, u64> = HashMap::new();
        expected.insert("21172389_LCMV_L-segment_final".to_string(), 7194);
        expected.insert("21172389_LCMV_S-segment_final".to_string(), 3407);
        assert_eq!(paf_aln.seq_lengths, Some(expected));
    }

    #[test]
    fn paf_parser_create_new_fasta_input_not_provided_ok() {
        let test_cases = TestCases::new();
        let paf_aln =
            PafAlignment::from(test_cases.paf_test_file_ok, None, 0_u64, 0_f64, 0_u8).unwrap();
        assert_eq!(paf_aln.seq_lengths, None);
    }

    #[test]
    #[should_panic]
    fn paf_parser_create_new_fasta_input_not_exists_fail() {
        let test_cases = TestCases::new();
        let _ = PafAlignment::from(
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
            PafAlignment::from(test_cases.paf_test_file_ok, None, 0_u64, 0_f64, 0_u8).unwrap();
        let actual = paf_aln.target_intervals().unwrap();
        let actual_intervals_l = actual[0].1.intervals.clone();
        let actual_intervals_s = actual[1].1.intervals.clone();
        assert_eq!(actual_intervals_l, test_cases.paf_test_intervals_l_segment);
        assert_eq!(actual_intervals_s, test_cases.paf_test_intervals_s_segment);
    }

    #[test]
    fn paf_parser_compute_statistics_no_ref_ok() {
        let test_cases = TestCases::new();
        let paf_aln =
            PafAlignment::from(test_cases.paf_test_file_ok, None, 0_u64, 0_f64, 0_u8).unwrap();
        let actual_statistics = paf_aln.coverage_statistics(1).unwrap();
        assert_eq!(
            actual_statistics,
            test_cases.paf_test_coverage_statistics_no_ref
        );
    }

    #[test]
    fn paf_parser_compute_statistics_no_ref_no_tags_ok() {
        let test_cases = TestCases::new();
        let paf_aln =
            PafAlignment::from(test_cases.paf_test_file_ok, None, 0_u64, 0_f64, 0_u8).unwrap();
        let actual_statistics = paf_aln.coverage_statistics(0).unwrap();
        assert_eq!(
            actual_statistics,
            test_cases.paf_test_coverage_statistics_no_ref_no_tags
        );
    }

    #[test]
    // Weird edge case that probably doesn't occur ever, where there is a zero length sequence read from the reference sequence file
    fn paf_parser_compute_statistics_no_ref_seq_len_zero_ok() {
        let test_cases = TestCases::new();
        let mut paf_aln =
            PafAlignment::from(test_cases.paf_test_file_ok, None, 0_u64, 0_f64, 0_u8).unwrap();
        paf_aln.seq_lengths = Some(HashMap::from([
            ("21172389_LCMV_L-segment_final".to_string(), 0),
            ("21172389_LCMV_S-segment_final".to_string(), 0),
        ]));
        let actual_statistics = paf_aln.coverage_statistics(1).unwrap();
        assert_eq!(
            actual_statistics,
            test_cases.paf_test_coverage_statistics_no_ref
        );
    }

    #[test]
    fn paf_parser_compute_statistics_ref_ok() {
        let test_cases = TestCases::new();
        let paf_aln = PafAlignment::from(
            test_cases.paf_test_file_ok,
            Some(test_cases.paf_test_fasta_ok),
            0_u64,
            0_f64,
            0_u8,
        )
        .unwrap();
        let actual_statistics = paf_aln.coverage_statistics(1).unwrap();
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
    fn tabled_helper_format_tags_present_ok() {
        let expected = display_tags("test_tags");
        assert_eq!("test_tags".to_string(), expected);
    }

    #[test]
    fn tabled_helper_format_tags_absent_ok() {
        let expected = display_tags("");
        assert_eq!("-".to_string(), expected);
    }
}
