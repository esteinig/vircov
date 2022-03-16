use anyhow::Result;
use noodles::fasta;
use rust_lapper::{Interval, Lapper};
use std::collections::btree_map::Entry;
use std::collections::{HashMap, BTreeMap};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use thiserror::Error;
use itertools::Itertools;
use crate::covplot::CovPlot;
use crossterm::style::Color;

/*
========================
Custom error definitions
========================
*/

#[derive(Error, Debug)]
pub enum PafFileError {
    /// Indicates failure to read file
    #[error("failed to read file")]
    PafFileIOError(#[from] std::io::Error),
    /// Indicates failure to read a line into record
    #[error("failed to parse record")]
    PafRecordError(#[from] PafRecordError),
    /// Indicates a failure to get sequence length require for coverage plots
    #[error("failed to get sequence length for coverage plot")]
    PafCovPlotSeqLengthError(),
    /// Indicates failure with the coverage plot module
    #[error("failed to generate data for coverage plot")]
    PafCovPlotError(#[from] crate::covplot::CovPlotError),
}

#[derive(Error, Debug, PartialEq)]
pub enum PafRecordError {
    /// Indicates failure to parse a record due to insufficient fields
    #[error("record has too few fields")]
    PafRecordSizeError(),
    /// Indicates failure to parse an integer from a record string
    #[error("failed to parse an integer field")]
    PafRecordIntegerError(#[from] std::num::ParseIntError),
}

/*
=========================
PAF alignment file parser
=========================
*/

/// Paf struct
#[derive(Debug, Clone)]
pub struct PafFile {
    /// PafRecords parsed from file (PAF)
    records: Vec<PafRecord>,
    /// Reference sequence names and lengths from file (FASTA)
    seq_lengths: Option<HashMap<String, u64>>,
}

impl PafFile {
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
    ) -> Result<Self, PafFileError> {
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
    pub fn target_coverage_distribution(&self, verbosity: u64) -> Result<(), PafFileError> {
        for (target_name, targets) in self.target_intervals()? {
            // Bases of target sequence covered
            let target_cov_bp = targets.cov();
            
            let target_seq_len = match &self.seq_lengths {
                None => Some(&0),
                Some(seqs) => seqs.get(&target_name)
            };

            let target_seq_cov = match target_seq_len {
                None => 0_f64,
                Some(0) => 0_f64,
                Some(seq_len) => {
                    target_cov_bp as f64 / *seq_len as f64
                }
            };

            let target_seq_len_display = match target_seq_len {
                None => 0,
                Some(seq_len) => *seq_len
            };

            let mut merged_targets = targets.clone();
            merged_targets.merge_overlaps();

            let target_cov_n = merged_targets.intervals.len();

            let tags = match verbosity {
                1 => {
                    // Get the merged target interval data [start:stop:count]
                    let merged_target_intervals = merged_targets.depth().collect::<Vec<Interval<usize, usize>>>();
                    let merged_target_interval_strings = merged_target_intervals.iter().map(|i|{
                        let _start = i.start; let _stop = i.stop; let _count = targets.count(i.start, i.stop);
                        format!("{_start}:{_stop}:{_count}")
                    }).collect::<Vec<String>>();
                    merged_target_interval_strings.join(" ")
                },
                _ => {
                    "".to_string()
                }
            };
            

            // Get the number of unique reads in the alignments
            let reads: Vec<String> = targets.iter().map(|interval| interval.val.to_owned()).collect::<Vec<String>>();
            let unique_reads = reads.into_iter().unique().collect::<Vec<String>>();

            println!("{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}", &target_name, target_cov_n, unique_reads.len(), targets.len(), target_cov_bp, target_seq_len_display, target_seq_cov, tags);
        }

        Ok(())
    }
    /// Print the target coverage distributions to console in a reduced (approximate) text plot
    pub fn target_coverage_plots(&self, max_width: u64) -> Result<(), PafFileError>{


        println!();
        for (target_name, targets) in self.target_intervals()? {

            let target_seq_len = match &self.seq_lengths {
                None => Some(&0),
                Some(seqs) => seqs.get(&target_name)
            };

           let seq_length = match target_seq_len {
                None => return Err(PafFileError::PafCovPlotSeqLengthError()),
                Some(value) => *value
            };

            let covplot = CovPlot::new(targets, seq_length, max_width)?;

            covplot.to_console(target_name, seq_length, Color::Red)?;

        }

        Ok(())
    }

    /// Get target alignment interval lappers by target sequence
    fn target_intervals(&self) -> Result<Vec<(String, Lapper<usize, String>)>, PafFileError> {
        let mut target_intervals: BTreeMap<String, Vec<Interval<usize, String>>> = BTreeMap::new();
        for record in &self.records {
            match target_intervals.entry(record.tname.clone()) {
                Entry::Occupied(mut entry) => {
                    entry.get_mut().push(Interval {
                        start: record.tstart,
                        stop: record.tend,
                        val: record.qname.clone(),  // add the query read name here for reference, see if it affects performance
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
=================
PAF record struct
=================
*/

/// PAF record without tags
#[derive(Debug, Clone)]
#[allow(dead_code)]
struct PafRecord {
    /// Query sequence name.
    qname: String,
    /// Query sequence length.
    qlen: u64,
    /// Query start (0-based; BED-like; closed).
    qstart: usize,
    /// Query end (0-based; BED-like; open).
    qend: usize,
    /// Relative strand: "+" or "-".
    strand: String,
    /// Target sequence name.
    tname: String,
    /// Target sequence length.
    tlen: u64,
    /// Target start on original strand (0-based).
    tstart: usize,
    /// Target end on original strand (0-based).
    tend: usize,
    /// Number of matching bases in the mapping.
    mlen: u64,
    /// Alignment block length. Number of bases, including gaps, in the mapping.
    blen: u64,
    /// Mapping quality (0-255; 255 for missing).
    mapq: u8,
}

impl PafRecord {
    /// Populate a record from a record string, includes some error checks
    pub fn from_str(paf_str: String) -> Result<Self, PafRecordError> {
        let fields: Vec<&str> = paf_str.split("\t").collect();

        if fields.len() < 12 {
            return Err(PafRecordError::PafRecordSizeError());
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
        let expected_error = PafRecordError::PafRecordSizeError();
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
        PafAlignment
    ====================
    */

    #[test]
    fn paf_parser_create_new_no_filter_ok() {
        let test_cases = TestCases::new();
        let paf_aln = PafFile::from(test_cases.paf_test_file_ok, None, 0_u64, 0_f64, 0_u8).unwrap();
        assert_eq!(paf_aln.records.len(), 5);
    }
    #[test]
    #[should_panic]
    fn paf_parser_create_new_record_size_fail() {
        let test_cases = TestCases::new();
        // PafAlignmentError does not implement PartialEq to assure standard Error type for parsing, let test fail instead
        let _ = PafFile::from(
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
            PafFile::from(test_cases.paf_test_file_ok, None, 0_u64, 0_f64, 30_u8).unwrap();
        println!("{:?}", paf_aln.records);
        assert_eq!(paf_aln.records.len(), 2);
    }
    #[test]
    fn paf_parser_create_new_filter_min_len_ok() {
        let test_cases = TestCases::new();
        let paf_aln =
            PafFile::from(test_cases.paf_test_file_ok, None, 50_u64, 0_f64, 0_u8).unwrap();
        assert_eq!(paf_aln.records.len(), 3);
    }
    #[test]
    fn paf_parser_create_new_filter_min_cov_ok() {
        let test_cases = TestCases::new();
        let paf_aln = PafFile::from(test_cases.paf_test_file_ok, None, 0_u64, 0.5, 0_u8).unwrap();
        assert_eq!(paf_aln.records.len(), 3);
    }
}
