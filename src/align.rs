use crate::covplot::CovPlot;
use crate::utils::get_sanitized_fasta_writer;
use crate::utils::{get_grouped_segments, get_segment_selections};
use anyhow::Result;
use crossterm::style::Color;
use itertools::Itertools;
use noodles::fasta::{Reader as FastaReader, Record as FastaRecord};
use ordered_float::OrderedFloat;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::{bam, bam::record::Cigar, bam::HeaderView, bam::Read};
use rust_lapper::{Interval, Lapper};
use std::collections::btree_map::Entry;
use std::collections::{BTreeMap, HashMap};
use std::convert::TryInto;
use std::fs::File;
use std::io::Write;
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
    #[error("failed to parse a valid integer from record")]
    FloatError(#[from] std::num::ParseFloatError),
    /// Indicates failure to conduct grouping because no reference sequences were parsed
    #[error("failed to group outputs due to missing reference sequences")]
    GroupSequenceError,
    /// Indicates failure to plot coverage when data is grouped
    #[error("coverage plots are not enabled when grouping output")]
    GroupCovPlotError,
    /// Indicates failure to infer or identify file format from explicit option
    #[error("failed to parse a valid input format")]
    InputFormatError,
    /// Indicates failure when no grouping options are provided and selection is specified
    #[error("failed to use the group selection options becuase no grouping is specified")]
    GroupSelectSplitError,
    /// Indicates failure when no group select by is provided (should not occurr)
    #[error("failed to group select by reads or coverage")]
    GroupSelectByError,
    /// Indicates failure when no group select by is provided (should not occurr)
    #[error("failed to provide a negative segment field")]
    SegmentFieldNaNError,
    /// Indicates failure when no best value from the grouped coverage fields could be selected
    #[error("failed to select a best reference sequence")]
    GroupSelectReference,
    /// Indicates failure when no reference name could be selected from the grouped identifier
    #[error("failed to extract a reference name from the grouped identifier")]
    GroupSelectReferenceName,
    /// Indicates failure when no coverage value could be extracted from a group of coverage fields
    #[error("failed to extract the highest coverage value from the grouped fields")]
    GroupSelectCoverage,
    /// Indicates failure when zero-value reference sequences should be included, but no reference sequences were provided
    #[error("no reference sequences found, zero-value records cannot be included")]
    ZeroReferenceSequences,
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
pub struct CoverageFields {
    /// Name of the target sequence
    pub name: String,
    /// Number of non-overlapping alignment regions
    pub regions: u64,
    /// Number of unique reads aligned
    pub reads: u64,
    /// Number of alignments
    pub alignments: u64,
    /// Number of bases covered by alignments
    pub bases: u64,
    /// Length of the target sequence
    pub length: u64,
    /// Fractional coverage of the alignments
    pub coverage: f64,
    /// Descriptions of the reference sequence headers
    pub description: String,
    /// Tags for alignment regions in start:stop:aln format
    pub tags: String,
    /// Unique read identifiers for grouped operations
    pub unique_reads: Vec<String>,
}

impl CoverageFields {
    fn as_tag(&self) -> String {
        format!(
            "{:}|{:}|{:}|{:}|{:}|{:}|{:}|{:}",
            self.name,
            self.regions,
            self.reads,
            self.alignments,
            self.bases,
            self.length,
            self.coverage,
            self.description
        )
    }
    fn from_tag(tag: &str) -> Result<Self, ReadAlignmentError> {
        let fields: Vec<_> = tag.split('|').into_iter().collect();

        Ok(Self {
            name: fields[0].to_string(),
            regions: fields[1].parse::<u64>()?,
            reads: fields[2].parse::<u64>()?,
            alignments: fields[3].parse::<u64>()?,
            bases: fields[4].parse::<u64>()?,
            length: fields[5].parse::<u64>()?,
            coverage: fields[6].parse::<f64>()?,
            description: fields[7].to_string(),
            tags: "".to_string(),
            unique_reads: Vec::new(),
        })
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
    #[header("Reads")]
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

fn parse_exclude(exclude: Option<PathBuf>) -> Result<Option<Vec<String>>, ReadAlignmentError> {
    match exclude {
        Some(_file) => {
            let reader = File::open(_file).map(BufReader::new)?;

            let mut exclude_vec: Vec<String> = Vec::new();
            for line in reader.lines() {
                let content: String = line?;
                let content: &str = content.trim();

                if content.starts_with('|') || content.is_empty() {
                } else {
                    exclude_vec.push(content.to_lowercase())
                }
            }

            Ok(Some(exclude_vec))
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
    /// Reference sequence description strings to exclude
    pub target_exclude: Option<Vec<String>>,
}

impl ReadAlignment {
    pub fn new(
        // Reference sequence fasta file
        fasta: &Option<PathBuf>,
        // Reference sequence exclude file
        exclude: &Option<PathBuf>,
    ) -> Result<Self, ReadAlignmentError> {
        Ok(Self {
            target_intervals: Vec::new(),
            target_sequences: parse_fasta(fasta.clone())?,
            target_exclude: parse_exclude(exclude.clone())?,
        })
    }
    // Parse alignment by inferring format from file extension
    pub fn read(
        &mut self,
        // Path to alignment file [PAF or "-"]
        path: PathBuf,
        // Minimum query alignment length
        min_qaln_len: u64,
        // Minimum query coverage
        min_qaln_cov: f64,
        // Minimum mapping quality
        min_mapq: u8,
        // Explicit alignment format
        alignment_format: Option<String>,
    ) -> Result<&Self, ReadAlignmentError> {
        match alignment_format {
            Some(format) => match format.as_str() {
                "bam" => self.read_bam(path, min_qaln_len, min_qaln_cov, min_mapq),
                "paf" => self.read_paf(path, min_qaln_len, min_qaln_cov, min_mapq),
                _ => Err(ReadAlignmentError::InputFormatError),
            },
            None => match path.extension().map(|s| s.to_str()) {
                Some(Some("paf")) => self.read_paf(path, min_qaln_len, min_qaln_cov, min_mapq),
                Some(Some("bam") | Some("sam") | Some("cram")) => {
                    self.read_bam(path, min_qaln_len, min_qaln_cov, min_mapq)
                }
                _ => Err(ReadAlignmentError::InputFormatError),
            },
        }
    }
    // Parse alignments from file
    pub fn read_paf(
        &mut self,
        // Path to alignment file [PAF or "-"]
        path: PathBuf,
        // Minimum query alignment length
        min_qaln_len: u64,
        // Minimum query coverage
        min_qaln_cov: f64,
        // Minimum mapping quality
        min_mapq: u8,
    ) -> Result<&Self, ReadAlignmentError> {
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

        self.target_intervals = target_lappers;

        Ok(self)
    }
    // Parse alignments from file
    pub fn read_bam(
        &mut self,
        // Path to alignment file [SAM/BAM/CRAM or "-"]
        path: PathBuf,
        // Minimum query alignment length
        min_qaln_len: u64,
        // Minimum query coverage
        min_qaln_cov: f64,
        // Minimum mapping quality
        min_mapq: u8,
    ) -> Result<&Self, ReadAlignmentError> {
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

        self.target_intervals = target_lappers;

        Ok(self)
    }
    #[allow(clippy::too_many_arguments)]
    /// Compute coverage distribution by target sequence
    pub fn coverage_statistics(
        &self,
        regions: u64,
        seq_len: u64,
        coverage: f64,
        regions_coverage: Option<f64>,
        reads: u64,
        aligned: u64,
        group_by: &Option<String>,
        verbosity: u64,
        zero: bool,
    ) -> Result<Vec<CoverageFields>, ReadAlignmentError> {
        let mut coverage_fields: Vec<CoverageFields> = Vec::new();
        let mut included_references: Vec<String> = Vec::new();

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

            let target_description_decap = target_description.to_lowercase(); // exclude terms are parsed as lowercase (case-insensitive matching)

            let exclude_target_sequence = match &self.target_exclude {
                None => false,
                Some(exclude_list) => exclude_list
                    .iter()
                    .any(|x| target_description_decap.contains(x)),
            };

            let unique_read_ids: Vec<String> = targets
                .iter()
                .map(|interval| interval.val.to_string())
                .unique()
                .collect::<Vec<String>>();

            let unique_reads_n = unique_read_ids.len() as u64;

            let region_filter_passed: bool = match regions_coverage {
                Some(reg_cov_threshold) => {
                    // Apply the region filter only if the reference sequence
                    // coverage is below a coverage threshold - only apply this if
                    // the threshold is bigger than zero so that pipelines can still
                    // apply the option with a 0 value
                    if (reg_cov_threshold > 0.) && (target_cov < reg_cov_threshold) {
                        target_cov_n >= regions
                    } else {
                        // otherwise, if coverage is above the threshold, do
                        // not apply the region filter
                        true
                    }
                }
                _ => target_cov_n >= regions,
            };

            let reads_aligned = targets.len() as u64;

            if target_len >= seq_len
                && reads_aligned >= aligned
                && region_filter_passed
                && target_cov >= coverage
                && unique_reads_n >= reads
                && !exclude_target_sequence
            {
                coverage_fields.push(CoverageFields {
                    name: target_name.to_owned(),
                    regions: target_cov_n,
                    reads: unique_reads_n,
                    alignments: reads_aligned,
                    bases: target_cov_bp,
                    length: target_len,
                    coverage: target_cov,
                    description: target_description,
                    unique_reads: unique_read_ids,
                    tags: target_tags,
                });
                // For zero count checks below
                included_references.push(target_name.to_owned())
            }
        }

        if zero {
            // This will add all reference sequences as null fields
            // to the coverage field data - beware that this will
            // include these sequences in grouped output
            match &self.target_sequences {
                Some(ref_seqs) => {
                    for (ref_seq, record) in ref_seqs.iter() {
                        if !included_references.contains(ref_seq) {
                            let descr = match record.description() {
                                None => "-".to_string(),
                                Some(descr) => descr.to_owned(),
                            };
                            coverage_fields.push(CoverageFields {
                                name: ref_seq.to_owned(),
                                regions: 0,
                                reads: 0,
                                alignments: 0,
                                bases: 0,
                                length: record.sequence().len() as u64,
                                coverage: 0.0,
                                description: descr,
                                unique_reads: Vec::new(),
                                tags: "-".to_string(),
                            });
                        }
                    }
                }
                None => return Err(ReadAlignmentError::ZeroReferenceSequences),
            }
        }

        Ok(coverage_fields)
    }
    #[cfg(not(tarpaulin_include))]
    #[allow(clippy::too_many_arguments)]
    /// Print the coverage statistics to console, alternatively in pretty table format
    pub fn to_output(
        &self,
        coverage_fields: &mut Vec<CoverageFields>, // If grouped, these are grouped fields
        ungrouped_fields: Option<Vec<CoverageFields>>, // if grouped, ungrouped fields should be passed for select output writing
        table: bool,
        header: bool,
        group_sep: String,
        read_ids: Option<PathBuf>,
        read_ids_split: Option<PathBuf>,
        group_select_by: Option<String>,
        group_select_split: Option<PathBuf>,
        group_select_order: bool,
        group_select_data: Option<PathBuf>,
        segment_field: Option<String>,
        segment_field_nan: Option<String>,
    ) -> Result<(), ReadAlignmentError> {
        // If grouped and selected we skip as we want the  ungrouped summary
        // of the selected group members as primary output and the grouped data
        // as file output
        match group_select_split {
            None => self.write_table(coverage_fields, table, header, None)?,
            Some(_) => {}
        }

        // Non grouped data unique read outputs + grouped

        if let Some(path) = read_ids {
            let mut file_handle = File::create(&path)?;

            let all_unique_read_ids: Vec<String> = coverage_fields
                .iter()
                .flat_map(|field| field.unique_reads.to_owned()) // in case it is grouped
                .unique()
                .collect();

            for read_id in all_unique_read_ids {
                writeln!(file_handle, "{}", &read_id)?;
            }
        }

        if let Some(path) = read_ids_split {
            std::fs::create_dir_all(&path)?;
            for field in coverage_fields.iter_mut() {
                let sanitized_name = field.name.replace(' ', "_");
                let sanitized_name = sanitized_name.trim_matches(';'); // Virosaurus sanitize remainign header separator on seq id (weird format)

                let file_path = path.join(sanitized_name).with_extension("txt");
                let mut file_handle = File::create(file_path.as_path())?;
                for read_id in field.unique_reads.iter() {
                    writeln!(file_handle, "{}", &read_id)?
                }
            }
        }

        if let Some(path) = group_select_split.clone() {
            match &self.target_sequences {
                Some(ref_seqs) => {
                    std::fs::create_dir_all(&path)?;

                    // Flipped comparison function, highest to lowest sort
                    let coverage_fields = match group_select_order {
                        true => match group_select_by.clone() {
                            Some(value) => match value.as_str() {
                                "reads" => {
                                    coverage_fields.sort_by(|a, b| b.reads.cmp(&a.reads));
                                    coverage_fields
                                }
                                "coverage" => {
                                    coverage_fields.sort_by(|a, b| {
                                        b.coverage.partial_cmp(&a.coverage).unwrap()
                                    });
                                    coverage_fields
                                }
                                _ => return Err(ReadAlignmentError::GroupSelectByError),
                            },
                            None => coverage_fields,
                        },
                        false => coverage_fields,
                    };

                    // Store the reference identifiers of the selected grouped outputs
                    // to write a summary of the ungrouped field data for these
                    let mut selected_output_identifiers = Vec::new();

                    for (cov_field_idx, cov_field) in coverage_fields.iter().enumerate() {
                        // In the grouped data the fields are itself contained in the tag
                        // the following will produce a CoverageField for each tag using
                        // the `from_tag` parsing method
                        let tags = cov_field
                            .tags
                            .split("~~~")
                            .into_iter()
                            .filter_map(|x| match x {
                                "" => None,
                                _ => Some(CoverageFields::from_tag(x).unwrap()),
                            })
                            .collect::<Vec<CoverageFields>>();

                        // In grouped coverage field structs, the description contains a
                        // concatenated description of the grouped reference sequences;
                        // may need to modify later and add it to the group tags (tricky be)
                        let segmented = match segment_field.clone() {
                            Some(seg_field) => match segment_field_nan.clone() {
                                Some(seg_field_nan) => {
                                    let segment_tag_count =
                                        cov_field.description.matches(&seg_field).count();
                                    let segment_tag_na_count =
                                        cov_field.description.matches(&seg_field_nan).count();
                                    (segment_tag_count > 0)
                                        && (segment_tag_na_count != segment_tag_count)
                                }
                                None => return Err(ReadAlignmentError::SegmentFieldNaNError),
                            },
                            None => false,
                        };

                        // VIROSAURUS GENE IDENTIFIER - RETAIN FOR OPTION LATER

                        // let gene_split = (cov_field.tags.matches("GENE").count() > 0)
                        //     || (cov_field.name.matches("GENE").count() > 0);

                        // // Output all grouped reference sequences if segments are present (fix) or if the
                        // // GENE identifier is present (e.g. Virosaurus) - this will put all sequences into
                        // // a multi-FASTA

                        // if gene_split {
                        //     let name = &tags[0].name;
                        //     let mut writer = get_sanitized_fasta_writer(name, &path).expect(
                        //         "Could not get sanitized FASTA writer for GENE multi-FASTA",
                        //     );
                        //     for tag in &tags {
                        //         let seq = &ref_seqs[&tag.name];
                        //         writer.write_record(seq)?
                        //     }
                        // } else {

                        // This is the output name using the grouping variable (e.g. taxonomic identifier)
                        let cov_field_name = match cov_field.name.split_whitespace().next() {
                            Some(str) => str,
                            None => return Err(ReadAlignmentError::GroupSelectReferenceName),
                        };

                        // Filter reference sequences (contained in tags) by highest number
                        // of mapped reads or coverage
                        if segmented {
                            // If the grouped sequences are segmented, group unique segments (from identifier)
                            // and select the best segment of each group based on the group_select_by metric

                            let grouped_segments = get_grouped_segments(
                                tags,
                                segment_field.clone(),
                                group_sep.clone(),
                            )?;

                            let grouped_segment_selections =
                                get_segment_selections(grouped_segments, group_select_by.clone())?;

                            let seq_name = match group_select_order {
                                true => {
                                    format!("{:0>2}-{:}", cov_field_idx, cov_field_name)
                                }
                                false => cov_field.name.to_string(),
                            };

                            let mut writer = get_sanitized_fasta_writer(&seq_name, &path)
                                .expect("Could not get sanitized writer for segment multi-FASTA");

                            for (_, segment_field) in grouped_segment_selections {
                                let seq = &ref_seqs[&segment_field.name];
                                writer.write_record(seq)?;
                                selected_output_identifiers.push(segment_field.name);
                            }
                        } else {
                            // If not segmented, simply select the best reference sequence using the
                            // group_select_by metric and select the sequence from the header fields
                            // to write to FASTA
                            let max = match group_select_by.clone() {
                                Some(value) => match value.as_str() {
                                    "reads" => tags.iter().max_by_key(|x| x.reads),
                                    "coverage" => {
                                        tags.iter().max_by_key(|x| OrderedFloat(x.coverage))
                                    }
                                    _ => return Err(ReadAlignmentError::GroupSelectByError),
                                },
                                None => return Err(ReadAlignmentError::GroupSelectByError),
                            };

                            match max {
                                Some(field) => {
                                    let seq = &ref_seqs[&field.name];

                                    let seq_name = match group_select_order {
                                        true => {
                                            format!("{:0>2}-{:}", cov_field_idx, cov_field_name)
                                        }
                                        false => cov_field_name.to_string(),
                                    };

                                    let mut writer = get_sanitized_fasta_writer(&seq_name, &path)
                                        .expect("Could not get sanitized writer for single FASTA");
                                    writer.write_record(seq)?;
                                    selected_output_identifiers.push(field.name.clone());
                                }
                                _ => return Err(ReadAlignmentError::GroupSelectReference),
                            }
                        }
                    }

                    // We should now write a table of the selected identifiers from each group, so that the results are not group results at this stage
                    // (also written but these have the group summary values, not the individual alignment summaries of the selected sequences)

                    if let (Some(ungrouped_fields), Some(_)) =
                        (&ungrouped_fields, &group_select_split)
                    {
                        let selected_ungrouped_fields: Vec<CoverageFields> = ungrouped_fields
                            .iter()
                            .filter(|field| selected_output_identifiers.contains(&field.name))
                            .cloned()
                            .collect();
                        self.write_table(&selected_ungrouped_fields, table, header, None)?;
                        self.write_table(coverage_fields, table, header, group_select_data)?;
                    };
                }
                _ => return Err(ReadAlignmentError::GroupSequenceError),
            }
        }

        Ok(())
    }
    pub fn write_table(
        &self,
        coverage_fields: &[CoverageFields],
        table: bool,
        header: bool,
        file: Option<PathBuf>,
    ) -> Result<(), ReadAlignmentError> {
        match table {
            true => {
                let table_fields = coverage_fields
                    .iter()
                    .map(CoverageTableFields::from)
                    .collect::<Vec<CoverageTableFields>>();
                let _table = Table::new(table_fields)
                    .with(Modify::new(Column(7..)).with(MaxWidth::wrapping(32)))
                    .with(Style::modern());
                println!("{}", _table);
                Ok(())
            }
            false => {
                if header {
                    println!("sequence\tregions\treads\talignments\tbases\tlength\tcoverage\tdescription\ttags")
                }
                match file {
                    Some(file) => {
                        let mut file_handle = File::create(file)?;
                        for cov_fields in coverage_fields.iter() {
                            writeln!(
                                file_handle,
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
                            )?;
                        }
                        Ok(())
                    }
                    None => {
                        for cov_fields in coverage_fields.iter() {
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
                        Ok(())
                    }
                }
            }
        }
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

    #[allow(clippy::too_many_arguments)]
    pub fn group_output(
        &self,
        coverage_fields: &[CoverageFields],
        grouped_regions: u64,
        grouped_coverage: f64,
        grouped_alignments: u64,
        grouped_reads: u64,
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
                            entry.insert(vec![cov_field]);
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

            let mut ureads: Vec<String> = Vec::new();

            for field in &fields {
                grouped_fields.regions += field.regions;
                grouped_fields.alignments += field.alignments;
                grouped_fields.coverage += field.coverage;
                grouped_fields.bases += field.bases;

                if !grouped_fields.description.contains(&field.description) {
                    grouped_fields.description.push_str(&field.description);
                    grouped_fields.description.push(' ');
                }

                match field.tags.as_str() {
                    "-" => {
                        if !grouped_fields.tags.contains('-') {
                            grouped_fields.tags.push('-');
                        }
                    }
                    _ => {
                        // Separate tags by a weird delimiter (~~~) to include
                        // white space separated header descriptions in the tags
                        grouped_fields.tags.push_str(&field.as_tag());
                        grouped_fields.tags.push_str("~~~");
                    }
                }
                for read in &field.unique_reads {
                    ureads.push(read.to_string())
                }
            }

            // Trim the last grouped tag separator introduced by the field loop above (~~~)
            grouped_fields.tags = grouped_fields.tags.trim_end_matches("~~~").to_string();

            grouped_fields.bases /= fields.len() as u64;

            grouped_fields.coverage = match fields.iter().max_by_key(|x| OrderedFloat(x.coverage)) {
                Some(fie) => fie.coverage,
                None => return Err(ReadAlignmentError::GroupSelectCoverage),
            };

            let unique_reads_grouped: Vec<String> =
                ureads.iter().unique().map(|x| x.to_string()).collect();

            grouped_fields.reads = unique_reads_grouped.len() as u64;
            grouped_fields.unique_reads = unique_reads_grouped;

            if grouped_fields.regions >= grouped_regions
                && grouped_fields.coverage >= grouped_coverage
                && grouped_fields.alignments >= grouped_alignments
                && grouped_fields.reads >= grouped_reads
            {
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
/// PAF considers insertions but not deletions when computing
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

    #[allow(dead_code)]
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
}
