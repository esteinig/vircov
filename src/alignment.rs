use crate::covplot::CovPlot;
use crate::error::VircovError;
use crate::vircov::{AlignerConfig, Annotation, AnnotationConfig, FilterConfig, ReferenceConfig};


use anyhow::Result;
use crossterm::style::Color;
use itertools::Itertools;
use noodles::fasta::{Reader as FastaReader, Record as FastaRecord};
use ordered_float::OrderedFloat;


use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::{bam, bam::record::Cigar, bam::HeaderView, bam::Read};
use std::process::{ChildStdout, Command, Output, Stdio};
use std::str::from_utf8;
use std::convert::TryInto;

use rust_lapper::{Interval, Lapper};
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use serde::ser::SerializeSeq;
use serde::de::{Visitor, SeqAccess};
use std::collections::btree_map::Entry;
use std::collections::{BTreeMap, HashMap, HashSet};

use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use tabled::{settings::{Width, Style, object::Columns}, Table, Tabled};

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, clap::ValueEnum)]
pub enum SelectHighest {
    Reads,
    Coverage, 
    Alignments
}

/*
=====================
Table display helpers
=====================
*/

pub fn display_coverage(cov: &f64) -> String {
    format!("{:.2}", cov)
}

trait IntervalTag {
    fn as_tag(&self) -> String;
    fn from_tag(tag: &str) -> Result<Self, String>
    where
        Self: Sized;
}

type AlignmentInterval = Interval<usize, usize>;

impl IntervalTag for AlignmentInterval {
    fn as_tag(&self) -> String {
        format!("{}:{}:{}", self.start, self.stop, self.val)
    }
    fn from_tag(tag: &str) -> Result<Self, String> {
        let parts: Vec<&str> = tag.split(':').collect();
        if parts.len() != 3 {
            return Err(format!("Invalid interval tag: {}", tag));
        }

        let start = parts[0].parse().map_err(|_| "Invalid start value".to_string())?;
        let stop = parts[1].parse().map_err(|_| "Invalid stop value".to_string())?;
        let val = parts[2].parse().map_err(|_| "Invalid value".to_string())?;

        Ok(AlignmentInterval { start, stop, val })
    }
}


fn serialize_intervals<S>(intervals: &Vec<AlignmentInterval>, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    let mut seq = serializer.serialize_seq(Some(intervals.len()))?;
    for interval in intervals {
        seq.serialize_element(&interval.as_tag())?;
    }
    seq.end()
}


fn deserialize_intervals<'de, D>(deserializer: D) -> Result<Vec<AlignmentInterval>, D::Error>
where
    D: Deserializer<'de>,
{
    struct IntervalVisitor;

    impl<'de> Visitor<'de> for IntervalVisitor {
        type Value = Vec<AlignmentInterval>;

        fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
            formatter.write_str("a sequence of strings representing intervals")
        }

        fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
        where
            A: SeqAccess<'de>,
        {
            let mut intervals = Vec::new();

            while let Some(value) = seq.next_element::<String>()? {
                let interval = AlignmentInterval::from_tag(&value).map_err(serde::de::Error::custom)?;
                intervals.push(interval);
            }

            Ok(intervals)
        }
    }

    deserializer.deserialize_seq(IntervalVisitor)
}

/// A struct for computed output fields
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Coverage {
    /// Name of the target sequence
    pub reference: String,
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
    pub group: Option<String>,
    /// Descriptions of the reference sequence headers
    pub segment: Option<String>,
    /// Descriptions of the reference sequence headers
    pub name: Option<String>,
    /// Descriptions of the reference sequence headers
    pub description: String,
    #[serde(skip)]
    /// Tags for alignment regions in string format
    pub intervals: Vec<AlignmentInterval>,
    #[serde(skip)]
    /// Unique read identifiers for grouped operations
    pub read_id: Vec<String>,
}

impl std::fmt::Display for Coverage {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} | {}", self.reference, self.description)?;
        Ok(())
    }
}


/// A struct for computed output fields
#[derive(Debug, Clone, PartialEq)]
pub struct GroupedCoverage {
    /// Group of the target sequence
    pub group: String,
    /// Count of the 
    pub count: usize,
    /// Number of non-overlapping alignment regions
    pub total_regions: u64,
    /// Number of unique reads aligned to this group
    pub total_reads: u64,
    /// Number of alignments
    pub total_alignments: u64,
    /// Number of bases covered by alignments
    pub total_bases: u64,
    /// Fractional coverage of the alignments
    pub mean_coverage: f64,
    /// Fractional coverage of the alignments
    pub max_coverage: f64,
    /// Coverage tags for alignment regions 
    pub coverage: Vec<Coverage>,
    /// Unique read identifiers for grouped operations
    pub read_id: Vec<String>,
}

impl std::fmt::Display for GroupedCoverage {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "{} (n = {})", self.group, self.count)?;
        for cov in &self.coverage {
            writeln!(f, "{}", cov)?;
        }
        Ok(())
    }
}

impl GroupedCoverage {
    pub fn from_coverage(group: &str, coverage: Vec<&Coverage>) -> Result<Self, VircovError> {

        let mut grouped_coverage = Self {
            group: group.to_string(),
            count: coverage.len(),
            total_regions: 0,
            total_reads: 0,
            total_alignments: 0,
            total_bases: 0,
            mean_coverage: 0.,
            max_coverage: 0.,
            coverage: Vec::new(),
            read_id: Vec::new(),
        };
        let num_cov = coverage.len();
        let mut ureads: Vec<String> = Vec::new();


        grouped_coverage.max_coverage = match coverage.iter().max_by_key(|x| OrderedFloat(x.coverage)) {
            Some(fie) => fie.coverage,
            None => return Err(VircovError::GroupSelectCoverage),
        };

        for field in coverage {
            grouped_coverage.total_regions += field.regions;
            grouped_coverage.total_alignments += field.alignments;
            grouped_coverage.mean_coverage += field.coverage;
            grouped_coverage.total_bases += field.bases;

            grouped_coverage.coverage.push(field.clone());

            for read in &field.read_id {
                ureads.push(read.to_string())
            }
        }

        grouped_coverage.mean_coverage /= num_cov as f64;

        let unique_reads_grouped: Vec<String> = ureads.iter().unique().map(|x| x.to_string()).collect();

        grouped_coverage.total_reads = unique_reads_grouped.len() as u64;
        grouped_coverage.read_id = unique_reads_grouped;

        Ok(grouped_coverage)
    }
}

/// A struct for computed output fields for display
#[derive(Debug, Clone, PartialEq, Tabled)]
pub struct CoverageTableFields {
    #[tabled(rename="Reference")]
    /// Name of the target sequence
    name: String,
    #[tabled(rename="Regions")]
    /// Number of non-overlapping alignment regions
    regions: u64,
    #[tabled(rename="Reads")]
    /// Number of unique reads aligned
    reads: u64,
    #[tabled(rename="Alignments")]
    /// Number of alignments
    alignments: u64,
    #[tabled(rename="Bases")]
    /// Number of bases covered by alignments
    bases: u64,
    #[tabled(rename="Length")]
    /// Length of the target sequence
    length: u64,
    #[tabled(rename="Coverage")]
    #[tabled(display_with = "display_coverage")]
    /// Fractional coverage of the alignments
    coverage: f64,
    /// Descriptions of the reference sequence headers
    #[tabled(rename="Description")]
    description: String,
    /// Tags for alignment regions in start:stop:aligned format
    #[tabled(rename="Intervals")]
    intervals: String,
}

impl CoverageTableFields {
    /// Takes a coverage field and subsets it to fields
    /// suitable for tabled output
    pub fn from(coverage_fields: &Coverage) -> Self {
        Self {
            name: coverage_fields.reference.clone(),
            regions: coverage_fields.regions,
            reads: coverage_fields.reads,
            alignments: coverage_fields.alignments,
            bases: coverage_fields.bases,
            length: coverage_fields.length,
            coverage: coverage_fields.coverage*100.,
            description: coverage_fields.description.clone(),
            intervals: coverage_fields.intervals.iter().map(|interval| interval.as_tag()).join(" "),
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
pub fn parse_reference_fasta(
    fasta: Option<PathBuf>,
) -> Result<Option<HashMap<String, FastaRecord>>, VircovError> {
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

pub fn parse_exclude_file(exclude: Option<PathBuf>) -> Result<Option<Vec<String>>, VircovError> {
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

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, clap::ValueEnum)]
pub enum Preset {
    LrHq,
    Splice,
    SpliceHq,
    Asm,        // minigraph + minimap
    Asm5,
    Asm10,
    Asm20,
    Sr,         // minigraph + minimap
    Lr,         // minigraph 
    MapPb,
    MapHifi,
    MapOnt,
    AvaPb,
    AvaOnt,
}
impl std::fmt::Display for Preset {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Preset::LrHq => write!(f, "lr:hq"),
            Preset::Splice => write!(f, "splice"),
            Preset::SpliceHq => write!(f, "splice:hq"),
            Preset::Asm => write!(f, "asm"),
            Preset::Asm5 => write!(f, "asm5"),
            Preset::Asm10 => write!(f, "asm10"),
            Preset::Asm20 => write!(f, "asm20"),
            Preset::Sr => write!(f, "sr"),
            Preset::Lr => write!(f, "lr"),
            Preset::MapPb => write!(f, "map-pb"),
            Preset::MapHifi => write!(f, "map-hifi"),
            Preset::MapOnt => write!(f, "map-ont"),
            Preset::AvaPb => write!(f, "ava-pb"),
            Preset::AvaOnt => write!(f, "ava-ont")
        }
    }
}



#[derive(Serialize, Deserialize, Clone, Debug, clap::ValueEnum)]
pub enum AlignmentFormat {
    Sam,
    Bam,
    Cram,
    Paf,
    Gaf
}


/// Enum representing the available aligners.
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, clap::ValueEnum)]
pub enum Aligner {
    #[serde(rename="bowtie2")]
    Bowtie2,
    #[serde(rename="minimap2")]
    Minimap2,
    #[serde(rename="strobealign")]
    Strobealign,
    #[serde(rename="minimap2-rs")]
    #[cfg(feature = "mm2")]
    Minimap2Rs
}
impl Aligner {
    // Used for identification of pre-built-indices
    pub fn short_name(&self) -> &str {
        match self {
            Aligner::Bowtie2 => "bt",
            Aligner::Minimap2 => "mm",
            Aligner::Strobealign => "st",
            #[cfg(feature = "mm2")]
            Aligner::Minimap2Rs => "mm"
        }
    }
}
impl std::fmt::Display for Aligner {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Aligner::Bowtie2 => write!(f, "bowtie2"),
            Aligner::Minimap2 => write!(f, "minimap2"),
            Aligner::Strobealign => write!(f, "strobealign"),
            #[cfg(feature = "mm2")]
            Aligner::Minimap2Rs => write!(f, "minimap2-rs"),
        }
    }
}


pub struct VircovAligner {
    pub config: AlignerConfig,
    pub reference: ReferenceConfig,
    pub filter: FilterConfig
}
impl VircovAligner {
    pub fn from(config: &AlignerConfig, reference: &ReferenceConfig, filter: &FilterConfig) -> Self {
        Self {
            config: config.clone(), reference: reference.clone(), filter: filter.clone()
        }
    }
    pub fn check_aligner_dependency(&self, aligner: &Aligner) -> Result<(), VircovError> {
        let command = match aligner {
            Aligner::Minimap2 => "minimap2 --version",
            Aligner::Bowtie2 => "bowtie2 --version",
            Aligner::Strobealign => "strobealign --version",
            #[cfg(feature = "mm2")]
            Aligner::Minimap2Rs => return Ok(())
        };
        self.run_version_command(command).map_err(|_| VircovError::AlignerDependencyMissing(aligner.clone()))?;
        Ok(())
    }
    pub fn run_aligner(&self) -> Result<Option<ReadAlignment>, VircovError> {
        match self.config.aligner {
            Aligner::Minimap2 => self.run_minimap2(),
            Aligner::Bowtie2 => self.run_bowtie2(),
            Aligner::Strobealign => self.run_strobealign(),
            #[cfg(feature = "mm2")]
            Aligner::Minimap2Rs => self.run_minimap2_rs()?,
        }
    }
    fn run_version_command(&self, command: &str) -> Result<Output, VircovError> {
        let output = Command::new("sh")
            .arg("-c")
            .arg(command)
            .output()
            .map_err(|e| VircovError::CommandExecutionFailed(command.to_string(), e.to_string()))?;

        if !output.status.success() {
            return Err(VircovError::CommandFailed(command.to_string(), output.status.code().unwrap_or(-1)));
        }

        Ok(output)
    }
    fn run_minimap2(&self) -> Result<Option<ReadAlignment>, VircovError> {
        let aligner_args = self.config.args.as_deref().unwrap_or("");
        let aligner_preset = self.config.preset.clone().ok_or(VircovError::MissingMinimap2Preset)?;

        
        let output = match &self.config.output {
            Some(output) => format!("> {}", output.display()),
            None => String::from("")
        };
        
        let secondary_arg = if self.config.secondary {
            "--secondary=yes" // TODO: documentation default behaviour of Bowtie2
        } else {
            "--secondary=no"
        };


        let cmd = if self.config.paired_end {
            format!(
                "minimap2 -ax {aligner_preset} {secondary_arg} -t {} {aligner_args} '{}' '{}' '{}' | samtools view -@ {} -hbF 12 - | samtools sort -@ {} - {output}",
                self.config.threads,
                self.config.index.display(),
                self.config.input[0].display(),
                self.config.input[1].display(),
                self.config.threads,
                self.config.threads,
            )
        } else {
            format!(
                "minimap2 -ax {aligner_preset} {secondary_arg} -t {} {aligner_args} '{}' '{}' | samtools view -@ {} -hbF 12 - | samtools sort -@ {} - {output}",
                self.config.threads,
                self.config.index.display(),
                self.config.input[0].display(),
                self.config.threads,
                self.config.threads,
            )
        };

        self.run_command(&cmd)?;


        let alignment = match &self.config.output {
            Some(output) => {
                Some(ReadAlignment::from_bam(
                    &output, &self.reference, &self.filter
                )?)
            },
            None => None
        };

        Ok(alignment)
    }
    fn run_bowtie2(&self) -> Result<Option<ReadAlignment>, VircovError> {
        let aligner_args = self.config.args.as_deref().unwrap_or("");

        let output = match &self.config.output {
            Some(output) => format!("> {}", output.display()),
            None => String::from("")
        };

        let index = match self.config.create_index {
            true => {
                let index_cmd = format!("bowtie2-build {} {}", self.config.index.display(), self.config.index.display());
                self.run_command_no_stdout(&index_cmd)?;
                &self.config.index
            },
            false => &self.config.index
        };

        let secondary_arg = if self.config.secondary {
            "" // TODO: documentation default behaviour of Bowtie2
        } else {
            "-k 0"
        };

        let cmd = if self.config.paired_end {
            format!(
                "bowtie2 {secondary_arg} -x '{}' -1 '{}' -2 '{}' --mm -p {} {aligner_args} | samtools view -@ {} -hbF 12 - | samtools sort -@ {} - {output} ",
                index.display(),
                self.config.input[0].display(),
                self.config.input[1].display(),
                self.config.threads,
                self.config.threads,
                self.config.threads,
            )
        } else {
            format!(
                "bowtie2 {secondary_arg} -x '{}' -U '{}' --mm -p {} {aligner_args} | samtools view -@ {} -hbF 12 - | samtools sort -@ {} - {output}",
                index.display(),
                self.config.input[0].display(),
                self.config.threads,
                self.config.threads,
                self.config.threads,
            )
        };

        self.run_command(&cmd)?;

        let alignment = match &self.config.output {
            Some(output) => {
                Some(ReadAlignment::from_bam(
                    &output, &self.reference, &self.filter
                )?)
            },
            None => None
        };

        Ok(alignment)

    }
    fn run_strobealign(&self) -> Result<Option<ReadAlignment>, VircovError> {
        let aligner_args = self.config.args.as_deref().unwrap_or("");

        let output = match &self.config.output {
            Some(output) => format!("> {}", output.display()),
            None => String::from("")
        };

        let secondary_arg = if self.config.secondary {
            "-N 50"  // TODO: documentation default behaviour of Strobealign
        } else {
            "-N 0"
        };

        let cmd = if self.config.paired_end {
            format!(
                "strobealign {secondary_arg} -t {} {aligner_args} '{}' '{}' '{}' | samtools view -@ {} -hbF 12 - | samtools sort -@ {} - {output}",
                self.config.threads,
                self.config.index.display(),
                self.config.input[0].display(),
                self.config.input[1].display(),
                self.config.threads,
                self.config.threads,
            )
        } else {
            format!(
                "strobealign {secondary_arg} -t {} {aligner_args} '{}' '{}' | samtools view -@ {} -hbF 12 - | samtools sort -@ {} - {output}",
                self.config.threads,
                self.config.index.display(),
                self.config.input[0].display(),
                self.config.threads,
                self.config.threads,
            )
        };
        self.run_command(&cmd)?;


        let alignment = match &self.config.output {
            Some(output) => {
                Some(ReadAlignment::from_bam(
                    &output, &self.reference, &self.filter
                )?)
            },
            None => None
        };

        Ok(alignment)
    }
    fn run_command(&self, cmd: &str) -> Result<(), VircovError> {
        log::debug!("Running command: {}", cmd);

        let status = Command::new("sh")
            .arg("-c")
            .arg(cmd)
            .stderr(Stdio::null())
            .status()
            .map_err(|e| VircovError::CommandExecutionFailed(cmd.to_string(), e.to_string()))?;

        if !status.success() {
            return Err(VircovError::CommandFailed(cmd.to_string(), status.code().unwrap_or(-1)));
        }

        Ok(())
    }
    fn run_command_no_stdout(&self, cmd: &str) -> Result<(), VircovError> {
        log::debug!("Running command: {}", cmd);

        let status = Command::new("sh")
            .arg("-c")
            .arg(cmd)
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .map_err(|e| VircovError::CommandExecutionFailed(cmd.to_string(), e.to_string()))?;

        if !status.success() {
            return Err(VircovError::CommandFailed(cmd.to_string(), status.code().unwrap_or(-1)));
        }

        Ok(())
    }

}

/// Paf struct
#[derive(Debug, Clone)]
pub struct ReadAlignment {
    pub target_intervals: TargetIntervals,
    pub target_sequences: Option<HashMap<String, FastaRecord>>,
    pub target_exclude: Option<Vec<String>>,
    pub reference: ReferenceConfig,
    pub filter: FilterConfig,
}
impl ReadAlignment {
    pub fn new(
        target_intervals: TargetIntervals,
        reference: &ReferenceConfig,
        filter: &FilterConfig,
    ) -> Result<Self, VircovError> {

        let target_sequences = parse_reference_fasta(
            reference.fasta.clone()
        )?;

        let target_exclude = parse_exclude_file(
            reference.exclude.clone()
        )?;

        Ok(Self {
            target_intervals,
            target_sequences,
            target_exclude,
            reference: reference.clone(),
            filter: filter.clone(),
        })
    }
    pub fn from(
        alignment: &PathBuf,
        reference: &ReferenceConfig,
        filter: &FilterConfig,
        alignment_format: Option<AlignmentFormat>,
    ) -> Result<Self, VircovError> {
        match alignment_format {
            Some(format) => match format {
                AlignmentFormat::Paf | AlignmentFormat::Gaf => ReadAlignment::from_paf(
                    alignment, reference, filter
                ),
                
                AlignmentFormat::Sam | AlignmentFormat::Bam | 
                AlignmentFormat::Cram  => ReadAlignment::from_bam(
                    alignment, reference, filter
                ),
                #[cfg(not(feature = "htslib"))]
                _ =>  Err(VircovError::AlignmentInputFormatInvalid),
            },
            None => match alignment.extension().map(|s| s.to_str()) {
                Some(Some("paf")) | 
                Some(Some("paf.gz")) | Some(Some("paf.xz")) | 
                Some(Some("paf.bz")) | Some(Some("paf.bz2")) => ReadAlignment::from_paf(
                    alignment, reference, filter
                ),
                Some(Some("gaf")) | 
                Some(Some("gaf.gz")) | Some(Some("gaf.xz")) | 
                Some(Some("gaf.bz")) | Some(Some("gaf.bz2")) => ReadAlignment::from_paf(
                    alignment, reference, filter
                ),
                
                Some(Some("sam")) | Some(Some("bam")) | Some(Some("cram")) => ReadAlignment::from_bam(
                    alignment, reference, filter
                ),
                _ => Err(VircovError::AlignmentInputFormatNotRecognized),
            }
        }
    }
    pub fn from_paf(
        path: &PathBuf,
        reference: &ReferenceConfig,
        filter: &FilterConfig
    ) -> Result<Self, VircovError> {

        let reader: Box<dyn BufRead> = if path.to_str() == Some("-") {
            Box::new(BufReader::new(std::io::stdin()))
        } else {
            let (reader, _) = niffler::from_path(path)?;
            Box::new(BufReader::new(reader))
        };
        
        let mut target_intervals: BTreeMap<String, Vec<Interval<usize, String>>> = BTreeMap::new();

        for result in reader.lines() {
            let record: PafRecord = PafRecord::from_str(result?)?;
            if record.query_aligned_length() >= filter.min_query_length
                && record.query_coverage() >= filter.min_query_coverage
                && record.mapq >= filter.min_mapq
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

        Self::new(target_lappers, reference, filter)
    }
    pub fn from_paf_stdout(
        stdout: ChildStdout,
        reference: &ReferenceConfig,
        filter: &FilterConfig
    ) -> Result<Self, VircovError> {

        let reader = BufReader::new(stdout);
        
        let mut target_intervals: BTreeMap<String, Vec<Interval<usize, String>>> = BTreeMap::new();

        for result in reader.lines() {
            let record: PafRecord = PafRecord::from_str(result?)?;
            if record.query_aligned_length() >= filter.min_query_length
                && record.query_coverage() >= filter.min_query_coverage
                && record.mapq >= filter.min_mapq
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

        Self::new(target_lappers, reference, filter)
    }
    
    pub fn from_bam(
        path: &PathBuf,
        reference: &ReferenceConfig,
        filter: &FilterConfig
    ) -> Result<Self, VircovError> {

        let mut reader = if path.to_str() == Some("-") {
            bam::Reader::from_stdin()?
        } else {
            bam::Reader::from_path(path)?
        };
        

        let header_view = reader.header().to_owned();

        let mut target_intervals: BTreeMap<String, Vec<Interval<usize, String>>> = BTreeMap::new();

        for result in reader.records() {
            let record = result?;
            if record.is_unmapped() {
                continue;
            }

            let bam_record = BamRecord::from(&record, &header_view)?;

            if bam_record.qalen >= filter.min_query_length
                && bam_record.query_coverage() >= filter.min_query_coverage
                && bam_record.mapq >= filter.min_mapq
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

        Self::new(target_lappers, reference, filter)
    }

    /// Compute coverage distribution by target sequence
    pub fn coverage(
        &self,
        annotation_config: &AnnotationConfig,
        tags: bool,
        zero: bool,
    ) -> Result<Vec<Coverage>, VircovError> {

        let mut coverage_fields: Vec<Coverage> = Vec::new();
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

            let target_tags = match tags {
                true => {
                    merged_targets
                        .depth()
                        .collect::<Vec<Interval<usize, usize>>>()
                }
                false => Vec::new(),
            };

            let target_description = match target_record {
                None => "-".to_string(),
                Some(fasta_record) => match fasta_record.description() {
                    None => "-".to_string(),
                    Some(descr) => descr.to_owned(),
                },
            };
            
            // exclude terms are parsed as lowercase (case-insensitive matching)
            let target_description_decap = target_description.to_lowercase(); 

            // blacklist application
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

            let region_filter_passed: bool = match self.filter.min_scan_regions_coverage {
                Some(reg_cov_threshold) => {
                    // Apply the region filter only if the reference sequence
                    // coverage is below a coverage threshold - only apply this if
                    // the threshold is bigger than zero so that pipelines can still
                    // apply the option with a 0 value
                    if (reg_cov_threshold > 0.) && (target_cov < reg_cov_threshold) {
                        target_cov_n >= self.filter.min_scan_regions
                    } else {
                        // Otherwise, if coverage is above the threshold, do
                        // not apply the region filter
                        true
                    }
                }
                _ => target_cov_n >= self.filter.min_scan_regions,
            };

            let annotation = Annotation::from(&target_description, annotation_config);
            
            let reads_aligned = targets.len() as u64;

            if target_len >= self.filter.min_reference_length
                && reads_aligned >= self.filter.min_scan_alignments
                && region_filter_passed
                && target_cov >= self.filter.min_scan_coverage
                && unique_reads_n >= self.filter.min_scan_reads
                && !exclude_target_sequence
            {
                coverage_fields.push(Coverage {
                    reference: target_name.to_owned(),
                    regions: target_cov_n,
                    reads: unique_reads_n,
                    alignments: reads_aligned,
                    bases: target_cov_bp,
                    length: target_len,
                    coverage: target_cov,
                    description: target_description,
                    group: annotation.group,
                    segment: annotation.segment,
                    name: annotation.name,
                    read_id: unique_read_ids,
                    intervals: target_tags,
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

                            let annotation = Annotation::from(&descr, annotation_config);

                            coverage_fields.push(Coverage {
                                reference: ref_seq.to_owned(),
                                regions: 0,
                                reads: 0,
                                alignments: 0,
                                bases: 0,
                                length: record.sequence().len() as u64,
                                coverage: 0.0,
                                description: descr,
                                group: annotation.group,
                                segment: annotation.segment,
                                name: annotation.name,
                                read_id: Vec::new(),
                                intervals: Vec::new(),
                            });
                        }
                    }
                }
                None => return Err(VircovError::ZeroReferenceSequences),
            }
        }

        Ok(coverage_fields)
    }
    pub fn write_reads(coverage: &[Coverage], output: PathBuf, ref_id: bool) -> Result<(), VircovError> {

        let mut file_handle = File::create(&output)?;

        if ref_id {
            let mut read_ids = HashSet::new();
            for cov in coverage {
                for read_id in &cov.read_id {
                    read_ids.insert((read_id, &cov.reference));
                }
            }
            for (read_id, ref_id) in read_ids {
                writeln!(file_handle, "{}\t{}", &read_id, ref_id)?;
            }
        } else {
            let mut read_ids = HashSet::new();
            for cov in coverage {
                for read_id in &cov.read_id {
                    read_ids.insert(read_id);
                }
            }
            for read_id in read_ids {
                writeln!(file_handle, "{}", &read_id)?;
            }
        }
        

        Ok(())
    }
    pub fn write_grouped_reads(&self, coverage_fields: &[GroupedCoverage], output: PathBuf) -> Result<(), VircovError> {

        std::fs::create_dir_all(&output)?;
        for field in coverage_fields.iter() {
            let sanitized_name = field.group.replace(' ', "__");

            let file_path = output.join(sanitized_name).with_extension("txt");
            let mut file_handle = File::create(file_path.as_path())?;
            for read_id in field.read_id.iter() {
                writeln!(file_handle, "{}", &read_id)?
            }
        }
        Ok(())
    }

    pub fn write_tsv(coverage: &[Coverage], output: &PathBuf, header: bool) -> Result<(), VircovError> {

        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(header)
            .from_path(&output)
            .unwrap();
        
        for rec in coverage {
            writer.serialize(rec)?;
        }
        Ok(())

    }

    pub fn print_coverage_table(coverage: &mut [Coverage], sort: bool) {

        if sort {
            coverage.sort_by(|a, b| {
                b.coverage.partial_cmp(&a.coverage)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
        }
        let table_fields = coverage
            .iter()
            .map(CoverageTableFields::from)
            .collect::<Vec<CoverageTableFields>>();
        
        let mut table = Table::new(table_fields);

        table.modify(
            Columns::new(7..), 
            Width::wrap(32).keep_words()
        ).with(
            Style::modern()
        );

        eprintln!("{}", table);
    }
    /// Print the coverage distributions to console as a text-based coverage plot
    pub fn coverage_plots(
        &self,
        coverage: &[Coverage],
        max_width: u64,
    ) -> Result<(), VircovError> {
        println!();
        for (target_name, targets) in &self.target_intervals {
            let target_record = match &self.target_sequences {
                None => None,
                Some(seqs) => seqs.get(target_name),
            };
            let target_len = match target_record {
                None => return Err(VircovError::CovPlotSeqLengthError()),
                Some(fasta_record) => fasta_record.sequence().len() as u64,
            };

            if coverage
                .iter()
                .map(|x| &x.reference)
                .any(|x| x == target_name)
            {
                let covplot = CovPlot::new(targets, target_len, max_width)?;
                covplot.to_console(target_name, target_len, Color::Red)?;
            }
        }

        Ok(())
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
    pub fn from(record: &bam::Record, header: &HeaderView) -> Result<Self, VircovError> {
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
    pub fn from_str(paf: String) -> Result<Self, VircovError> {
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
