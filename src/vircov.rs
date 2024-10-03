use crate::alignment::{parse_reference_fasta, Aligner, AlignmentFormat, Coverage, GroupedCoverage, Preset, ReadAlignment, SelectHighest, VircovAligner};
use crate::consensus::{ConsensusAssembler, ConsensusRecord, VircovConsensus};
use crate::error::VircovError;
use crate::terminal::{CoverageArgs, RunArgs};
use crate::utils::{get_file_component, FileComponent};

use itertools::Itertools;
use anyhow::Result;
use indexmap::IndexMap;
use noodles::fasta::Record;
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};
use tabled::settings::object::Columns;
use tabled::settings::{Style, Width};
use tabled::{Table, Tabled};
use std::collections::btree_map::Entry;
use std::collections::{BTreeMap, HashMap};
use std::fs::{create_dir_all, remove_dir_all, remove_file};
use std::path::PathBuf;
use std::io::{BufRead, BufReader};
use rayon::iter::{ParallelIterator, IntoParallelIterator};

/// Vircov application structure.
pub struct Vircov {
    config: VircovConfig,
    references: Option<HashMap<String, Record>>
}
impl Vircov {
    pub fn from(config: VircovConfig) -> Result<Self, VircovError> {

        Ok(Self { 
            references: parse_reference_fasta(
                config.reference.fasta.clone()
            )?,
            config, 
        })
    }
    pub fn run(
        &self, 
        tsv: &PathBuf,
        parallel: usize, 
        threads: usize, 
        consensus: bool, 
        keep: bool,
        table: bool,
        select_by: SelectHighest,
        alignment: Option<PathBuf>,
        alignment_format: Option<AlignmentFormat>
    ) -> Result<(), VircovError> {

        log::info!("Alignment scan against reference database ({})", self.config.alignment.aligner);
        let read_alignment = match alignment {
            Some(path) => self.alignment(&path, alignment_format)?,
            None => self.align()?
        };

        log::info!("Coverage metrics computation from alignment scan");
        let coverage = read_alignment.coverage(false, false)?;

        log::info!(
            "Reference sequence group determination using '{}' header field", 
            self.config.reference.annotation.group.clone().unwrap_or("None".to_string())
        );
        let grouped = self.group(&coverage)?;

        log::info!("Reference genome selection by highest '{}'", select_by);
        let selections = self.select(
            &grouped, 
            SelectHighest::Coverage,
            None
        )?;

        let group_read_files = if self.config.alignment.remap_group_reads {
            log::info!("Reference grouping reads will be used for remapping");
            Some(self.extract_group_reads(
                &grouped,
                &self.config.outdir.clone()
            )?)
        } else {
            log::info!("All input reads will be used for remapping");
            None
        };
        

        log::info!("Starting remapping and consensus assembly for each genome ({})", self.config.alignment.aligner);
        let (consensus, remap) = self.remap(
            &selections, 
            group_read_files,
            &self.config.outdir.clone(), 
            parallel, 
            threads, 
            consensus, 
            keep
        )?;

        let scan = selections.values().flatten().cloned().collect();
        let summary = self.summary(consensus, scan, remap)?;

        log::info!("Writing output table to: {}", tsv.display());
        if table { summary.print_table(true); }
        summary.write_tsv(&tsv, true)?;

        if !keep {
            remove_dir_all(&self.config.outdir)?;
        }

        Ok(())
    }
    pub fn run_coverage(
        &self,
        tsv: &PathBuf,
        tags: bool,
        zero: bool,
        table: bool,
        read_id: Option<PathBuf>,
        alignment: Option<PathBuf>,
        alignment_format: Option<AlignmentFormat>
    ) -> Result<(), VircovError> {

        let read_alignment = match alignment {
            Some(path) => self.alignment(&path, alignment_format)?,
            None => self.align()?
        };

        let mut coverage = read_alignment.coverage(tags, zero)?;

        if table {
            ReadAlignment::print_coverage_table(&mut coverage, true);
        }

        ReadAlignment::write_tsv(&coverage, tsv, true)?;

        if let Some(path) = read_id {
            ReadAlignment::write_reads(&coverage, path, false)?;
        }

        Ok(())
    }
    pub fn alignment(
        &self, 
        path: &PathBuf, 
        alignment_format: Option<AlignmentFormat>
    ) -> Result<ReadAlignment, VircovError> {
        ReadAlignment::from(
            path,
            &self.config.reference, 
            &self.config.filter, 
            alignment_format
        )
    }
    pub fn align(
        &self,
    ) -> Result<ReadAlignment, VircovError> {

        let aligner = VircovAligner::from(
            &self.config.alignment, 
            &&self.config.reference, 
            &self.config.filter
        );

        aligner.check_aligner_dependency(&aligner.config.aligner)?;

        let read_alignment = aligner
            .run_aligner()?
            .ok_or(VircovError::AlignerStdoutMisconfigured)?;

        Ok(read_alignment)
    }
    pub fn group(
        &self,
        coverage: &[Coverage]
    ) -> Result<Vec<GroupedCoverage>, VircovError> {
        let mut coverage_groups: BTreeMap<String, Vec<&Coverage>> = BTreeMap::new();

        for cov in coverage {
            if let Some(group) = &cov.group {
                match coverage_groups.entry(group.to_string()) {
                    Entry::Occupied(mut entry) => {
                        entry.get_mut().push(cov);
                    }
                    Entry::Vacant(entry) => {
                        entry.insert(vec![cov]);
                    }
                }
            } else {
                return Err(VircovError::GroupAnnotationMissing(cov.reference.to_owned()))
            }
        }

        let mut grouped_coverage_all: Vec<GroupedCoverage> = Vec::new();
        for (group, coverage) in coverage_groups {
            
            let group = group.trim().replace(
                self.config.reference.annotation.group
                    .as_ref()
                    .unwrap(),
                 ""
            ); // TODO

            let grouped_coverage_data = GroupedCoverage::from_coverage(&group, coverage)?;

            if grouped_coverage_data.total_regions >= self.config.filter.min_grouped_regions
                && grouped_coverage_data.mean_coverage >= self.config.filter.min_grouped_mean_coverage
                && grouped_coverage_data.total_alignments >= self.config.filter.min_grouped_alignments
                && grouped_coverage_data.total_reads >= self.config.filter.min_grouped_reads
            {
                grouped_coverage_all.push(grouped_coverage_data);
            }
        }

        Ok(grouped_coverage_all)
    }
    pub fn extract_group_reads(
        &self,
        grouped_coverage: &Vec<GroupedCoverage>,
        outdir: &PathBuf
    ) -> Result<HashMap<String, Vec<PathBuf>>, VircovError> {  // group, read fastq

        let mut group_reads_fastq = HashMap::new();
        for group_coverage in grouped_coverage {
            
            let output = self.config.alignment.get_output_read_paths(
                &group_coverage.group, 
                "fastq.gz", 
                outdir
            )?;

            group_coverage.write_group_reads(
                self.config.alignment.input.clone(), 
                output.clone()
            )?;

            // Same keys as the group selection coverage structs for remapping
            group_reads_fastq.insert(group_coverage.group.clone(), output);
        }

        Ok(group_reads_fastq)
    }
    pub fn select(
        &self,
        grouped_coverage: &Vec<GroupedCoverage>, // If grouped, these are grouped fields
        select_by: SelectHighest,
        outdir: Option<PathBuf>,
    ) -> Result<HashMap<String, Vec<Coverage>> , VircovError> {
        
        let refs = self.references
            .as_ref()
            .ok_or(VircovError::GroupSequenceError)?;

        if let Some(ref path) = outdir {
            std::fs::create_dir_all(path)?;
        }

        let mut selected_reference_coverage: HashMap<String, Vec<Coverage>> = HashMap::new();

        for group_coverage in grouped_coverage {
            let mut grouped_segments: BTreeMap<String, Vec<Coverage>> = BTreeMap::new();

            for cov in &group_coverage.coverage {

                if let Some(segment) = &cov.segment {

                    if segment != self.config.reference.annotation.segment_na.as_ref().ok_or(VircovError::NoSegmentFieldProvided)? {
                        match grouped_segments.entry(segment.clone()) {
                            Entry::Occupied(mut entry) => {
                                entry.get_mut().push(cov.clone());
                            }
                            Entry::Vacant(entry) => {
                                entry.insert(vec![cov.clone()]);
                            }
                        }
                    }   
                }
            }

            let segmented = !grouped_segments.is_empty();

            // Filter reference sequences (contained in tags) by highest number of (unique) reads, alignments or coverage

            if segmented {

                let mut selected_segments: BTreeMap<String, Coverage> = BTreeMap::new();

                for (segment, cov_fields) in grouped_segments {
                    let selected = match select_by {
                        SelectHighest::Reads => cov_fields.iter().max_by_key(|x| x.reads),
                        SelectHighest::Coverage => cov_fields.iter().max_by_key(|x| OrderedFloat(x.coverage)),
                        SelectHighest::Alignments => cov_fields.iter().max_by_key(|x| x.alignments)
                    };

                    match selected {
                        Some(selected) => {
                            selected_segments
                                .entry(segment)
                                .or_insert_with(|| selected.clone());
                        }
                        None => return Err(VircovError::GroupSelectReference),
                    }
                }

                if let Some(path) = &outdir {

                    let file_handle = std::fs::File::create(
                        path.join(
                            format!("{}.fasta", group_coverage.group)
                        )
                    )?;
                    let mut writer = noodles::fasta::Writer::new(file_handle);

                    for (_, segment_cov) in selected_segments {
                        let seq = &refs[&segment_cov.reference];
                        writer.write_record(seq)?;
                        
                        selected_reference_coverage.entry(group_coverage.group.clone())
                            .and_modify(|x| x.push(segment_cov.clone()))
                            .or_insert(vec![segment_cov]);
                    }
                } else {
                    for (_, segment_cov) in selected_segments {
                        selected_reference_coverage.entry(group_coverage.group.clone())
                            .and_modify(|x| x.push(segment_cov.clone()))
                            .or_insert(vec![segment_cov]);
                    }
                }

            } else {

                // If not segmented, simply select the best reference sequence using the
                // selection metric and select the sequence from the header fields
                // to write to FASTA

                let max = match select_by {
                    SelectHighest::Reads => group_coverage.coverage.iter().max_by_key(|x| x.reads),
                    SelectHighest::Coverage => group_coverage.coverage.iter().max_by_key(|x| OrderedFloat(x.coverage)),
                    SelectHighest::Alignments => group_coverage.coverage.iter().max_by_key(|x| x.alignments)
                };

                match max {
                    Some(field) => {
                        let seq = &refs[&field.reference];

                        if let Some(path) = &outdir {
                            
                            let file_handle = std::fs::File::create(
                                path.join(
                                    format!("{}.fasta", group_coverage.group)
                                )
                            )?;
                            let mut writer = noodles::fasta::Writer::new(file_handle);
                            writer.write_record(seq)?;
                        }

                        selected_reference_coverage.entry(group_coverage.group.clone())
                            .and_modify(|x| x.push(field.clone()))
                            .or_insert(vec![field.clone()]);
                    }
                    _ => return Err(VircovError::GroupSelectReference),
                }
            }

        }

        Ok(selected_reference_coverage)
    }
    pub fn remap(
        &self, 
        selections: &HashMap<String, Vec<Coverage>>, 
        group_read_files: Option<HashMap<String, Vec<PathBuf>>>,
        outdir: &PathBuf,
        parallel: usize, 
        threads: usize,
        consensus: bool,
        keep: bool
    ) -> Result<(Vec<ConsensusRecord>, Vec<Coverage>), VircovError> {

        let refs = self.references
            .as_ref()
            .ok_or(VircovError::GroupSequenceError)?;

        rayon::ThreadPoolBuilder::new()
            .num_threads(parallel)
            .build()
            .expect("Failed to create thread pool")
            .install(|| -> Result<(Vec<ConsensusRecord>, Vec<Coverage>), VircovError> {

            let results: Vec<Result<(Vec<Vec<ConsensusRecord>>, Vec<Coverage>), VircovError>> = selections
                .into_par_iter()
                .map(|(group, coverage)| -> Result<(Vec<Vec<ConsensusRecord>>, Vec<Coverage>), VircovError> {
                    
                    let remap_id = group;

                    let group_read_files = match &group_read_files {
                        Some(group_read_map) => group_read_map.get(remap_id),
                        None => None
                    };

                    log::info!("Reference group '{}' launched remapping stage (n = {})", group, coverage.len());

                    let remap_reference = outdir.join(remap_id).with_extension("fasta");
                    
                    self.extract_remap_sequences(
                        &remap_reference, 
                        coverage, 
                        refs
                    )?;

                    let bam = outdir.join(remap_id).with_extension("bam");

                    let remap_aligner = VircovAligner::from(
                        &self.config.alignment.remap_config(
                            &remap_reference, 
                            None, 
                            group_read_files.cloned(),
                            Some(bam.clone()), 
                            threads
                        ),
                        &self.config.reference.remap_config(
                            &remap_reference
                        ),
                        &self.config.filter,
                    );

                    let remap_coverage = remap_aligner
                        .run_aligner()?
                        .unwrap()
                        .coverage( true, false)?;

                    let mut consensus_records = Vec::new();

                    if consensus {
                        for ref_cov in &remap_coverage {
                            if ref_cov.coverage >= self.config.filter.min_remap_coverage {
    
                                // Reference output must have extension '.fa' matching auto-extension from iVar
                                let (filter_reference, output_name) = match &ref_cov.segment {
                                    Some(segment) => {
                                        let seg = self.config.reference.annotation.segment_name_file(segment);
                                        (
                                            Some(ref_cov.reference.clone()), 
                                            format!("{group}.{seg}.consensus.fa")
                                        )
                                    }, 
                                    None => (
                                        None, 
                                        format!("{group}.nan.consensus.fa")
                                    )
                                };
    
                                log::info!("Creating consensus assembly with iVar: {output_name}");
                                let consensus = VircovConsensus::new(
                                    ConsensusConfig::with_default(
                                        &bam.clone(), 
                                        &outdir.join(output_name), 
                                        filter_reference,
                                        Some(format!("{} {}", ref_cov.reference, ref_cov.description)),
                                        self.config.consensus.min_quality,
                                        self.config.consensus.min_frequency,
                                        self.config.consensus.min_depth
                                    )
                                )?;
    
                                let records = consensus.assemble()?;
    
                                consensus_records.push(records)
                            }
                        }
                    }

                    if !keep {
                        remove_file(&bam)?;
                        remove_file(&remap_reference)?;

                        let bai = bam.with_extension("bai");
                        if bai.exists() { remove_file(bai)? };
                        
                    }

                    Ok((consensus_records, remap_coverage))
                })
                .collect();

            let mut remap_data = Vec::new();
            let mut consensus_data = Vec::new();

            for result in results {
                let (consensus, remap) = result?;
                
                for consensus_record in consensus.into_iter().flatten() {
                    consensus_data.push(consensus_record)
                }
                for cov in remap.into_iter() {
                    remap_data.push(cov)
                }
            }

            Ok((consensus_data, remap_data))
        })

    }
    fn extract_remap_sequences(
        &self, 
        remap_fasta: &PathBuf, 
        coverage: &Vec<Coverage>, 
        refs: &HashMap<String, Record>
    ) -> Result<(), VircovError> {


        let file_handle = std::fs::File::create(&remap_fasta)?;
        let mut writer = noodles::fasta::Writer::new(file_handle);

        let refseqs: Vec<&noodles::fasta::Record> = coverage
            .iter()
            .map(|ref_cov| {
                refs.get(&ref_cov.reference)
                    .ok_or(VircovError::AlignmentInputFormatInvalid)
            })
            .collect::<Result<_, _>>()?;

        for seq in refseqs {
            writer.write_record(seq)?;
        }

        Ok(())
    }
    pub fn summary(
        &self, 
        consensus: Vec<ConsensusRecord>, 
        scan: Vec<Coverage>, 
        remap: Vec<Coverage>
    ) -> Result<VircovSummary, VircovError> {

        VircovSummary::new(
            None, 
            None, 
            None,
            scan
                .iter()
                .map(|cov| AlignmentRecord::from(
                    cov.clone(), 
                    false, 
                    None
                ))
                .collect(),
            remap.iter().map(|cov| AlignmentRecord::from(
                cov.clone(), 
                false, 
                None
            )).collect(),
            consensus,
            &&self.config.reference.annotation
        )
    }
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VircovConfig {
    pub outdir: PathBuf,
    pub alignment: AlignmentConfig,
    pub reference: ReferenceConfig,
    pub filter: FilterConfig,
    pub coverage: CoverageConfig,
    pub consensus: ConsensusConfig,
    pub subtype: SubtypeConfig
} 
impl VircovConfig {
    pub fn from_run_args(args: &RunArgs) -> Result<Self, VircovError> {

        let outdir = match &args.workdir {
            None => PathBuf::from("."),
            Some(dir) => dir.to_owned()
        };

        if !outdir.exists() && !outdir.is_dir() {
            create_dir_all(&outdir)?
        }

        Ok(Self {
            outdir: outdir.clone(),
            alignment: AlignmentConfig::with_default(
                &args.input, 
                &args.index, 
                args.aligner.clone(), 
                args.preset.clone(), 
                false, 
                args.secondary,
                Some(outdir.join("scan.bam")),
                args.scan_threads,
                !args.remap_all
            ),
            reference: ReferenceConfig::with_default(
                Some(args.reference.clone())
            ),
            consensus: ConsensusConfig::with_default_from_args(
                args.min_consensus_depth,
                args.min_consensus_frequency,
                args.min_consensus_depth
            ),
            filter: FilterConfig::with_default(args.min_remap_coverage),
            coverage: CoverageConfig {},
            subtype: SubtypeConfig {}
        })
    }
    pub fn from_coverage_args(args: &CoverageArgs) -> Result<Self, VircovError> {

        let outdir = match &args.workdir {
            None => PathBuf::from("."),
            Some(dir) => dir.to_owned()
        };

        if !outdir.exists() && !outdir.is_dir() {
            create_dir_all(&outdir)?
        }

        Ok(Self {
            outdir: outdir.clone(),
            alignment: AlignmentConfig::with_default(
                &args.input, 
                &args.index, 
                args.aligner.clone(), 
                args.preset.clone(), 
                false, 
                args.secondary,
                Some(outdir.join("scan.bam")),
                args.threads,
                false
            ),
            reference: ReferenceConfig::with_default(
                args.reference.clone()
            ),
            filter: FilterConfig::default(),
            consensus: ConsensusConfig::default(), // not used
            coverage: CoverageConfig {},
            subtype: SubtypeConfig {}
        })
    }
}
impl Default for VircovConfig {
    fn default() -> Self {
        Self {
            outdir: PathBuf::from("."),
            alignment: AlignmentConfig::default(),
            reference: ReferenceConfig::default(),
            filter: FilterConfig::default(),
            consensus: ConsensusConfig::default(),
            coverage: CoverageConfig {},
            subtype: SubtypeConfig {}
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignmentConfig {
    pub aligner: Aligner,
    pub preset: Option<Preset>,
    pub args: Option<String>,
    pub input: Vec<PathBuf>,
    pub paired_end: bool,
    pub index: PathBuf,
    pub create_index: bool,
    pub secondary: bool,
    pub output: Option<PathBuf>,
    pub threads: usize,
    pub remap_group_reads: bool
}
impl AlignmentConfig {
    pub fn remap_config(
        &self, 
        index: &PathBuf, 
        args: Option<String>,
        input: Option<Vec<PathBuf>>,
        output: Option<PathBuf>, 
        threads: usize,
    ) -> Self {
        
        let mut config = self.clone();

        // Input read files are provided explicitly
        // for example from group-extracted read files

        // Otherwise user configured input files 
        // for the initial alignment (all reads)
        // are used

        if let Some(input) = input {
            config.input = input;
        };


        config.create_index = true;
        config.args = args;
        config.index = index.clone();
        config.threads = threads;
        config.output = output;

        config
    }
    pub fn with_default(
        input: &Vec<PathBuf>, 
        index: &PathBuf, 
        aligner: Aligner, 
        preset: Option<Preset>, 
        create_index: bool, 
        secondary: bool, 
        output: Option<PathBuf>, 
        threads: usize,
        remap_group_reads: bool
    ) -> Self {
        Self {
            input: input.clone(),
            index: index.clone(),
            aligner,
            preset,
            output,
            create_index,
            secondary,
            threads,
            remap_group_reads,
            ..Default::default()
        }
    }
    pub fn get_output_read_paths(&self, id: &str, ext: &str, outdir: &PathBuf) -> Result<Vec<PathBuf>, VircovError> {

        // Paired-end short reads
        if self.paired_end {
            return Ok(vec![
                outdir.join(format!("{id}_R1.{ext}")),
                outdir.join(format!("{id}_R2.{ext}"))
            ])
        } else {
            return Ok(vec![
                outdir.join(format!("{id}.{ext}")),
            ])
        }

    }
}
impl Default for AlignmentConfig {
    fn default() -> Self {
        Self {
            aligner: Aligner::Minimap2,
            preset: Some(Preset::Sr),
            args: None,
            input: Vec::from([PathBuf::from("smoke_R1.fastq.gz"), PathBuf::from("smoke_R2.fastq.gz")]),
            paired_end: true,
            create_index: false,
            secondary: false,
            index: PathBuf::from("reference.fasta"),
            output: Some(PathBuf::from("test.sam")),
            threads: 8,
            remap_group_reads: true
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReferenceConfig {
    pub fasta: Option<PathBuf>,
    pub exclude: Option<PathBuf>,
    pub annotation: AnnotationConfig,
}
impl ReferenceConfig {
    pub fn remap_config(
        &self, 
        fasta: &PathBuf, 
    ) -> Self {
        let mut config = self.clone();
        config.fasta = Some(fasta.to_owned());
        config
    }
    pub fn with_default(fasta: Option<PathBuf>) -> Self {
        Self {
            fasta: fasta.clone(),
            ..Default::default()
        }
    }
}
impl Default for ReferenceConfig {
    fn default() -> Self {
        Self {
            fasta: None,
            exclude: None,
            annotation: AnnotationConfig::virosaurus()
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FilterConfig {
    pub min_query_length: u64,
    pub min_query_coverage: f64,
    pub min_mapq: u8,
    pub min_scan_alignments: u64,
    pub min_scan_regions: u64,
    pub min_scan_coverage: f64,
    pub min_scan_reads: u64,
    pub min_scan_reference_length: u64,
    pub min_scan_regions_coverage: Option<f64>,
    pub min_grouped_regions: u64,
    pub min_grouped_mean_coverage: f64,
    pub min_grouped_alignments: u64,
    pub min_grouped_reads: u64,
    pub min_remap_coverage: f64
}
impl Default for FilterConfig {
    fn default() -> Self {
        Self {
            min_query_length: 0,
            min_query_coverage: 0.0,
            min_mapq: 0,
            min_scan_alignments: 0,
            min_scan_regions: 0,
            min_scan_coverage: 0.,
            min_scan_reads: 0,
            min_scan_reference_length: 0,
            min_scan_regions_coverage: None,
            min_grouped_regions: 0,
            min_grouped_mean_coverage: 0.0,
            min_grouped_alignments: 0,
            min_grouped_reads: 0,
            min_remap_coverage: 0.0
        }
    }
}
impl FilterConfig {
    pub fn with_default(
        min_remap_coverage: f64
    ) -> Self {
        Self {
            min_query_length: 0,
            min_query_coverage: 0.0,
            min_mapq: 0,
            min_scan_alignments: 0,
            min_scan_regions: 0,
            min_scan_coverage: 0.,
            min_scan_reads: 0,
            min_scan_reference_length: 0,
            min_scan_regions_coverage: None,
            min_grouped_regions: 0,
            min_grouped_mean_coverage: 0.0,
            min_grouped_alignments: 0,
            min_grouped_reads: 0,
            min_remap_coverage
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoverageConfig {
}


#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct ConsensusConfig {
    pub alignment: PathBuf,
    pub output: PathBuf,
    pub assembler: ConsensusAssembler,
    pub reference: Option<String>,
    pub header: Option<String>,
    pub args: Option<String>,
    pub mpileup: Option<String>,
    pub min_quality: usize,
    pub min_frequency: f64,
    pub min_depth: usize,
    pub missing: String,
}
impl ConsensusConfig {
    pub fn with_default(alignment: &PathBuf, output: &PathBuf, reference: Option<String>, header: Option<String>, min_quality: usize, min_frequency: f64, min_depth: usize) -> Self {
        Self {
            alignment: alignment.clone(),
            output: output.clone(),
            header: header.clone(),
            reference,
            min_quality,
            min_frequency,
            min_depth,
            ..Default::default()
        }
    }
    pub fn with_default_from_args(min_quality: usize, min_frequency: f64, min_depth: usize) -> Self {
        Self {
            min_quality,
            min_frequency,
            min_depth,
            ..Default::default()
        }
    }
}
impl Default for ConsensusConfig {
    fn default() -> Self {
        Self {
            alignment: PathBuf::from("remap.sam"),
            output: PathBuf::from("remap.consensus.fa"), // must be .fa for Ivar
            assembler: ConsensusAssembler::Ivar,
            reference: None,
            header: None,
            args: None,
            mpileup: None,
            min_quality: 20,
            min_frequency: 0.75,
            min_depth: 10, 
            missing: String::from("N"),

        }
    }
}


pub struct Annotation {
    pub group: Option<String>,
    pub name: Option<String>,
    pub segment: Option<String>
}
impl Annotation {
    pub fn from(description: &str, options: &AnnotationConfig) -> Self {

        let sep: String = options.sep.to_string();

        let fields: Vec<&str> = description.split(&sep).map(|field| field.trim()).collect();

        let taxid = match &options.group {
            Some(field_taxid) => {
                let taxid_fields: Vec<&str> = fields.iter().filter(|field| field.starts_with(field_taxid)).map(|x| *x).collect();
                match taxid_fields.first() {
                    Some(f) => Some(f.trim().trim_start_matches(field_taxid).to_string()),
                    _ => None
                }
            },
            None => None
        };

        let taxon = match &options.name {
            Some(field_taxon) => {
                let taxon_fields: Vec<&str> = fields.iter().filter(|field| field.starts_with(field_taxon)).map(|x| *x).collect();
                match taxon_fields.first() {
                    Some(f) => Some(f.trim().trim_start_matches(field_taxon).to_string()),
                    _ => None
                }
            },
            None => None
        };

        let segment = match &options.segment {
            Some(field_segment) => {
                let segment_fields: Vec<&str> = fields.iter().filter(|field| field.starts_with(field_segment)).map(|x| *x).collect();
                match segment_fields.first() {
                    Some(segment) => Some(segment.trim_start_matches(field_segment).to_string()),
                    None => None
                }
            },
            None => None
        };

        Self {
            group: taxid, name: taxon, segment
        }       

    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnnotationConfig {
    pub sep: String, 
    pub group: Option<String>, 
    pub name: Option<String>, 
    pub segment: Option<String>, 
    pub segment_na: Option<String>
}
impl AnnotationConfig {
    pub fn segment_name_file(&self, segment: &str) -> String {
        match &self.segment_na {
            Some(segment_nan) => {
                if segment == segment_nan {
                    "nan".to_string()
                } else {
                    segment.to_string()
                }
            },
            None => "nan".to_string()
        }
    }
    pub fn virosaurus() -> Self {
        Self {
            sep: String::from(";"),
            group: Some(String::from("taxid=")),
            name: Some(String::from("usual name=")),
            segment: Some(String::from("segment=")),
            segment_na: Some(String::from("N/A"))
        }
    }
    pub fn with_default(
        &self, 
        sep: String, 
        group: Option<String>, 
        segment: Option<String>, 
        segment_na: Option<String>
    ) -> Self {
        Self {
            sep,
            group,
            segment,
            segment_na,
            ..Default::default()
        }
    }
}
impl Default for AnnotationConfig {
    fn default() -> Self {
        Self {
            sep: ";".to_string(),
            group: None,
            name: None,
            segment: None,
            segment_na: None
        }
    }
}

pub fn display_option_string(opt: &Option<String>) -> String 
{
    match opt {
        Some(s) => format!("{}", s),
        None => format!(""),
    }
}

pub fn display_option_u64(opt: &Option<u64>) -> String 
{
    match opt {
        Some(s) => format!("{}", s),
        None => format!(""),
    }
}


pub fn display_option_f64(opt: &Option<f64>) -> String 
{
    match opt {
        Some(s) => format!("{:.2}", s),
        None => format!(""),
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Tabled)]
pub struct VircovRecord {
    #[tabled(skip)]
    #[serde(skip_serializing_if = "Option::is_none")]
    pub id: Option<String>,
    #[tabled(skip)]
    #[serde(skip_serializing_if = "Option::is_none")]
    pub index: Option<String>,
    #[tabled(skip)]
    #[serde(skip_serializing_if = "Option::is_none")]
    pub aligner: Option<Aligner>,
    #[tabled(skip)]
    #[serde(skip_serializing_if = "Option::is_none")]
    pub assembler: Option<ConsensusAssembler>,
    #[tabled(rename="Taxid")]
    #[tabled(display_with = "display_option_string")]
    pub taxid: Option<String>,
    #[tabled(rename="Name")]
    #[tabled(display_with = "display_option_string")]
    pub name: Option<String>,
    #[tabled(rename="Segment")]
    #[tabled(display_with = "display_option_string")]
    pub segment: Option<String>,
    #[tabled(rename="Reference")]
    pub reference: String,
    #[tabled(rename="Length")]
    pub reference_length: u64,
    #[tabled(skip)]
    #[tabled(rename="Scan Regions")]
    pub scan_regions: u64,
    #[tabled(skip)]
    #[tabled(rename="Scan Reads")]
    pub scan_reads: u64,
    #[tabled(skip)]
    #[tabled(rename="Scan Alignments")]
    pub scan_alignments: u64,
    #[tabled(skip)]
    #[tabled(rename="Scan Bases")]
    pub scan_bases_covered: u64,
    #[tabled(skip)]
    #[tabled(rename="Scan Coverage")]
    pub scan_coverage: f64,
    #[tabled(rename="Regions")]
    #[tabled(display_with = "display_option_u64")]
    pub remap_regions: Option<u64>,
    #[tabled(rename="Reads")]
    #[tabled(display_with = "display_option_u64")]
    pub remap_reads: Option<u64>,
    #[tabled(rename="Alignments")]
    #[tabled(display_with = "display_option_u64")]
    pub remap_alignments: Option<u64>,
    #[tabled(rename="Bases")]
    #[tabled(display_with = "display_option_u64")]
    pub remap_bases_covered: Option<u64>,
    #[tabled(rename="Coverage")]
    #[tabled(display_with = "display_option_f64")]
    pub remap_coverage: Option<f64>,
    // pub remap_mean_depth: Option<f64>,
    #[tabled(skip)]
    #[tabled(rename="Consensus Length")]
    #[tabled(display_with = "display_option_u64")]
    pub consensus_length: Option<u64>,
    #[tabled(skip)]
    #[tabled(rename="Consensus Missing")]
    #[tabled(display_with = "display_option_u64")]
    pub consensus_missing: Option<u64>,
    #[tabled(rename="Completeness")]
    #[tabled(display_with = "display_option_f64")]
    pub consensus_completeness: Option<f64>,
    #[tabled(skip)]
    pub reference_description: String,
}
impl VircovRecord {
    pub fn from(
        index: Option<PathBuf>, 
        aligner: Option<Aligner>, 
        assembler: Option<ConsensusAssembler>,
        annotation: Annotation, 
        scan_record: AlignmentRecord, 
        remap_record: Option<AlignmentRecord>, 
        consensus_record: Option<ConsensusRecord>, 
    ) -> Result<Self, VircovError> {
        
        let (
            remap_regions,
            remap_reads,
            remap_alignments,
            remap_bases_covered,
            remap_coverage
        ) = match remap_record { 
            Some(r) =>(
                Some(r.regions),
                Some(r.reads), 
                Some(r.alignments),
                Some(r.bases),
                Some(r.coverage*100.)
            ),
            None => (None, None, None, None, None)
        };

        let (
            consensus_length, consensus_missing, consensus_completeness
        ) = match consensus_record {
            Some(r) => (Some(r.length), Some(r.missing), Some(r.completeness)),
            None => (None, None, None)
        };


        Ok(Self {
            id: None,
            index: match index { 
                Some(index) => Some(get_file_component(
                    &index, 
                    FileComponent::FileName
                )?), 
                None => None 
            },
            aligner: aligner.clone(),
            assembler,
            reference: scan_record.reference,
            reference_length: scan_record.length,
            scan_regions: scan_record.regions,
            scan_reads: scan_record.reads,
            scan_alignments: scan_record.alignments,
            scan_bases_covered: scan_record.bases,
            scan_coverage: scan_record.coverage*100.,
            remap_regions,
            remap_reads,
            remap_alignments,
            remap_bases_covered,
            remap_coverage,
            consensus_length,
            consensus_missing,
            consensus_completeness,
            taxid: annotation.group,
            name: annotation.name,
            segment: annotation.segment,
            reference_description: scan_record.description,
        })
    }
    pub fn from_coverage(
        coverage: Coverage,
        intervals: bool,
    ) -> Self {
        let scan_record = AlignmentRecord::from(coverage, intervals, None);
        Self {
            id: None,
            index: None,
            aligner: None,
            assembler: None,
            reference: scan_record.reference,
            reference_length: scan_record.length,
            scan_regions: scan_record.regions,
            scan_reads: scan_record.reads,
            scan_alignments: scan_record.alignments,
            scan_bases_covered: scan_record.bases,
            scan_coverage: scan_record.coverage*100.,
            remap_regions: None,
            remap_reads: None,
            remap_alignments: None,
            remap_bases_covered: None,
            remap_coverage: None,
            consensus_length: None,
            consensus_missing: None,
            consensus_completeness: None,
            taxid: None,
            name: None,
            segment: None,
            reference_description: scan_record.description
        }

    }

}



#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct AlignmentRecord {
    pub reference: String,
    pub regions: u64,
    pub reads: u64,
    pub alignments: u64,
    pub bases: u64,
    pub length: u64,
    pub coverage: f64,
    pub description: String,
    pub intervals: String,
    pub rpm: Option<f64>
}

impl AlignmentRecord {
    pub fn from(coverage: Coverage, intervals: bool, total_reads: Option<u64>) -> Self {
        
        Self {
            reference: coverage.reference,
            regions: coverage.regions,
            reads: coverage.reads,
            alignments: coverage.alignments,
            bases: coverage.bases,
            length: coverage.length,
            coverage: coverage.coverage,
            description: coverage.description,
            intervals: if intervals { 
                coverage.intervals
                    .iter()
                    .map(|tag| {
                        format!("{}:{}:{}", tag.start, tag.stop, tag.val)
                    })
                    .join(" ")
                } else { 
                    String::new() 
                },
            rpm: match total_reads {
                Some(total) => Some(total as f64 / coverage.alignments as f64),  // reads doesn't cover paired-end (unique read ids aligned)
                None => None
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct MatchedRecords {
    scan: AlignmentRecord,
    remap: Option<AlignmentRecord>
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VircovSummary {
    pub records: Vec<VircovRecord>
}
impl VircovSummary {
    pub fn new(
        index: Option<PathBuf>, 
        aligner: Option<Aligner>,
        assembler: Option<ConsensusAssembler>, 
        scan_records: Vec<AlignmentRecord>,
        remap_records: Vec<AlignmentRecord>, 
        consensus_records: Vec<ConsensusRecord>, 
        annotation_options: &AnnotationConfig,
    )-> Result<Self, VircovError> {

        // Group the records by reference sequence identifier
        let mut reference_records = Vec::new();

        for scan_record in scan_records {

            let scan_annotation = Annotation::from(&scan_record.description, &annotation_options);

            // We can have multiple results with the same reference  if the reference was segmented. We therefore extract the 
            // segment description from the record annotations and append the extracted segment

            let matching_records: Vec<&AlignmentRecord> = remap_records.iter().filter(|remap_record| {
                let mut scan_id = scan_record.reference.clone();
                let mut remap_id = remap_record.reference.clone();
                if let Some(segment) =  scan_annotation.segment.clone() {
                    scan_id.push_str(&segment);
                    remap_id.push_str(&segment);
                };
                scan_id == remap_id
            }).collect();

            let matched_records = match matching_records.len() {
                0 => MatchedRecords { scan: scan_record, remap: None },
                1 => MatchedRecords { scan: scan_record, remap: Some(matching_records[0].clone()) },
                _ => {
                    log::error!("Found multiple matching records for remap reference with identifier: {}", &scan_record.reference);
                    return Err(VircovError::RemapMatchNotIdentifiable(scan_record.reference.clone()))
                }
            };

            reference_records.push(matched_records);

        }
        

        let mut summary_records = Vec::new();
        for matched_records in reference_records {
            
            let annotation = Annotation::from(&matched_records.scan.description, &annotation_options);

            let consensus_record = match &matched_records.remap {
                None => None,
                Some(remap_record) => {
                    
                    let matching_records: Vec<&ConsensusRecord> = consensus_records.iter().filter(|record | {
                        let mut consensus_id = record.id.clone();
                        let mut remap_id = remap_record.reference.clone();
                        if let Some(segment) = annotation.segment.clone() {
                            consensus_id.push_str(&segment);
                            remap_id.push_str(&segment);
                        };
                        consensus_id == remap_id
                    }).collect();
        

                    match matching_records.len() {
                        0 => None,
                        1 => Some(matching_records[0].clone()),
                        _ => {
                            log::error!("Found multiple matching records for consensus reference with identifier: {}", &remap_record.reference);
                            return Err(VircovError::ConsensusSequenceMatchNotIdentifiable(remap_record.reference.clone()))
                        }
                    }
                }
            };


            summary_records.push(
                VircovRecord::from(
                    index.clone(),
                    aligner.clone(), 
                    assembler.clone(),
                    annotation, 
                    matched_records.scan, 
                    matched_records.remap,
                    consensus_record
                )?
            )
        }

        summary_records.sort_by(|a, b| b.scan_reads.cmp(&a.scan_reads));

        Ok(Self {
            records: summary_records
        })
    }
    pub fn from_tsv(input: &PathBuf, header: bool) -> Result<Self, VircovError> {

        let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').has_headers(header).from_path(&input)?;
        let mut records = Vec::new();
        for rec in reader.deserialize() {
            records.push(rec?);
        }

        Ok(Self { records })

    }
    pub fn write_tsv(&self, output: &PathBuf, header: bool) -> Result<(), VircovError> {

        let mut writer = csv::WriterBuilder::new().delimiter(b'\t').has_headers(header).from_path(&output)?;
        for rec in &self.records {
            writer.serialize(rec)?;
        }
        writer.flush()?;
        Ok(())

    }
    pub fn concatenate(input: &Vec<PathBuf>, output: &PathBuf, min_completeness: f64, file_id: bool) -> Result<(), VircovError> {

        let mut records = Vec::new();
        for file in input {
            let filename = get_file_component(&file, FileComponent::FileStem)?;
            let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').has_headers(true).from_path(&file).unwrap();
            for rec in reader.deserialize() {
                let mut record: VircovRecord = rec?;

                if file_id {
                    record.id = Some(filename.clone());
                }

                if let Some(completeness) = record.consensus_completeness {
                    if completeness >= min_completeness {
                        records.push(record)
                    }
                }
                
            }
        }
        records.sort_by(|a, b| {
            b.consensus_completeness.partial_cmp(&a.consensus_completeness)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let mut writer = csv::WriterBuilder::new().delimiter(b'\t').has_headers(true).from_path(&output).unwrap();
        for record in records {
            writer.serialize(record)?;
        }
        writer.flush()?;

        Ok(())

    }
    pub fn print_table(&self, sort: bool) {

        let records = if sort {
            let mut records = self.records.clone();

            records.sort_by(|a, b| {
                b.remap_coverage.partial_cmp(&a.remap_coverage)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });

            records
        } else {
            self.records.clone()
        };
        
        let mut table = Table::new(&records);
    
        table.modify(
            Columns::new(7..), 
            Width::wrap(32).keep_words()
        ).with(
            Style::modern()
        );

        eprintln!("{}", table);
}
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SubtypeConfig {
}

/// Helper function to read lines from a file into a vector of strings.
pub fn read_lines_to_vec(filename: &PathBuf) -> Result<Vec<String>> {
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
pub fn get_supported_subtypes() -> HashMap<&'static str, IndexMap<&'static str, Vec<&'static str>>> {
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
