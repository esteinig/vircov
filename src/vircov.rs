use crate::alignment::{parse_reference_fasta, Aligner, AlignmentFormat, ReferenceCoverage, CoverageBin, DepthCoverageRecord, Preset, ReadAlignment, SelectHighest, VircovAligner};
use crate::annotation::{Annotation, AnnotationConfig, AnnotationPreset};
use crate::consensus::{ConsensusAssembler, ConsensusRecord, VircovConsensus};
use crate::error::VircovError;
use crate::haplotype::{Haplotyper, VircovHaplotype};
use crate::subtype::concatenate_fasta_files;
use crate::terminal::{CoverageArgs, RunArgs};
use crate::utils::{get_file_component, write_tsv, FileComponent};

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
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::{create_dir_all, remove_dir_all, remove_file};
use std::path::PathBuf;
use std::io::{BufRead, BufReader};
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};

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
        remap_threads: usize, 
        consensus: bool, 
        haplotype: bool,
        keep: bool,
        table: bool,
        select_by: SelectHighest,
        alignment: Option<PathBuf>,
        alignment_format: Option<AlignmentFormat>,
        remap_args: Option<String>,
        remap_filter_args: Option<String>,
        remap_exclude_bins: Option<Vec<String>>,
        include_scans: bool
    ) -> Result<(), VircovError> {

        log::info!("Alignment scan against reference database ({})", self.config.alignment.aligner);
        let read_alignment = match alignment {
            Some(path) => self.alignment(&path, alignment_format)?,
            None => self.align()?
        };

        log::info!("Coverage metrics computation from alignment scan");
        let coverage = read_alignment.coverage(false, false)?;

        log::info!(
            "Reference alignment binning using '{}' field", 
            self.config.reference.annotation.bin
        );

        let coverage_bins = self.bin_alignments(
            &coverage, 
            remap_exclude_bins
        )?;

        log::info!("Reference genome selection by highest '{}'", select_by);
        let selections = self.select_bin_reference(
            &coverage_bins, 
            select_by,
            None
        )?;
    
        log::info!("Starting bin remapping and consensus assembly ({})", self.config.alignment.aligner);
        let (consensus, remap, depth) = self.remap_bin_reference(
            &selections, 
            remap_args,
            remap_filter_args,
            &self.config.outdir.clone(), 
            parallel, 
            remap_threads, 
            consensus, 
            haplotype,
            keep
        )?;

        log::info!("Starting bin consensus assembly remapping ({})", self.config.alignment.aligner);
        let consensus_remap = self.remap_bin_consensus(&consensus)?;
        
        let scan = if include_scans {
            coverage
        } else {
            selections.values().flatten().cloned().collect()
        };

        let summary = self.summary(scan, remap, consensus, consensus_remap, depth)?;

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
            &self.config.outdir,
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
    pub fn bin_alignments(
        &self,
        coverage: &[ReferenceCoverage],
        exclude_bins: Option<Vec<String>>
    ) -> Result<Vec<CoverageBin>, VircovError> {

        let mut coverage_bins: BTreeMap<String, Vec<&ReferenceCoverage>> = BTreeMap::new();

        for cov in coverage {
            if let Some(bin) = &cov.bin {

                if let Some(ref exclude_bins) = exclude_bins {
                    if exclude_bins.contains(bin) {
                        continue;
                    }
                }

                match coverage_bins.entry(bin.to_string()) {
                    Entry::Occupied(mut entry) => {
                        entry.get_mut().push(cov);
                    }
                    Entry::Vacant(entry) => {
                        entry.insert(vec![cov]);
                    }
                }
            } else {
                return Err(VircovError::BinAnnotationMissing(cov.reference.to_owned()))
            }
        }

        let mut binned_coverage_all: Vec<CoverageBin> = Vec::new();
        for (bin, coverage) in coverage_bins {
            
            let bin = bin.trim().replace(
                &self.config.reference.annotation.bin_field(),""
            ); // TODO

            let binned_coverage_data = CoverageBin::from_coverage(&bin, coverage)?;

            if binned_coverage_data.total_regions >= self.config.filter.min_grouped_regions
                && binned_coverage_data.mean_coverage >= self.config.filter.min_grouped_mean_coverage
                && binned_coverage_data.total_alignments >= self.config.filter.min_grouped_alignments
                && binned_coverage_data.total_reads >= self.config.filter.min_grouped_reads
            {
                binned_coverage_all.push(binned_coverage_data);
            }
        }

        Ok(binned_coverage_all)
    }
    pub fn extract_bin_reads(
        &self,
        coverage_bins: &Vec<CoverageBin>,
        outdir: &PathBuf,
        threads: usize
    ) -> Result<HashMap<String, Vec<PathBuf>>, VircovError> {  // group, read fastq

        let result = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("Failed to create thread pool")
            .install(|| -> HashMap<String, Vec<PathBuf>> {

            let bin_reads = coverage_bins.par_iter().map(|coverage_bin| -> (String, Vec<PathBuf>) {
                
                log::info!("Extracting reads for remapping of bin: {}", coverage_bin.sanitized_id());

                let output = self.config.alignment.get_output_read_paths(
                    &coverage_bin.sanitized_id(), 
                    "fastq", 
                    outdir
                ).expect(&format!("Failed to get read paths for reference bin remapping: {}", coverage_bin.sanitized_id()));

                coverage_bin.write_group_reads(
                    self.config.alignment.input.clone(), 
                    output.clone()
                ).expect(&format!("Failed to write reads for reference bin remapping: {}", coverage_bin.sanitized_id()));

                // Same keys as the group selection coverage structs for remapping
                (coverage_bin.sanitized_id().clone(), output)
            }).collect::<Vec<_>>();

            HashMap::from_iter(bin_reads)
        });

        Ok(result)
    }
    pub fn select_bin_reference(
        &self,
        coverage_bins: &Vec<CoverageBin>, // If grouped, these are grouped fields
        select_by: SelectHighest,
        outdir: Option<PathBuf>,
    ) -> Result<HashMap<String, Vec<ReferenceCoverage>> , VircovError> {
        
        let refs = self.references
            .as_ref()
            .ok_or(VircovError::BinSequenceError)?;

        if let Some(ref path) = outdir {
            std::fs::create_dir_all(path)?;
        }

        let mut selected_reference_coverage: HashMap<String, Vec<ReferenceCoverage>> = HashMap::new();

        for coverage_bin in coverage_bins {
            let mut grouped_segments: BTreeMap<String, Vec<ReferenceCoverage>> = BTreeMap::new();

            for cov in &coverage_bin.coverage {

                if let Some(segment) = &cov.segment {

                    if segment != &self.config.reference.annotation.segment_nan {
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

                let mut selected_segments: BTreeMap<String, ReferenceCoverage> = BTreeMap::new();

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
                        None => return Err(VircovError::BinSelectReference),
                    }
                }

                if let Some(path) = &outdir {

                    let file_handle = std::fs::File::create(
                        path.join(
                            format!("{}.fasta", coverage_bin.sanitized_id())
                        )
                    )?;
                    let mut writer = noodles::fasta::Writer::new(file_handle);

                    for (_, segment_cov) in selected_segments {
                        let seq = &refs[&segment_cov.reference];
                        writer.write_record(seq)?;
                        
                        selected_reference_coverage.entry(coverage_bin.sanitized_id().clone())
                            .and_modify(|x| x.push(segment_cov.clone()))
                            .or_insert(vec![segment_cov]);
                    }
                } else {
                    for (_, segment_cov) in selected_segments {
                        selected_reference_coverage.entry(coverage_bin.sanitized_id().clone())
                            .and_modify(|x| x.push(segment_cov.clone()))
                            .or_insert(vec![segment_cov]);
                    }
                }

            } else {

                // If not segmented, simply select the best reference sequence using the
                // selection metric and select the sequence from the header fields
                // to write to FASTA

                let max = match select_by {
                    SelectHighest::Reads => coverage_bin.coverage.iter().max_by_key(|x| x.reads),
                    SelectHighest::Coverage => coverage_bin.coverage.iter().max_by_key(|x| OrderedFloat(x.coverage)),
                    SelectHighest::Alignments => coverage_bin.coverage.iter().max_by_key(|x| x.alignments)
                };

                match max {
                    Some(field) => {
                        let seq = &refs[&field.reference];

                        if let Some(path) = &outdir {
                            
                            let file_handle = std::fs::File::create(
                                path.join(
                                    format!("{}.fasta", coverage_bin.sanitized_id())
                                )
                            )?;
                            let mut writer = noodles::fasta::Writer::new(file_handle);
                            writer.write_record(seq)?;
                        }

                        selected_reference_coverage.entry(coverage_bin.sanitized_id().clone())
                            .and_modify(|x| x.push(field.clone()))
                            .or_insert(vec![field.clone()]);
                    }
                    _ => return Err(VircovError::BinSelectReference),
                }
            }
        }

        Ok(selected_reference_coverage)
    }
    pub fn remap_bin_reference(
        &self, 
        selections: &HashMap<String, Vec<ReferenceCoverage>>, 
        remap_args: Option<String>,
        remap_filter_args: Option<String>,
        outdir: &PathBuf,
        parallel: usize, 
        threads: usize,
        consensus: bool,
        haplotype: bool,
        keep: bool
    ) -> Result<(Vec<ConsensusRecord>, Vec<ReferenceCoverage>, Vec<DepthCoverageRecord>), VircovError> {

        let refs = self.references
            .as_ref()
            .ok_or(VircovError::BinSequenceError)?;

        rayon::ThreadPoolBuilder::new()
            .num_threads(parallel)
            .build()
            .expect("Failed to create thread pool")
            .install(|| -> Result<(Vec<ConsensusRecord>, Vec<ReferenceCoverage>, Vec<DepthCoverageRecord>), VircovError> {

            let results: Vec<Result<(Vec<Vec<ConsensusRecord>>, Vec<ReferenceCoverage>, Vec<DepthCoverageRecord>), VircovError>> = selections
                .into_par_iter()
                .map(|(bin, coverage)| -> Result<(Vec<Vec<ConsensusRecord>>, Vec<ReferenceCoverage>, Vec<DepthCoverageRecord>), VircovError> {
                    
                    let remap_id = bin;

                    let bin_read_files = if self.config.alignment.remap_bin_reads {
                        
                        // Recreate the coverage bin here for read set access
                        let coverage_bin = CoverageBin::from_coverage(
                            bin, 
                            coverage.iter().collect()
                        )?;

                        log::info!("Extracting reads for remapping of bin: {}", coverage_bin.sanitized_id());
                        
                        let output = self.config.alignment.get_output_read_paths(
                            &coverage_bin.sanitized_id(), 
                            "fastq", 
                            outdir
                        )?;
        
                        coverage_bin.write_group_reads(
                            self.config.alignment.input.clone(), 
                            output.clone()
                        )?;
                        Some(output)
                    } else {
                        log::info!("All input reads will be used for remapping");
                        None
                    };

                    log::info!("Reference bin '{}' remapping (n = {})", bin, coverage.len());

                    let remap_reference = outdir.join(remap_id).with_extension("fasta");
                    
                    self.extract_remap_sequences(
                        &remap_reference, 
                        coverage, 
                        refs
                    )?;

                    let bam = outdir.join(remap_id).with_extension("bam");

                    let remap_aligner = VircovAligner::from(
                        &self.config.outdir,
                        &self.config.alignment.remap_config(
                            None, // uses remap reference
                            remap_args.clone(),
                            remap_filter_args.clone(),
                            bin_read_files.clone(),
                            Some(bam.clone()), 
                            false,
                            threads
                        ),
                        &self.config.reference.remap_config(
                            &remap_reference
                        ),
                        &self.config.filter,
                    );

                    let remap_coverage: Vec<ReferenceCoverage> = remap_aligner
                        .run_aligner()?
                        .unwrap()
                        .coverage(true, false)?;
                    
                    // Output read sets based on remap coverage
                    
                    let mut bin_read_ids = HashSet::new();
                    for cov in &remap_coverage {
                        for read_id in &cov.read_id {
                            bin_read_ids.insert(read_id);
                        }
                    }

                    write_tsv(
                        &bin_read_ids.iter().collect(), 
                        &self.config.outdir.join(format!("{remap_id}.reads.tsv.gz"))
                    )?;

                    // Remove binned reads after alignment - take up a lot of disk space
                    if let Some(bin_reads) = bin_read_files {
                        for file in bin_reads {
                            remove_file(file).map_err(|e| {
                                log::error!("Error removing bin read file");
                                e
                            })?;
                        }
                    }

                    let mut consensus_records = Vec::new();
                    let mut depth_records = Vec::new();

                    for ref_cov in &remap_coverage {

                        // Reference output must have extension '.fa' matching auto-extension from iVar
                        let (filter_reference, consensus_name, seg_name) = match &ref_cov.segment {
                            Some(segment) => {
                                let seg = self.config.reference.annotation.segment_name_file(segment);
                                (
                                    Some(ref_cov.reference.clone()), 
                                    format!("{bin}.{seg}.consensus.fa"),
                                    seg
                                )
                            }, 
                            None => (
                                None, 
                                format!("{bin}.nan.consensus.fa"),
                                String::from("nan"),
                            )
                        };
                        
                        log::info!("Extracting depth of coverage information for bin '{}'", bin);
                        let depth = outdir.join(format!("{bin}.{seg_name}.depth"));
                        let depth_coverage = remap_aligner.run_depth_coverage(  // runs only if min_depth_coverage is Some and > 0
                            &bam, 
                            &ref_cov.reference,
                            &depth, 
                            self.config.alignment.min_depth_coverage,
                            ref_cov.segment.clone()
                        )?;

                        depth_records.push(depth_coverage);


                        log::info!("Remap coverage for '{bin}': {:.2}", ref_cov.coverage*100.0);

                        if haplotype && ref_cov.coverage >= self.config.filter.min_remap_coverage {

                            let vircov_haplotype = VircovHaplotype::new(
                                    outdir, 
                                    HaplotypeConfig::with_default(
                                    &bam, 
                                    &remap_reference, 
                                    filter_reference.clone(),
                                    self.config.haplotype.haplotyper.clone()  // set from run args
                                )
                            )?;
                            vircov_haplotype.haplotype()?;
                        }

                        if consensus && ref_cov.coverage >= self.config.filter.min_remap_coverage {

                            log::info!("Creating consensus assembly from bin alignment '{bin}': {consensus_name}");
                            let vircov_consensus = VircovConsensus::new(
                                ConsensusConfig::with_default(
                                    &bam.clone(), 
                                    &outdir.join(consensus_name), 
                                    filter_reference.clone(),
                                    Some(format!("{} {}", ref_cov.reference, ref_cov.description.replace("'", ""))),
                                    self.config.consensus.min_quality,
                                    self.config.consensus.min_frequency,
                                    self.config.consensus.min_depth,
                                    self.config.consensus.max_depth
                                )
                            )?;

                            let records = vircov_consensus.assemble()?;

                            consensus_records.push(records)
                        }
                    }

                    if !keep {
                        remove_file(&bam)?;
                        remove_file(&remap_reference)?;

                        let bai = bam.with_extension("bai");
                        if bai.exists() { remove_file(bai)? };
                        
                    }
                    
                    Ok((consensus_records, remap_coverage, depth_records))
                })
                .collect();
            

            let mut remap_data = Vec::new();
            let mut consensus_data = Vec::new();
            let mut depth_data = Vec::new();

            for (i, result) in results.into_iter().enumerate() {

                let (consensus, remap, depth_coverage) = result.map_err(|e| {
                    log::error!("Error occurring here at index: {i}");
                    e
                })?;
                
                for consensus_record in consensus.into_iter().flatten() {
                    consensus_data.push(consensus_record)
                }
                for cov in remap.into_iter() {
                    remap_data.push(cov)
                }
                for depth in depth_coverage.into_iter() {
                    depth_data.push(depth)
                }
            }

            Ok((consensus_data, remap_data, depth_data))
        })

    }

    pub fn remap_bin_consensus(
        &self,
        consensus: &Vec<ConsensusRecord>
    ) -> Result<Vec<ReferenceCoverage>, VircovError> {

        log::info!("Remapping input reads against detected consensus sequences");

        // Remap input reads against the consensus sequences of detected viruses
        let fasta: Vec<PathBuf> = consensus
            .into_iter()
            .filter(|r| r.completeness >= self.config.filter.min_consensus_completeness)
            .map(|r| r.fasta.to_path_buf())
            .collect();

        let consensus_db = self.config.outdir.join("consensus.fasta");
        concatenate_fasta_files(&fasta, &consensus_db)?;

        let bam = self.config.outdir.join("consensus.bam");

        let remap_aligner = VircovAligner::from(
            &self.config.outdir,
            &self.config.alignment.remap_config(
                None, // uses consensus_db reference
                None,
                Some("-F 4".to_string()),
                Some(self.config.alignment.input.clone()),
                Some(bam.clone()), 
                false,
                self.config.alignment.threads
            ),
            &self.config.reference.remap_config(
                &consensus_db
            ),
            &FilterConfig::with_default_mapq(
                if self.config.alignment.aligner == Aligner::Bowtie2 { 40 } else { 60 },
                self.config.filter.min_consensus_completeness
            ),
        );

        let coverage: Vec<ReferenceCoverage> = remap_aligner
            .run_aligner()?
            .unwrap()
            .coverage(true, false)?;
        
        Ok(coverage)

    }
    fn extract_remap_sequences(
        &self, 
        remap_fasta: &PathBuf, 
        coverage: &Vec<ReferenceCoverage>, 
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
        scan: Vec<ReferenceCoverage>, 
        remap: Vec<ReferenceCoverage>,
        consensus: Vec<ConsensusRecord>, 
        consensus_remap: Vec<ReferenceCoverage>,
        depth: Vec<DepthCoverageRecord>
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
            depth,
            consensus,
        consensus_remap
            .iter()
            .map(|cov| AlignmentRecord::from(cov.clone(), false, None))
            .collect(),
            &&self.config.reference.annotation
        )
    }
}

// fn identify_unique_reads(coverages: &Vec<ReferenceCoverage>) -> Vec<(String, Option<String>, usize, usize)> {
//     use std::collections::{HashMap, HashSet};

//     // Step 1: Count occurrences of each read ID across all references
//     let mut read_counts: HashMap<String, usize> = HashMap::new();
//     for coverage in coverages.iter() {
//         for read_id in &coverage.read_id {
//             *read_counts.entry(read_id.clone()).or_insert(0) += 1;
//         }
//     }

//     // Step 2: Identify reads that appear only once across all references
//     let unique_reads: HashSet<String> = read_counts
//         .into_iter()
//         .filter(|&(_, count)| count == 1)
//         .map(|(read_id, _)| read_id)
//         .collect();

//     // Step 3: Create new ReferenceCoverage structs with unique reads
//     let filtered_coverages: Vec<(String, Option<String>, usize, usize)> = coverages
//         .iter()
//         .map(|coverage| {
//             // Retain only read IDs that are unique to this reference
//             let filtered_read_ids: Vec<String> = coverage
//                 .read_id
//                 .iter()
//                 .filter(|read_id| unique_reads.contains(*read_id))
//                 .cloned()
//                 .collect();

//             (coverage.reference.clone(), coverage.bin.clone(), coverage.read_id.len(), filtered_read_ids.len())
//         })
//         .collect();

//     filtered_coverages
// }

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VircovConfig {
    pub outdir: PathBuf,
    pub alignment: AlignmentConfig,
    pub reference: ReferenceConfig,
    pub filter: FilterConfig,
    pub consensus: ConsensusConfig,
    pub haplotype: HaplotypeConfig,
    pub subtype: SubtypeConfig
} 
impl VircovConfig {
    pub fn from_run_args(args: &RunArgs) -> Result<Self, VircovError> {

        let outdir = match &args.workdir {
            None => PathBuf::from("vircov_data"),
            Some(dir) => dir.to_owned()
        };

        if !outdir.exists() && !outdir.is_dir() {
            create_dir_all(&outdir)?
        }

        Ok(Self {
            outdir: outdir.clone(),
            alignment: AlignmentConfig::with_default(
                &args.input, 
                args.index.clone(), 
                args.aligner.clone(), 
                args.preset.clone(), 
                args.secondary,
                Some(outdir.join("scan.bam")),
                args.scan_threads,
                args.scan_args.clone(),
                args.scan_filter_args.clone(),
                !args.remap_all,
                args.min_depth_coverage
            ),
            reference: ReferenceConfig::with_default(
                Some(args.reference.clone()),
                args.annotation_preset.clone()
            ),
            consensus: ConsensusConfig::with_default_from_args(
                args.consensus_max_depth,
                args.consensus_min_quality,
                args.consensus_min_frequency,
                args.consensus_min_depth,
            ),
            filter: FilterConfig::with_default(args.min_remap_coverage, args.consensus_min_completeness),
            haplotype: HaplotypeConfig::from_run_args(args.haplotyper.clone()),
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
                args.index.clone(), 
                args.aligner.clone(), 
                args.preset.clone(), 
                args.secondary,
                Some(outdir.join("scan.bam")),
                args.threads,
                args.args.clone(),
                None,
                false,
                None
            ),
            reference: ReferenceConfig::with_default(
                args.reference.clone(),
                args.annotation_preset.clone()
            ),
            filter: FilterConfig::default(),
            consensus: ConsensusConfig::default(), // not used
            haplotype: HaplotypeConfig::default(), // not used
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
            haplotype: HaplotypeConfig::default(),
            subtype: SubtypeConfig {}
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignmentConfig {
    pub aligner: Aligner,
    pub preset: Option<Preset>,
    pub args: Option<String>,
    pub filter_args: Option<String>,
    pub input: Vec<PathBuf>,
    pub paired_end: bool,
    pub index: Option<PathBuf>,
    pub secondary: bool,
    pub output: Option<PathBuf>,
    pub threads: usize,
    pub remap_bin_reads: bool,
    pub min_depth_coverage: Option<usize>,
}
impl AlignmentConfig {
    pub fn remap_config(
        &self, 
        index: Option<PathBuf>, 
        args: Option<String>,
        filter_args: Option<String>,
        input: Option<Vec<PathBuf>>,
        output: Option<PathBuf>, 
        secondary: bool,
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

        config.args = args;
        config.index = index.clone();
        config.threads = threads;
        config.output = output;
        config.filter_args = filter_args;

        config.secondary = secondary;

        config
    }
    pub fn with_default(
        input: &Vec<PathBuf>, 
        index: Option<PathBuf>, 
        aligner: Aligner, 
        preset: Option<Preset>, 
        secondary: bool, 
        output: Option<PathBuf>, 
        threads: usize,
        args: Option<String>,
        filter_args: Option<String>,
        remap_bin_reads: bool,
        min_depth_coverage: Option<usize>
    ) -> Self {
        Self {
            input: input.clone(),
            index: index.clone(),
            aligner,
            preset,
            output,
            secondary,
            threads,
            args,
            filter_args,
            remap_bin_reads,
            min_depth_coverage,
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
            filter_args: None,
            input: Vec::from([PathBuf::from("smoke_R1.fastq.gz"), PathBuf::from("smoke_R2.fastq.gz")]),
            paired_end: true,
            secondary: false,
            index: Some(PathBuf::from("reference.fasta")),
            output: Some(PathBuf::from("test.sam")),
            threads: 8,
            remap_bin_reads: true,
            min_depth_coverage: Some(10)
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
    pub fn with_default(fasta: Option<PathBuf>, preset: AnnotationPreset) -> Self {
        Self {
            fasta: fasta.clone(),
            annotation: AnnotationConfig::from_preset(preset),
            ..Default::default()
        }
    }
}
impl Default for ReferenceConfig {
    fn default() -> Self {
        Self {
            fasta: None,
            exclude: None,
            annotation: AnnotationConfig::default()
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
    pub min_remap_coverage: f64,
    pub min_consensus_completeness: f64
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
            min_remap_coverage: 0.0,
            min_consensus_completeness: 0.0
        }
    }
}
impl FilterConfig {
    pub fn with_default(
        min_remap_coverage: f64,
        min_consensus_completeness: f64
    ) -> Self {
        Self {
            min_remap_coverage,
            min_consensus_completeness,
            ..Default::default()
        }
    }
    pub fn with_default_mapq(
        min_mapq: u8,
        min_consensus_completeness: f64
    ) -> Self {
        Self {
            min_mapq,
            min_consensus_completeness,
            ..Default::default()
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct HaplotypeConfig {
    pub haplotyper: Haplotyper,
    pub alignment: PathBuf,
    pub fasta: PathBuf,
    pub reference: Option<String>,
    pub threads: usize,
    pub min_var_freq: f64,
    pub snp_count_filter: usize,
    pub snp_density: f64,
    pub ploidy_sensitivity: u8,
    pub min_consensus_quality: usize,
    pub min_consensus_freq: f64,
    pub min_consensus_depth: usize
}
impl HaplotypeConfig {
    pub fn with_default(
        alignment: &PathBuf,
        fasta: &PathBuf,
        reference: Option<String>,
        haplotyper: Haplotyper,
    ) -> Self {
        Self {
            alignment: alignment.to_path_buf(),
            fasta: fasta.to_path_buf(),
            reference: reference.to_owned(),
            haplotyper: haplotyper.clone(),
            ..Default::default()
        }
    }
    pub fn from_run_args(
        haplotyper: Haplotyper
    ) -> Self {
        Self {
            haplotyper,
            ..Default::default()
        }
    }
}
impl Default for HaplotypeConfig {
    fn default() -> Self {
        Self {
            haplotyper: Haplotyper::Floria,
            alignment: PathBuf::from(""),
            fasta: PathBuf::from(""),
            reference: None,            // ivar segment reference subset
            threads: 8,
            min_var_freq: 0.05,         // not used currently
            snp_count_filter: 5,
            snp_density: 0.00001,
            ploidy_sensitivity: 3,
            min_consensus_quality: 10,
            min_consensus_freq: 0.0,
            min_consensus_depth: 1,
        }
    }
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
    pub max_depth: usize,
    pub min_quality: usize,
    pub min_frequency: f64,
    pub min_depth: usize,
    pub missing: String,
}
impl ConsensusConfig {
    pub fn with_default(
        alignment: &PathBuf, 
        output: &PathBuf, 
        reference: Option<String>, 
        header: Option<String>, 
        min_quality: usize, 
        min_frequency: f64, 
        min_depth: usize,
        max_depth: usize
    ) -> Self {
        Self {
            alignment: alignment.clone(),
            output: output.clone(),
            header: header.clone(),
            reference,
            min_quality,
            min_frequency,
            min_depth,
            max_depth,
            ..Default::default()
        }
    }
    pub fn with_default_from_args(max_depth: usize, min_quality: usize, min_frequency: f64, min_depth: usize) -> Self {
        Self {
            min_quality,
            min_frequency,
            min_depth,
            max_depth,
            ..Default::default()
        }
    }
}
impl Default for ConsensusConfig {
    fn default() -> Self {
        Self {
            alignment: PathBuf::from(""),
            output: PathBuf::from("remap.consensus.fa"), // must be .fa for Ivar
            assembler: ConsensusAssembler::Ivar,
            reference: None,
            header: None,
            args: None,
            mpileup: None,
            max_depth: 10000,
            min_quality: 20,
            min_frequency: 0.75,
            min_depth: 10, 
            missing: String::from("N"),

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

pub fn display_option_usize(opt: &Option<usize>) -> String 
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
    #[tabled(rename="Bin")]
    #[tabled(display_with = "display_option_string")]
    pub bin: Option<String>,
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
    #[tabled(display_with = "display_option_f64")]
    pub remap_depth: Option<f64>,
    #[tabled(rename="Coverage (depth >= X)")]
    #[tabled(display_with = "display_option_f64")]
    pub remap_depth_coverage: Option<f64>,
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
    #[tabled(rename="Consensus Alignments")]
    pub consensus_alignments_mapq: Option<u64>,
    #[tabled(skip)]
    #[tabled(rename="Consensus Coverage")]
    #[tabled(display_with = "display_option_u64")]
    pub consensus_coverage_mapq: Option<f64>,
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
        remap_depth_record: Option<DepthCoverageRecord>,
        consensus_remap_record: Option<AlignmentRecord>,
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

        let (mean_depth, mean_depth_coverage_percent) = match remap_depth_record {
            Some(record) => (Some(record.mean_depth), Some(record.mean_depth_coverage*100.0)),
            None => (None, None)
        };

        let (consensus_alignments_mapq, consensus_coverage_mapq) = match consensus_remap_record  {
            Some(record) => (Some(record.alignments), Some(record.coverage*100.0)),
            None => (None, None)
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
            remap_depth: mean_depth,
            remap_depth_coverage: mean_depth_coverage_percent,
            consensus_length,
            consensus_missing,
            consensus_completeness,
            consensus_alignments_mapq,
            consensus_coverage_mapq,
            bin: annotation.bin,
            name: annotation.name,
            segment: annotation.segment,
            reference_description: scan_record.description,
        })
    }
    pub fn from_coverage(
        coverage: ReferenceCoverage,
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
            remap_depth: None,
            remap_depth_coverage: None,
            consensus_length: None,
            consensus_missing: None,
            consensus_completeness: None,
            consensus_alignments_mapq: None,
            consensus_coverage_mapq: None,
            bin: None,
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
    pub fn from(coverage: ReferenceCoverage, intervals: bool, total_reads: Option<u64>) -> Self {
        
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
        depth_records: Vec<DepthCoverageRecord>,
        consensus_records: Vec<ConsensusRecord>, 
        consensus_remap_records: Vec<AlignmentRecord>,
        annotation_options: &AnnotationConfig,
    )-> Result<Self, VircovError> {

        // Group the records by reference sequence identifier
        let mut reference_records = Vec::new();

        for scan_record in scan_records {

            let scan_annotation = Annotation::from(&scan_record.description, &annotation_options);

            // We can have multiple results with the same reference if the reference was segmented. We therefore extract the 
            // segment description from the record annotations and append the extracted segment

            let matching_records: Vec<&AlignmentRecord> = remap_records.iter().filter(|remap_record| {
                let mut scan_id = scan_record.reference.clone();
                let mut remap_id = remap_record.reference.clone();
                if let Some(segment) = scan_annotation.segment.clone() {
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

            let consensus_remap_record = consensus_remap_records.iter().find(|remap| {
                let mut remap_id = remap.reference.clone();
                let mut scan_id = matched_records.scan.reference.clone();
                if let Some(segment) = annotation.segment.clone() {
                    scan_id.push_str(&segment);
                    remap_id.push_str(&segment);
                }
                scan_id == remap_id
            });

            // Find the corresponding depth coverage record
            let depth_record = depth_records.iter().find(|depth_record| {
                let mut scan_id = matched_records.scan.reference.clone();
                if let Some(segment) = annotation.segment.clone() {
                    scan_id.push_str(&segment); // record.id field for depth_record has this done when it is created
                };
                scan_id == depth_record.id
            });

            

            summary_records.push(
                VircovRecord::from(
                    index.clone(),
                    aligner.clone(), 
                    assembler.clone(),
                    annotation, 
                    matched_records.scan, 
                    matched_records.remap,
                    depth_record.cloned(),
                    consensus_remap_record.cloned(),
                    consensus_record,
                )?
            )
        }

        summary_records.sort_by(|a, b| b.scan_alignments.cmp(&a.scan_alignments));

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
    pub fn filter_table(
        input: &PathBuf,
        output: &PathBuf,
        id_filter: Option<Vec<String>>,
        consensus_completeness: Option<f64>,
        consensus_coverage_mapq: Option<f64>,
        remap_coverage: Option<f64>,
        remap_depth_coverage: Option<f64>,
        scan_alignments: Option<u64>,
        bin_filter: Option<String>,
    ) -> Result<VircovSummary, VircovError> {

        // Read records from the input TSV file
        let summary = VircovSummary::from_tsv(input, true)?;
        let mut filtered_records = Vec::new();

        for record in summary.records {

            // Apply each filter condition
            let id_match = id_filter.as_ref().map_or(true, |ids| {
                record.id.as_ref().map_or(false, |rec_id| ids.contains(rec_id))
            });
            let completeness_match = consensus_completeness.map_or(true, |min| {
                record.consensus_completeness.map_or(false, |comp| comp >= min)
            });
            let coverage_match = remap_coverage.map_or(true, |min_cov| {
                record.remap_coverage.map_or(false, |cov| cov >= min_cov)
            });
            let consensus_coverage_mapq_match = consensus_coverage_mapq.map_or(true, |min_cov| {
                record.consensus_coverage_mapq.map_or(false, |cov| cov >= min_cov)
            });

            let coverage_depth_match = remap_depth_coverage.map_or(true, |min_cov| {
                record.remap_depth_coverage.map_or(false, |cov| cov >= min_cov)
            });
            let alignments_match = scan_alignments.map_or(true, |min_align| {
                record.scan_alignments >= min_align
            });
            let bin_match = bin_filter.as_ref().map_or(true, |bin| {
                record.bin.as_ref().map_or(false, |rec_bin| rec_bin == bin)
            });

            // Only include records that meet all criteria
            if id_match && completeness_match && coverage_match && alignments_match && bin_match && coverage_depth_match && consensus_coverage_mapq_match {
                filtered_records.push(record);
            }
        }

        let summary = VircovSummary {
            records: filtered_records,
        };

        VircovSummary::write_tsv(&summary, output, true)?;

        Ok(summary)
    }
    pub fn filter_samples(
        input: &PathBuf,
        output: &PathBuf,
        bin: Option<Vec<String>>,
        exclude_bin: Option<Vec<String>>,
        strict: bool,
    ) -> Result<Self, VircovError> {
        let summary = VircovSummary::from_tsv(input, true)?;
    
        let mut retained_samples = Vec::new();
    
        let target_bins = bin.unwrap_or_default();
        let excluded_bins = exclude_bin.unwrap_or_default();
    
        if !target_bins.is_empty() {
            // Retain all records of samples that match the target bins based on the `strict` flag
            let sample_ids: Vec<String> = summary
                .records
                .iter()
                .filter_map(|record| {
                    if let Some(ref record_bin) = record.bin {
                        if strict {
                            // Strict mode: Retain only samples containing ALL target bins
                            if target_bins.iter().all(|bin| record_bin.contains(bin)) {
                                return record.id.clone();
                            }
                        } else {
                            // Non-strict mode: Retain samples containing ANY target bin
                            if target_bins.iter().any(|bin| record_bin.contains(bin)) {
                                return record.id.clone();
                            }
                        }
                    }
                    None
                })
                .collect();
    
            // Retain all records of these samples
            for record in &summary.records {
                if let Some(ref record_id) = record.id {
                    if sample_ids.contains(record_id) {
                        // Exclude records with excluded bins
                        if let Some(ref record_bin) = record.bin {
                            if excluded_bins.contains(record_bin) {
                                continue;
                            }
                        }
                        retained_samples.push(record.clone());
                    }
                }
            }
        }
    
        let filtered_summary = VircovSummary {
            records: retained_samples,
        };
    
        VircovSummary::write_tsv(&filtered_summary, output, true)?;
    
        Ok(filtered_summary)
    }    

    pub fn concatenate(input: &Vec<PathBuf>, output: &PathBuf, min_completeness: Option<f64>, file_id: bool, file_dir: bool) -> Result<(), VircovError> {

        let mut records = Vec::new();
        for file in input {
            let filename = get_file_component(&file, FileComponent::FileStem)?;
            let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').has_headers(true).from_path(&file).unwrap();
            for rec in reader.deserialize() {
                let mut record: VircovRecord = rec?;

                if file_id {
                    record.id = Some(filename.clone());
                } else if file_dir {
                    let parent_dir = file.parent().unwrap_or(&file).to_path_buf();
                    let dirname = get_file_component(&parent_dir, FileComponent::FileName)?;
                    record.id = Some(dirname)
                }

                if let Some(min_completeness) = min_completeness {
                    if let Some(completeness) = record.consensus_completeness {
                        if completeness >= min_completeness {
                            records.push(record)
                        }
                    }
                } else {
                    records.push(record)
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
