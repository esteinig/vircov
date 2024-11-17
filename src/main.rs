use crate::terminal::{App, Commands};
use crate::utils::init_logger;
use crate::subtype::SubtypeDatabase;
use crate::subtype::SubtypeSummary;
use crate::vircov::{Vircov, VircovConfig, VircovSummary};
use crate::annotation::{AnnotationConfig, AnnotationPreset, DatabaseAnnotation};

use rayon::prelude::*;
use anyhow::Result;
use clap::Parser;
use terminal::ToolsCommands;
use vircov::{get_supported_subtypes, read_lines_to_vec};


mod alignment;
mod covplot;
mod error;
mod subtype;
mod terminal;
mod utils;
mod vircov;
mod consensus;
mod annotation;
mod haplotype;

/// Vircov application
///
/// Run the application from arguments provided
/// by the command-line interface
fn main() -> Result<()> {

    init_logger();

    let terminal = App::parse();

    match &terminal.command {
        Commands::Run(args) => {

            let config = VircovConfig::from_run_args(args)?;
            let vircov = Vircov::from(config)?;

            vircov.run(
                &args.output, 
                args.remap_parallel, 
                args.remap_threads, 
                !args.consensus_disabled, 
                args.haplotype,
                args.keep, 
                args.table,
                args.select_by.clone(),
                None, 
                None,
                args.remap_args.clone(),
                args.remap_filter_args.clone(),
                args.remap_exclude_bins.clone(),
                args.include_scans
            )?;

        },
        Commands::Coverage(args) => {
            
            let config = VircovConfig::from_coverage_args(args)?;
            let vircov = Vircov::from(config)?;

            vircov.run_coverage(
                &args.output, 
                args.intervals, 
                args.zero, 
                args.table,
                args.read_id.clone(),
                None, 
                None
            )?;
            
        }
        Commands::Abundance(args) => {
            unimplemented!("Abundance estimates have not been implemented yet")
        },
        Commands::Subtype(args) => {

            log::info!("Preparing database for subtyping...");
            let subtype_db = SubtypeDatabase::from(
                &args.database, 
                &args.workdir
            )?;
            
            log::info!("Database prepared, computing identity and coverage metrics...");
            let collected: Vec<SubtypeSummary> = args.input.clone()
                .into_par_iter()
                .flat_map(|fasta| {
                    subtype_db.subtype(
                        &fasta, 
                        args.min_cov, 
                        args.min_cov_aa, 
                        args.min_cov_prot, 
                        None, 
                        args.threads
                    )
                        .unwrap_or_else(|err| {
                            log::error!("Error processing input file: {}", err);
                            Vec::new()
                        })
                })
                .collect();
        

            log::info!("Creating ranked genomic neighbor typing metrics...");
            subtype_db.create_ranked_tables(
                &args.output, 
                &collected, 
                &args.metric, 
                args.ranks, 
                true, 
                args.with_genotype, 
                subtype_db.protein
            )?;
            
            if !args.keep {
                subtype_db.remove_workdir()?;
            };

        }
        Commands::Tools( subcommand) => {
            match subcommand {
                ToolsCommands::ProcessNcbi(args) => {
            

                    // Messy but working in a basic way. Need better structure for the processor
                    let subtypes = get_supported_subtypes();
                    
                    let virus_subtypes = match subtypes.get(&args.virus.as_str()) {
                        Some(subtype_map) => subtype_map,
                        None => unimplemented!("No entry found for the given virus key."),
                    };
        
                    subtype::process_ncbi_genotypes(&args.input, &args.output, &virus_subtypes)?;
                }
                ToolsCommands::ProcessGisaid(args) => {
                    
                    let db_file = args.outdir.join("db.csv");
                    let db_nuc_file = args.outdir.join("db_nuc.fasta");
        
                    subtype::process_gisaid_genotypes(&args.clades, &args.fasta, db_nuc_file, db_file, args.segment.as_deref())?;
                }
                ToolsCommands::FilterDatabase(args) => {
                    
                    let accessions = match (&args.accessions, &args.accession_file){
                        (_, Some(file)) => read_lines_to_vec(&file).map_err(|_| subtype::SubtypeDatabaseError::AccessionFileError)?,
                        (Some(accessions), _) => accessions.clone(),
                        (None, None) => Vec::new()
                    };
        
                    let (fasta_out, meta_out) = match (&args.output_fasta, &args.output_genotypes) {
                        (Some(fasta), None) => (fasta.clone(), args.genotypes.with_extension("_filtered.csv")),
                        (Some(fasta), Some(meta)) => (fasta.clone(), meta.clone()),
                        (None, Some(meta)) => (args.fasta.with_extension("_filtered.fasta"), meta.clone()),
                        (None, None) => (args.fasta.with_extension("_filtered.fasta"), args.genotypes.with_extension("_filtered.csv"))
                    };
                    
                    subtype::filter_database(
                        args.genotypes.clone(), 
                        args.fasta.clone(), 
                        meta_out, 
                        fasta_out, 
                        args.min_length, 
                        args.remove_duplicates, 
                        accessions
                    )?;
                }
                ToolsCommands::ValidateGenotypes(args) => {
                    subtype::validate_genotypes( &args.genotypes, &args.fasta)?;
                }
                ToolsCommands::ConcatOutput(args) => {
                    VircovSummary::concatenate(&args.input, &args.output, args.min_completeness, args.file_id, args.file_dir)?;
                },
                ToolsCommands::FilterOutput(args) => {
                    VircovSummary::filter_table(
                        &args.input, 
                        &args.output, 
                        args.id.clone(), 
                        args.min_consensus_completeness, 
                        args.min_consensus_coverage_mapq,
                        args.min_remap_coverage, 
                        args.min_remap_depth_coverage, 
                        args.min_scan_alignments, 
                        args.bin.clone()
                    )?;
                },

                ToolsCommands::FilterSample(args) => {
                    VircovSummary::filter_samples(
                        &args.input, 
                        &args.output, 
                        args.bin.clone(),
                        args.exclude_bin.clone()
                    )?;
                },
                ToolsCommands::AnnotateDatabase(args) => {

                    let dba = DatabaseAnnotation::new(
                        &args.annotations, 
                        AnnotationConfig::from_preset(args.preset.clone())
                    )?;

                    dba.annotate(
                        &args.fasta, 
                        &args.output, 
                        args.skipped.clone()
                    )?;

                }
            }
        }
        
    }

    Ok(())
}
