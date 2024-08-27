use crate::error::VircovError;
use crate::terminal::{App, Commands};
use crate::utils::init_logger;
use crate::subtype::SubtypeDatabase;
use crate::subtype::SubtypeSummary;
use crate::vircov::{Vircov, VircovConfig, VircovSummary};

use rayon::prelude::*;
use anyhow::Result;
use clap::Parser;
use indexmap::IndexMap;
use std::collections::HashMap;
use std::path::PathBuf;
use std::io::{BufRead, BufReader};

mod alignment;
mod covplot;
mod error;
mod subtype;
mod terminal;
mod utils;
mod vircov;
mod consensus;

/// Vircov application
///
/// Run the application from arguments provided
/// by the command-line interface
#[cfg(not(tarpaulin_include))]
fn main() -> Result<(), VircovError> {


    init_logger();

    let terminal = App::parse();

    match &terminal.command {
        Commands::Run(args) => {

            let config = VircovConfig::from_run_args(args)?;
            let vircov = Vircov::from(config)?;

            vircov.run(
                &args.output, 
                args.parallel, 
                args.remap_threads, 
                !args.no_consensus, 
                args.keep, 
                args.table,
                None, 
                None
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
        Commands::Abundance(args) => {},
        Commands::Subtype(args) => {

            let subtype_db = SubtypeDatabase::from(&args.database, &args.workdir)?;
            
            
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
        

            subtype_db.create_ranked_tables(&args.output, &collected, &args.metric, args.ranks, true, args.with_genotype, subtype_db.protein)?;
            
            if !args.keep {
                subtype_db.remove_workdir()?;
            };

        }
        Commands::ProcessNcbi(args) => {
            

            // Messy but working in a basic way. Need better structure for the processor
            let subtypes = HashMap::from([
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
            ]);
            
            let virus_subtypes = match subtypes.get(&args.virus.as_str()) {
                Some(subtype_map) => subtype_map,
                None => unimplemented!("No entry found for the given virus key."),
            };

            subtype::process_ncbi_genotypes(&args.input, &args.output, &virus_subtypes)?;
        }
        Commands::ProcessGisaid(args) => {
            
            let db_file = args.outdir.join("db.csv");
            let db_nuc_file = args.outdir.join("db_nuc.fasta");

            subtype::process_gisaid_genotypes(&args.clades, &args.fasta, db_nuc_file, db_file, args.segment.as_deref())?;
        }
        Commands::FilterDatabase(args) => {
            
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
        Commands::ValidateGenotypes(args) => {
            subtype::validate_genotypes( &args.genotypes, &args.fasta)?;
        }
        Commands::Concat(args) => {
            VircovSummary::concatenate(&args.input, &args.output, args.min_completeness, args.file_id)?;
        }
    }

    Ok(())
}

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