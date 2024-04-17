use crate::align::{ReadAlignment, ReadAlignmentError};
use crate::error::VircovError;
use crate::terminal::{App, Commands};
use crate::utils::init_logger;
use crate::subtype::SubtypeDatabase;
use crate::subtype::SubtypeSummary;

use rayon::prelude::*;
use anyhow::Result;
use clap::Parser;
use indexmap::IndexMap;
use std::collections::HashMap;

mod align;
mod covplot;
mod error;
mod subtype;
mod terminal;
mod utils;

/// Vircov application
///
/// Run the application from arguments provided
/// by the command-line interface
#[cfg(not(tarpaulin_include))]
fn main() -> Result<(), VircovError> {



    init_logger();

    let terminal = App::parse();

    match &terminal.command {
        Commands::Coverage(args) => {
            let verbose = match args.group_select_split {
                Some(_) => 2, // for group refseq selection we need the tags
                None => args.verbose,
            };

            let mut align = ReadAlignment::new(&args.fasta, &args.exclude)?;

            let align = align.read(
                args.alignment.clone(),
                args.min_len,
                args.min_cov,
                args.min_mapq,
                args.alignment_format.clone(),
            )?;

            let mut data = align.coverage_statistics(
                args.regions,
                args.seq_len,
                args.coverage,
                args.regions_coverage,
                args.reads,
                args.aligned,
                &args.group_by,
                verbose,
                args.zero,
            )?;

            match &args.group_by {
                None => {
                    if args.group_select_split.is_some() {
                        return Err(VircovError::ReadAlignment(
                            ReadAlignmentError::GroupSelectSplitError,
                        ));
                    };

                    align.to_output(
                        &mut data,
                        None,
                        args.table,
                        args.header,
                        args.group_sep.clone(),
                        args.read_ids.clone(),
                        args.read_ids_split.clone(),
                        None,
                        None,
                        false,
                        None,
                        args.segment_field.clone(),
                        args.segment_field_nan.clone(),
                    )?;
                }
                Some(group_field) => {
                    match align.target_sequences {
                        None => {
                            return Err(VircovError::ReadAlignment(
                                ReadAlignmentError::GroupSequenceError,
                            ))
                        }
                        Some(_) => {
                            match args.covplot {
                                true => {
                                    return Err(VircovError::ReadAlignment(
                                        ReadAlignmentError::GroupCovPlotError,
                                    ))
                                }
                                false => {
                                    // Make a clone of the original alignment data
                                    let ungrouped_data = Some(data.clone());

                                    // If reference sequences have been provided, continue with grouping outputs
                                    let mut grouped_data = align.group_output(
                                        &data,
                                        args.group_regions,
                                        args.group_coverage,
                                        args.group_aligned,
                                        args.group_reads,
                                        group_field.to_string(),
                                        args.group_sep.clone(),
                                    )?;
                                    align.to_output(
                                        &mut grouped_data,
                                        ungrouped_data,
                                        args.table,
                                        args.header,
                                        args.group_sep.clone(),
                                        args.read_ids.clone(),
                                        args.read_ids_split.clone(),
                                        args.group_select_by.clone(),
                                        args.group_select_split.clone(),
                                        args.group_select_order,
                                        args.group_select_data.clone(),
                                        args.segment_field.clone(),
                                        args.segment_field_nan.clone(),
                                    )?;
                                }
                            };
                        }
                    }
                }
            };

            match args.covplot {
                true => align.coverage_plots(&data, args.width)?,
                false => {}
            }
        }
        Commands::Subtype(args) => {

            let subtype_db = SubtypeDatabase::from(&args.database, &args.workdir)?;
            
            
            log::info!("Database prepared, computing identity and coverage metrics...");
            let collected: Vec<SubtypeSummary> = args.input.clone()
                .into_par_iter()
                .flat_map(|fasta| {
                    subtype_db.subtype(&fasta, args.min_cov, args.min_cov_aa, args.min_cov_prot, None, args.threads)
                        .unwrap_or_else(|err| {
                            log::error!("Error processing input file: {}", err);
                            Vec::new()
                        })
                })
                .collect();
        

            subtype_db.create_ranked_tables(&args.output, &collected, &args.metric, args.ranks, true, false)?;
            
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
                    ("HMPV-A", Vec::from(["/A", "type A", "A/"])),
                    ("HMPV-B", Vec::from(["/B", "type B", "B/"])),
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
    }

    Ok(())
}
