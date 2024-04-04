use crate::align::{ReadAlignment, ReadAlignmentError};
use crate::error::VircovError;
use crate::terminal::{App, Commands};
use crate::utils::init_logger;

use anyhow::Result;
use clap::Parser;
use subtype::SubtypeDatabase;

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
    use std::fs::create_dir_all;

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
            /* Outline of initial basic subtyping for consensus genomes

            Subtyping based on best matches (nucleotide, protein) against reference database,
            auto-curated with Cipher for specific panels of viruses e.g. respiratory

            1. Cipher NCBI virus database construction with headers

                Create a FASTA database of the desired combination of sequences from NCBI Viruses,
                using nucleotide sequences (default) and protein sequences (appropriately renamed),
                the create DIAMOND and BLAST reference indices for each - on the fly in Vircov or
                in Cipher?

            2. Implement DIAMOND and BLAST searches against the reference indices.

            3. Adopt BLAST ANI and DIAMON AAI distance measures and summarise as table.

            4. Maybe use kNN or k-MNN to define nearest neighbors?

            5. Maybe output graphs for Netview?

            6. Phylogeny - can of worms...

            7. Taxid database selection or leave it to workflow?

            */

            let subtype_db = SubtypeDatabase::from(&args.database, &args.outdir)?;

            for fasta in &args.input {
                subtype_db.subtype(fasta, &args.outdir, None, args.threads)?;
            }
        }
    }

    Ok(())
}
