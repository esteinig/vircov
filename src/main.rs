use std::ops::DerefMut;

use crate::align::{ReadAlignment, ReadAlignmentError};
use crate::cli::Cli;
use crate::cli::Commands::Vircov;
use anyhow::Result;
use structopt::StructOpt;

mod align;
mod cli;
mod covplot;

/// Vircov application
///
/// Run the application from arguments provided
/// by the command line interface
#[cfg(not(tarpaulin_include))]
fn main() -> Result<(), ReadAlignmentError> {
    let args = Cli::from_args();
    match args.commands {
        Vircov {
            fasta,
            alignment,
            exclude,
            min_len,
            min_cov,
            min_mapq,
            alignment_format,
            regions,
            seq_len,
            coverage,
            group_by,
            reads,
            verbose,
            table,
            read_ids,
            covplot,
            group_regions,
            group_coverage,
            group_sep,
            width,
            read_ids_split,
            group_select_by,
            group_select_split,
        } => {
            let verbose = match group_select_split {
                Some(_) => 2, // for group refseq selection we need the tags
                None => verbose,
            };

            let mut align = ReadAlignment::new(&fasta, &exclude)?;

            let align = align.from(alignment, min_len, min_cov, min_mapq, alignment_format)?;

            let data =
                align.coverage_statistics(regions, seq_len, coverage, reads, &group_by, verbose)?;

            match group_by {
                None => {
                    align.to_output(&data, table, read_ids, read_ids_split, None, None)?;
                }
                Some(group_field) => {
                    match align.target_sequences {
                        None => return Err(ReadAlignmentError::GroupSequenceError()),
                        Some(_) => {
                            match covplot {
                                true => return Err(ReadAlignmentError::GroupCovPlotError()),
                                false => {
                                    // If reference sequences have been provided, continue with grouping outputs
                                    let grouped_data = align.group_output(
                                        &data,
                                        group_regions,
                                        group_coverage,
                                        group_field,
                                        group_sep,
                                    )?;
                                    align.to_output(
                                        &grouped_data,
                                        table,
                                        read_ids,
                                        read_ids_split,
                                        group_select_by,
                                        group_select_split,
                                    )?;
                                }
                            };
                        }
                    }
                }
            };

            match covplot {
                true => align.coverage_plots(&data, width)?,
                false => {}
            }
        }
    }

    Ok(())
}
