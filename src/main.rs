use anyhow::Result;
use structopt::StructOpt;

use crate::align::{ReadAlignment, ReadAlignmentError};
use crate::cli::Cli;
use crate::cli::Commands::Vircov;

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
        } => {
            let mut _alignment = ReadAlignment::new(fasta, exclude)?;

            let alignment =
                _alignment.from(alignment, min_len, min_cov, min_mapq, alignment_format)?;

            let data = alignment
                .coverage_statistics(regions, seq_len, coverage, reads, &group_by, verbose)?;

            match group_by {
                None => {
                    alignment.to_output(&data, table, read_ids)?;
                }
                Some(group_field) => {
                    match alignment.target_sequences {
                        None => return Err(ReadAlignmentError::GroupSequenceError()),
                        Some(_) => {
                            match covplot {
                                true => return Err(ReadAlignmentError::GroupCovPlotError()),
                                false => {
                                    // If reference sequences have been provided, continue with grouping outputs
                                    let grouped_data = alignment.group_output(
                                        &data,
                                        group_regions,
                                        group_coverage,
                                        group_field,
                                        group_sep,
                                    )?;
                                    alignment.to_output(&grouped_data, table, read_ids)?;
                                }
                            };
                        }
                    }
                }
            };

            match covplot {
                true => alignment.coverage_plots(&data, width)?,
                false => {}
            }
        }
    }
    Ok(())
}
