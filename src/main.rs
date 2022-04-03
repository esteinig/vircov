use anyhow::Result;
use structopt::StructOpt;

use crate::align::{ReadAlignment, ReadAlignmentError};
use crate::cli::Cli;

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

    let alignment = match args.paf {
        Some(path) => {
            ReadAlignment::from_paf(path, args.fasta, args.min_len, args.min_cov, args.min_mapq)?
        }
        None => match args.bam {
            None => return Err(ReadAlignmentError::FileInputError()),
            Some(path) => ReadAlignment::from_bam(
                path,
                args.fasta,
                args.min_len,
                args.min_cov,
                args.min_mapq,
            )?,
        },
    };

    let data = alignment.coverage_statistics(
        args.regions,
        args.seq_len,
        args.coverage,
        &args.group_by,
        args.verbose,
    )?;

    match args.group_by {
        None => {
            alignment.to_console(&data, args.table)?;
        }
        Some(group_field) => {
            match alignment.target_sequences {
                None => return Err(ReadAlignmentError::GroupSequenceError()),
                Some(_) => {
                    match args.covplot {
                        true => return Err(ReadAlignmentError::GroupCovPlotError()),
                        false => {
                            // If reference sequences have been provided, continue with grouping outputs
                            let grouped_data = alignment.group_output(
                                &data,
                                args.regions,
                                args.coverage,
                                group_field,
                                args.group_sep,
                            )?;
                            alignment.to_console(&grouped_data, args.table)?;
                        }
                    };
                }
            }
        }
    };

    match args.covplot {
        true => alignment.coverage_plots(&data, args.width)?,
        false => {}
    }

    Ok(())
}
