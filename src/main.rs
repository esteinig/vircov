use crate::align::{ReadAlignment, ReadAlignmentError};
use crate::cli::Cli;
use anyhow::Result;
use structopt::StructOpt;

mod align;
mod cli;
mod covplot;
mod utils;

/// Vircov application
///
/// Run the application from arguments provided
/// by the command line interface
#[cfg(not(tarpaulin_include))]
fn main() -> Result<(), ReadAlignmentError> {
    let args = Cli::from_args();

    let verbose = match args.group_select_split {
        Some(_) => 2, // for group refseq selection we need the tags
        None => args.verbose,
    };

    let mut align = ReadAlignment::new(&args.fasta, &args.exclude)?;

    let align = align.read(
        args.alignment,
        args.min_len,
        args.min_cov,
        args.min_mapq,
        args.alignment_format,
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

    match args.group_by {
        None => {
            if args.group_select_split.is_some() {
                return Err(ReadAlignmentError::GroupSelectSplitError);
            };

            align.to_output(
                &mut data,
                None,
                args.table,
                args.header,
                args.group_sep,
                args.read_ids,
                args.read_ids_split,
                None,
                None,
                false,
                None,
                args.segment_field,
                args.segment_field_nan,
            )?;
        }
        Some(group_field) => {
            match align.target_sequences {
                None => return Err(ReadAlignmentError::GroupSequenceError),
                Some(_) => {
                    match args.covplot {
                        true => return Err(ReadAlignmentError::GroupCovPlotError),
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
                                group_field,
                                args.group_sep.clone(),
                            )?;
                            align.to_output(
                                &mut grouped_data,
                                ungrouped_data,
                                args.table,
                                args.header,
                                args.group_sep.clone(),
                                args.read_ids,
                                args.read_ids_split,
                                args.group_select_by,
                                args.group_select_split,
                                args.group_select_order,
                                args.group_select_data,
                                args.segment_field,
                                args.segment_field_nan,
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

    Ok(())
}
