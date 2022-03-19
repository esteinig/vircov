use anyhow::Result;
use structopt::StructOpt;

use crate::cli::Cli;
use crate::alignment::Alignment;

mod cli;
mod covplot;
mod alignment;

/// Vircov application
///
/// Run the application from arguments provided
/// by the command line interface
#[cfg(not(tarpaulin_include))]
fn main() -> Result<()> {
    let args = Cli::from_args();

    let paf = Alignment::from_paf(
        args.path,
        args.fasta,
        args.min_len,
        args.min_cov,
        args.min_mapq,
    )?;

    let data = paf.coverage_statistics(args.regions, args.seq_len, args.verbose)?;
    paf.to_console(&data, args.table)?;

    match args.covplot {
        true => paf.coverage_plots(&data, args.width)?,
        false => {}
    }

    Ok(())
}
