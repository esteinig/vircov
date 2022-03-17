use anyhow::Result;
use structopt::StructOpt;

use crate::cli::Cli;
use crate::paf::PafAlignment;

mod cli;
mod covplot;
mod paf;

/// Vircov application
///
/// Run the application from arguments provided
/// by the command line interface
#[cfg(not(tarpaulin_include))]
fn main() -> Result<()> {
    let args = Cli::from_args();

    let paf = PafAlignment::from(
        args.path,
        args.fasta,
        args.min_len,
        args.min_cov,
        args.min_mapq,
    )?;

    let data = paf.coverage_statistics(args.verbose)?;
    paf.to_console(data, args.table)?;

    match args.covplot {
        true => paf.coverage_plots(args.width)?,
        false => {}
    }

    Ok(())
}
