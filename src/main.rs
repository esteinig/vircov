use anyhow::Result;
use structopt::StructOpt;

use crate::cli::Cli;
use crate::paf::PafFile;
use crate::covplot::CovPlot;

mod cli;
mod paf;
mod covplot;

/// Vircov application
///
/// Run the application from arguments provided
/// by the command line interface
fn main() -> Result<()> {
    env_logger::init();

    let args = Cli::from_args();

    let paf = PafFile::from(
        args.paf,
        args.fasta,
        args.min_len,
        args.min_cov,
        args.min_mapq,
    )?;
    
    paf.target_coverage_distribution(args.verbose)?;
    paf.target_coverage_plots(100)?;
    
    Ok(())
}
