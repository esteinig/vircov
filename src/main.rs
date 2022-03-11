use anyhow::Result;
use structopt::StructOpt;

use crate::cli::Cli;
use crate::cli::Commands::Flag;
use crate::paf::PafAlignment;
use crate::vircov::Vircov;

mod cli;
mod paf;
mod vircov;

/// Vircov application
///
/// Run the application from arguments provided
/// by the command line interface
fn main() -> Result<()> {
    env_logger::init();

    let args = Cli::from_args();
    let vircov = Vircov::new();

    match args.commands {
        Flag {
            paf,
            min_len,
            min_cov,
            min_mapq,
        } => {
            let paf_align = PafAlignment::from(&paf, min_len, min_cov, min_mapq)?;
        }
    }

    Ok(())
}
