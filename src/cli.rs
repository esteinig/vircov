use std::ffi::{OsStr, OsString};
use std::path::PathBuf;
use structopt::StructOpt;

// Vircov command-line client
#[derive(Debug, StructOpt)]
#[structopt(name = "vircov")]
pub struct Cli {
    #[structopt(subcommand)]
    pub commands: Commands,
}

#[derive(Debug, StructOpt)]
pub enum Commands {
    /// Flag positive coverage distribution from read alignment
    Flag {
        /// Read alignment in PAF format
        #[structopt(
            short = "p",
            long= "paf", 
            parse(try_from_os_str = check_file_exists)
        )]
        paf: PathBuf,
        /// Minimum query alignment length
        #[structopt(short = "l", long = "min_len", default_value = "0")]
        min_len: u64,
        /// Minimum query alignment coverage
        #[structopt(short = "c", long = "min_cov", default_value = "0")]
        min_cov: f64,
        /// Minimum mapping quality
        #[structopt(short = "q", long = "min_mapq", default_value = "0")]
        min_mapq: u8,
    },
}

fn check_file_exists(file: &OsStr) -> Result<PathBuf, OsString> {
    let path = PathBuf::from(file);
    if path.exists() {
        Ok(path)
    } else {
        Err(OsString::from(format!("{:?} does not exist", path)))
    }
}
