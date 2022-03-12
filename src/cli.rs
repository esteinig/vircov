use std::ffi::{OsStr, OsString};
use std::path::PathBuf;
use structopt::StructOpt;

// Vircov command-line client
#[derive(Debug, StructOpt)]
#[structopt(name = "vircov")]
pub struct Cli {
    /// Alignment in PAF format
    #[structopt(
        short = "p",
        long = "paf", 
        parse(try_from_os_str = check_file_exists)
    )]
    pub paf: PathBuf,
    /// Reference sequences in FASTA format
    #[structopt(
        short = "f",
        long = "fasta", 
        parse(try_from_os_str = check_file_exists)
    )]
    pub fasta: Option<PathBuf>,
    /// Minimum lenmgth of the aligned query sequence
    #[structopt(short = "l", long = "min-len", default_value = "0")]
    pub min_qaln_len: u64,
    /// Minimum coverage of the aligned query sequence
    #[structopt(short = "c", long = "min-cov", default_value = "0")]
    pub min_qaln_cov: f64,
    /// Minimum mapping quality of the alignment
    #[structopt(short = "q", long = "min-mapq", default_value = "0")]
    pub min_mapq: u8,
}

fn check_file_exists(file: &OsStr) -> Result<PathBuf, OsString> {
    let path = PathBuf::from(file);
    if path.exists() {
        Ok(path)
    } else {
        Err(OsString::from(format!("{:?} does not exist", path)))
    }
}
