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
    #[structopt(short = "l", long = "min-len", default_value = "50")]
    pub min_len: u64,
    /// Minimum coverage of the aligned query sequence
    #[structopt(short = "c", long = "min-cov", default_value = "0.5")]
    pub min_cov: f64,
    /// Minimum mapping quality of the alignment
    #[structopt(short = "q", long = "min-mapq", default_value = "30")]
    pub min_mapq: u8,
    /// Verbose output statistics  
    ///
    /// Single flag (-v) adds whitespace separated tags in the last column,
    /// corresponding to the number of inferred alignment coverage blocks.
    /// Each tag specifies: start of block, end of block and number of
    /// alignments in the block separated by a semicolon [example: 3450:3500:9]
    #[structopt(short, long, parse(from_occurrences = parse_verbosity))]
    pub verbose: u64
}

fn check_file_exists(file: &OsStr) -> Result<PathBuf, OsString> {
    let path = PathBuf::from(file);
    if path.exists() {
        Ok(path)
    } else {
        Err(OsString::from(format!("{:?} does not exist", path)))
    }
}

/// Utility function to parse verbosity occurences
///
/// Up to one verbosity flags are allowed (-v), if more
/// are specified (-vv) the highest value is returned
pub fn parse_verbosity(v: u64) -> u64 {
    match v {
        0 | 1 => v,
        _ => 1,
    }
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn verbosity_exceeds_limit() {
        let passed_args = vec!["vircov", "-vv", "--paf", "tests/cases/test_ok.paf"];
        let args = Cli::from_iter_safe(passed_args);

        let actual = args.unwrap().verbose;
        let expected = 1;

        assert_eq!(actual, expected)
    }


    #[test]
    fn valid_verbosity_level() {
        let passed_args = vec!["vircov", "-v", "--paf", "tests/cases/test_ok.paf"];
        let args = Cli::from_iter_safe(passed_args);

        let actual = args.unwrap().verbose;
        let expected = 1;

        assert_eq!(actual, expected)
    }


    #[test]
    fn verbosity_from_occurrences() {
        assert_eq!(parse_verbosity(0), 0);
        assert_eq!(parse_verbosity(1), 1);
        assert_eq!(parse_verbosity(2), 1);
        assert_eq!(parse_verbosity(666), 1);
    }

}