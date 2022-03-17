use std::ffi::{OsStr, OsString};
use std::path::PathBuf;
use structopt::StructOpt;

// Vircov command-line client
#[derive(Debug, StructOpt)]
#[structopt(name = "vircov")]
pub struct Cli {
    /// Alignment file (PAF)
    ///
    /// Alignment input file, currently supports PAF (minimap2)
    #[structopt(
        parse(try_from_os_str = check_file_exists)
    )]
    pub path: PathBuf,
    /// Reference sequences in FASTA format
    ///
    /// If the input are alignments in PAF format, computation of the
    /// total coverage against each target sequence requires the sequence
    /// lengths extracted from the FASTA file used in the alignment.
    #[structopt(
        short = "f",
        long = "fasta", 
        parse(try_from_os_str = check_file_exists)
    )]
    pub fasta: Option<PathBuf>,
    /// Minimum length of the aligned query sequence
    ///
    /// Filters (&) alignments by minimum length of the aligned query sequence,
    /// which corresponds to the difference between query alignment end and
    /// start positions.
    #[structopt(short = "l", long = "min-len", default_value = "50")]
    pub min_len: u64,
    /// Minimum coverage of the aligned query sequence
    ///
    /// Filters (&) alignmentsby minimum proportion of the query sequence involved
    /// in the alignment which corresponds to the division of the length of the
    /// aligned query sequence by the length of the query sequence.
    #[structopt(short = "c", long = "min-cov", default_value = "0")]
    pub min_cov: f64,
    /// Minimum mapping quality of the alignment
    ///
    /// Filters (&) alignments by a minimum mapping quality.
    #[structopt(short = "q", long = "min-mapq", default_value = "30")]
    pub min_mapq: u8,
    /// Verbose output statistics  
    ///
    /// Single flag (-v) adds whitespace separated tags in the last column,
    /// corresponding to the number of inferred alignment coverage blocks.
    /// Each tag contains: start of block, end of block and number of
    /// alignments in the block, separated by semicolons (e.g.3450:3500:9)
    #[structopt(short, long, parse(from_occurrences = parse_verbosity))]
    pub verbose: u64,
    /// Pretty print output table  
    ///
    /// Output the coverage statistics as a pretty table.
    #[structopt(short, long)]
    pub table: bool,
    /// Prints coverage plots
    ///
    /// Output coverage plots below the coverage statistics.
    #[structopt(short = "p", long)]
    pub covplot: bool,
    /// Width of the coverage plots
    ///
    /// Adjusts the (approximate) width of the coverage plots by
    /// computing the bases covered by each coverage segment.
    #[structopt(short = "w", long, default_value = "100")]
    pub width: u64,
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
    #[should_panic]
    fn file_check_not_exist() {
        check_file_exists(OsStr::new("tests/cases/no_bueno.paf")).unwrap();
    }

    #[test]
    fn verbosity_exceeds_limit() {
        let passed_args = vec!["vircov", "-vv", "tests/cases/test_ok.paf"];
        let args = Cli::from_iter_safe(passed_args);
        let actual = args.unwrap().verbose;
        let expected = 1;
        assert_eq!(actual, expected)
    }

    #[test]
    fn valid_verbosity_level() {
        let passed_args = vec!["vircov", "-v", "tests/cases/test_ok.paf"];
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
