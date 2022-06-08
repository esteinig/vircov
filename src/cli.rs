use std::ffi::{OsStr, OsString};
use std::path::PathBuf;
use structopt::StructOpt;

/// Metacov Command Line client
#[derive(Debug, StructOpt)]
#[structopt(name = "lab-tools")]
pub struct Cli {
    #[structopt(subcommand)]
    pub commands: Commands,
}

#[derive(Debug, StructOpt)]
pub enum Commands {
    Vircov {
        /// Alignment file (SAM/BAM/CRAM/PAF)
        #[structopt(
            short, long, parse(try_from_os_str = check_file_exists), required = true
        )]
        alignment: PathBuf,
        /// bam: SAM/BAM/CRAM alignment; paf: PAF alignment
        ///
        /// Default is to attempt to infer the input alignment format automatically from the filename
        /// extension (bam|sam|cram|paf). This option is used to override that.
        #[structopt(
            short = "A",
            long,
            value_name = "bam|paf",
            possible_values = &["bam", "paf"],
            case_insensitive=true,
            hide_possible_values=true
        )]
        alignment_format: Option<String>,
        /// Reference sequences in FASTA format [default: None]
        ///
        /// If the input format is PAF, computation of the total coverage
        /// against each target sequence requires sequence lengths. Must
        /// be the same reference file as used for the alignment.
        #[structopt(
            short = "f",
            long = "fasta", 
            parse(try_from_os_str = check_file_exists)
        )]
        fasta: Option<PathBuf>,
        /// Minimum length of the aligned query sequence
        ///
        /// Filters (&) alignments by minimum length of the aligned query sequence,
        /// which corresponds to the difference between query alignment end and
        /// start positions.
        #[structopt(short = "L", long, default_value = "0")]
        min_len: u64,
        /// Minimum coverage of the aligned query sequence
        ///
        /// Filters (&) alignments by minimum proportion of the query sequence involved
        /// in the alignment which corresponds to division of the length of the
        /// aligned query sequence by the length of the query sequence.
        #[structopt(short = "C", long, default_value = "0")]
        min_cov: f64,
        /// Minimum mapping quality of the alignment
        ///
        /// Filters (&) alignments by a minimum mapping quality.
        #[structopt(short = "M", long, default_value = "0")]
        min_mapq: u8,
        /// Verbose output statistics  
        ///
        /// Single flag (-v) adds whitespace separated tags in the last column,
        /// corresponding to the number of inferred alignment coverage regions.
        /// Tag fields are separated by semicolons.
        ///
        /// When no grouping argument is given then each tag consists of
        ///     : start of coverage region
        ///     : end of coverage region
        ///     : number of alignments in the region
        ///
        /// When the output is grouped, then each tag consists of data from each
        /// output that has been grouped
        ///     : name of reference sequence
        ///     : number of regions
        ///     : number of reads in region
        ///     : number of alignments in region
        ///     : number of base pairs in region
        ///     : length of reference sequence
        ///     : coverage on reference sequence
        #[structopt(short, long, parse(from_occurrences = parse_verbosity))]
        verbose: u64,
        /// Group outputs by a field in the reference sequence description
        ///
        /// Grouping can help to summarize regions and other output statistics
        /// aross all members of the group. Value specified as group should be
        /// in the target sequence headers separated by a delimiter, for example
        /// `taxid=1101 | species=MeanVirus` where one can group by either the
        /// `taxid=` or `species=` fields (including comma). Alignment target
        /// sequence are split by the delimiter `|` and the field value is
        /// extracted for each alignment. Groups are summarized as follows:
        /// distinct regions are summed, alignments are summed, unique reads are
        /// recomputed across all members of the group, covered base pairs in the
        /// alignments and reference sequence lengths are set to 0, and coverage
        /// across all target sequences in the group is averaged.
        #[structopt(short, long)]
        group_by: Option<String>,
        /// Group field separator in the reference sequence description
        ///
        /// Often there are multiple fields in the header of the reference
        /// sequences (e.g. taxid=1101 | taxname=argh). A delimitor needs to
        /// be specified to signal the end of a specific grouping field
        /// (e.g. --group-sep "|" with --group-by "taxid=").
        #[structopt(long, default_value = ";")]
        group_sep: String,
        /// A file with a string per line to exclude alignments
        /// if the string occurs in the target sequence description
        ///
        /// This option can be used as a blacklist to filter out
        /// alignments of unwanted viruses, e.g. using taxonomy
        /// identifiers or species names in target sequence headers.
        #[structopt(short, long, parse(try_from_os_str = check_file_exists))]
        exclude: Option<PathBuf>,
        /// Prints pretty output table  
        ///
        /// Output the coverage statistics as a pretty table; may still get
        /// mangled on short terminals, as terminal width cannot currently
        /// be estimated in underlying library.
        #[structopt(short = "t", long = "table")]
        table: bool,
        /// Prints coverage plots
        ///
        /// Output coverage plots below the coverage statistics.
        #[structopt(short = "k", long = "cov-plot")]
        covplot: bool,
        /// Width of coverage plots
        ///
        /// Adjusts the (approximate) width of the coverage plots by
        /// computing the bases covered by each coverage segment.
        #[structopt(short = "w", long = "width", default_value = "100")]
        width: u64,
        /// Minimum reference sequence length
        ///
        /// Filters results by minimum reference sequence length
        /// which can help remove alignments against small genes
        /// or genome fragments.
        #[structopt(short = "l", long = "length", default_value = "0")]
        seq_len: u64,
        /// Minimum number of coverage regions
        ///
        /// Filters results by a minimum number of coverage regions, the
        /// primary output to determine a positive hit.
        #[structopt(short = "r", long = "regions", default_value = "0")]
        regions: u64,
        /// Minimum read threshold (unique reads in alignment)
        ///
        /// Filters results by a minimum reads in alignment; if results
        /// are grouped this is done before the grouping stage to weed
        /// out spurious alignments
        #[structopt(short = "u", long = "reads", default_value = "0")]
        reads: u64,
        /// Minimum coverage threshold (fraction)
        ///
        /// Filters results by a minimum coverage across the reference sequence
        #[structopt(short = "c", long = "coverage", default_value = "0")]
        coverage: f64,
        /// Minimum number of grouped coverage regions
        ///
        /// Filters results by a minimum number of grouped coverage regions, can be
        /// used in addition to the pre-grouping coverage region filter
        #[structopt(long, default_value = "0")]
        group_regions: u64,
        /// Minimum grouped coverage threshold (fraction)
        ///
        /// Filters results by a minimum average coverage
        /// across reference sequences if results are grouped
        #[structopt(long, default_value = "0")]
        group_coverage: f64,
        /// Output read identifiers of aligned reads to file
        ///
        /// Creates a file that contains the (unique) identifiers
        /// of all reads that passed the alignment and coverage
        /// filters (i.e. those involved in the table outputs)
        #[structopt(short = "O", long = "read-ids")]
        read_ids: Option<PathBuf>,
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

/// Utility function to parse verbosity occurences
///
/// Up to one verbosity flags are allowed (-v), if more
/// are specified (-vv) the highest value is returned
fn parse_verbosity(v: u64) -> u64 {
    match v {
        0 | 1 | 2 => v,
        _ => 2,
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
        let passed_args = vec!["metacov", "-vvv", "--paf", "tests/cases/test_ok.paf"];
        let args = Cli::from_iter_safe(passed_args);
        let actual = args.unwrap().verbose;
        let expected = 2;
        assert_eq!(actual, expected)
    }

    #[test]
    fn valid_verbosity_level() {
        let passed_args = vec!["metacov", "-v", "--paf", "tests/cases/test_ok.paf"];
        let args = Cli::from_iter_safe(passed_args);
        let actual = args.unwrap().verbose;
        let expected = 1;
        assert_eq!(actual, expected)
    }

    #[test]
    fn verbosity_from_occurrences() {
        assert_eq!(parse_verbosity(0), 0);
        assert_eq!(parse_verbosity(1), 1);
        assert_eq!(parse_verbosity(2), 2);
        assert_eq!(parse_verbosity(666), 2);
    }
}
