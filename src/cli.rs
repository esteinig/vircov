use std::ffi::{OsStr, OsString};
use std::path::PathBuf;
use structopt::StructOpt;

/// Metacov Command Line client
#[derive(Debug, StructOpt)]
#[structopt(name = "vircov")]
pub struct Cli {
    /// Alignment file (SAM/BAM/CRAM/PAF)
    #[structopt(
        short, long, parse(try_from_os_str = check_file_exists), required = true
    )]
    pub alignment: PathBuf, // TODO
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
    pub alignment_format: Option<String>,
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
    pub fasta: Option<PathBuf>,
    /// Minimum length of the aligned query sequence
    ///
    /// Filters (&) alignments by minimum length of the aligned query sequence,
    /// which corresponds to the difference between query alignment end and
    /// start positions.
    #[structopt(short = "L", long, default_value = "0")]
    pub min_len: u64,
    /// Minimum coverage of the aligned query sequence
    ///
    /// Filters (&) alignments by minimum proportion of the query sequence involved
    /// in the alignment which corresponds to division of the length of the
    /// aligned query sequence by the length of the query sequence.
    #[structopt(short = "C", long, default_value = "0")]
    pub min_cov: f64,
    /// Minimum mapping quality of the alignment
    ///
    /// Filters (&) alignments by a minimum mapping quality.
    #[structopt(short = "M", long, default_value = "0")]
    pub min_mapq: u8,
    /// Verbose output statistics  
    ///
    /// Single flag (-v) adds whitespace separated tags in the last column,
    /// corresponding to the number of inferred alignment coverage regions.
    /// Tag fields are separated by vertical bars (|).
    ///
    /// When no grouping argument is given then each tag consists of
    ///     | start of coverage region
    ///     | end of coverage region
    ///     | number of alignments in the region
    ///
    /// When the output is grouped, then each tag consists of data from each
    /// output that has been grouped
    ///     | name of reference sequence
    ///     | number of regions
    ///     | number of reads in region
    ///     | number of alignments in region
    ///     | number of base pairs in region
    ///     | length of reference sequence
    ///     | coverage on reference sequence
    ///     | reference sequence description header
    #[structopt(short, long, parse(from_occurrences = parse_verbosity))]
    pub verbose: u64,
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
    pub group_by: Option<String>,
    /// Group field separator in the reference sequence description
    ///
    /// Often there are multiple fields in the header of the reference
    /// sequences (e.g. taxid=1101 | taxname=argh). A delimitor needs to
    /// be specified to signal the end of a specific grouping field
    /// (e.g. --group-sep "|" with --group-by "taxid=").
    #[structopt(long, default_value = ";")]
    pub group_sep: String,
    /// A file with a string per line to exclude alignments
    /// if the string occurs in the target sequence description
    ///
    /// This option can be used as a blacklist to filter out
    /// alignments of unwanted viruses, e.g. using taxonomy
    /// identifiers or species names in target sequence headers.
    #[structopt(short, long, parse(try_from_os_str = check_file_exists))]
    pub exclude: Option<PathBuf>,
    /// Prints pretty output table  
    ///
    /// Output the coverage statistics as a pretty table; may still get
    /// mangled on short terminals, as terminal width cannot currently
    /// be estimated in underlying library.
    #[structopt(short = "t", long = "table")]
    pub table: bool,
    /// Prints coverage plots
    ///
    /// Output coverage plots below the coverage statistics.
    #[structopt(short = "k", long = "cov-plot")]
    pub covplot: bool,
    /// Width of coverage plots
    ///
    /// Adjusts the (approximate) width of the coverage plots by
    /// computing the bases covered by each coverage segment.
    #[structopt(short = "w", long = "width", default_value = "100")]
    pub width: u64,
    /// Minimum reference sequence length
    ///
    /// Filters results by minimum reference sequence length
    /// which can help remove alignments against small genes
    /// or genome fragments.
    #[structopt(short = "l", long = "length", default_value = "0")]
    pub seq_len: u64,
    /// Minimum number of coverage regions
    ///
    /// Filters results by a minimum number of coverage regions, the
    /// primary output to determine a positive hit.
    #[structopt(short = "r", long = "regions", default_value = "0")]
    pub regions: u64,
    /// Conditional coverage threshold for regions or grouped-regions to apply
    ///
    /// Applies the regions or group-regions filter only if the sequence alignment
    /// has <= regions-coverage threshold reference coverage. Setting this value to
    /// e.g. 0.6  only applies the regions filter to alignments with at most 60%  
    /// coverage against the reference sequence, thus allowing for high coverage hits
    /// with fewer distinct alignment regions to pass the basic filter. For example,
    /// a reference sequence with 100% coverage and 1 distinct alignment regions would
    /// not be filtered from the report output.
    #[structopt(short = "r", long = "regions-cov")]
    pub regions_cov: Option<f64>,
    /// Minimum read threshold (unique reads in alignment)
    ///
    /// Filters results by a minimum reads in alignment; if results
    /// are grouped this is done before the grouping stage to weed
    /// out spurious alignments
    #[structopt(short = "u", long = "reads", default_value = "0")]
    pub reads: u64,
    /// Minimum coverage threshold (fraction)
    ///
    /// Filters results by a minimum coverage across the reference sequence
    #[structopt(short = "c", long = "coverage", default_value = "0")]
    pub coverage: f64,
    /// Minimum number of grouped coverage regions
    ///
    /// Filters results by a minimum number of grouped coverage regions, can be
    /// used in addition to the pre-grouping coverage region filter
    #[structopt(long, default_value = "0")]
    pub group_regions: u64,
    /// Minimum grouped coverage threshold (fraction)
    ///
    /// Filters results by a minimum average coverage
    /// across reference sequences if results are grouped
    #[structopt(long, default_value = "0")]
    pub group_coverage: f64,
    /// Output read identifiers of all alignments to file
    ///
    /// Creates a file (.txt) that contains the identifiers
    /// of all reads that passed the alignment and coverage
    /// filters (those involved in the outputs)
    #[structopt(short = "O", long = "read-ids")]
    pub read_ids: Option<PathBuf>,
    /// Output read identifiers per genome or grouped genomes to file
    ///
    /// Creates a directory with a file (.txt) for each genome or group
    /// that contains the identifiers of all reads that passed the
    /// alignment and coverage filters (those involved in the outputs)
    #[structopt(short = "S", long = "read-ids-split")]
    pub read_ids_split: Option<PathBuf>,
    /// Output a single reference for a grouped coverage assessment
    ///
    /// Creates a directory with a file (.fasta) for a selected genome
    /// from a grouped (e.g. species / taxid) coverage assessment using
    /// a selection criterion as outlined in --group-select-by  
    #[structopt(short = "R", long = "group-select-split")]
    pub group_select_split: Option<PathBuf>,
    /// Select a representative genome from the groups by max reads or coverage
    #[structopt(
        short = "B", 
        long = "group-select-by",
        value_name = "coverage",
        possible_values = &["reads", "coverage"],
        case_insensitive=true,
        hide_possible_values=true
    )]
    pub group_select_by: Option<String>,
    /// Output selected sequences with a numeric prefix sorted by descending coverage
    #[structopt(short = "G", long = "group-select-order")]
    pub group_select_order: Option<bool>,
    /// Segment field identifier (e.g. "segment=")
    ///
    /// Use this value to identify segment fields in the referennce headers of grouped
    /// alignments to identify whether a group contains segmented genomes. If so
    /// the --group-select-split option is made "segment-aware". This means that
    /// the a representative segment is selected from all segments with the same
    /// identifier value (e.g. segment=L) using the metric in --group-select-by
    /// and all representative segments are output into a multi-fasta files into
    /// the directory specified by --group-select-split
    #[structopt(long)]
    pub segment_field: Option<String>,
    /// Segment field identifier negative (e.g. "segment=N/A")
    ///
    /// This value identifies reference sequences that are not segmented. When --group-select-by
    /// is activated and --segment-field is specified, selection of segments DOES NOT occur for
    /// grroups in which all alignments are negative.
    #[structopt(long)]
    pub segment_field_nan: Option<String>,
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
}
