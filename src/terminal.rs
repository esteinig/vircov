use clap::builder::TypedValueParser;
use clap::{ArgAction, Args, Parser, Subcommand};
use std::path::PathBuf;

/// Vircov: metagenomic diagnostics for viral genomes
#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(styles=get_styles())]
#[command(arg_required_else_help(true))]
#[clap(name = "vircov", version)]
pub struct App {
    #[clap(subcommand)]
    pub command: Commands,
}

#[derive(Debug, Subcommand)]
pub enum Commands {
    /// Coverage statistics and selections from database alignments
    Coverage(CoverageArgs),
    /// Subtyping of consensus assemblies against curated genotype collections
    Subtype(SubtypeArgs),
}

#[derive(Debug, Args)]
pub struct SubtypeArgs {
    /// Consensus assemblies for subtyping
    #[clap(long, short = 'i', num_args(0..))]
    pub input: Vec<PathBuf>,
    /// Databases for subtyping
    #[clap(long, short = 'd')]
    pub database: PathBuf,
    /// Output directory
    #[clap(long, short = 'd', default_value = ".")]
    pub outdir: PathBuf,
    /// Threads for Diamond and BLAST
    #[clap(long, short = 't', default_value = "2")]
    pub threads: u32,
}

#[derive(Debug, Args)]
pub struct CoverageArgs {
    /// Alignment file (SAM/BAM/CRAM/PAF)
    #[clap(
        short='i', long, value_parser = clap::builder::PathBufValueParser::new().map(|p| validate_file_path(&p)), required = true
    )]
    pub alignment: PathBuf,
    /// bam: SAM/BAM/CRAM alignment; paf: PAF alignment
    ///
    /// Default is to attempt to infer the input alignment format automatically from the filename
    /// extension (bam|sam|cram|paf). This option is used to override that.
    #[clap(
        long,
        value_name = "bam|paf",
        value_parser=clap::builder::PossibleValuesParser::new(["bam", "paf"]),
        ignore_case=true,
        hide_possible_values=true
    )]
    pub alignment_format: Option<String>,
    /// Reference sequences in FASTA format [default: None]
    ///
    /// If the input format is PAF, computation of the total coverage
    /// against each target sequence requires sequence lengths. Must
    /// be the same reference file as used for the alignment.
    #[clap(
        short = 'f',
        long = "fasta", 
        value_parser = clap::builder::PathBufValueParser::new().map(|p| validate_file_path(&p))
    )]
    pub fasta: Option<PathBuf>,
    /// Minimum length of the aligned query sequence
    ///
    /// Filters (&) alignments by minimum length of the aligned query sequence,
    /// which corresponds to the difference between query alignment end and
    /// start positions.
    #[clap(short = 'L', long, default_value = "0")]
    pub min_len: u64,
    /// Minimum coverage of the aligned query sequence
    ///
    /// Filters (&) alignments by minimum proportion of the query sequence involved
    /// in the alignment which corresponds to division of the length of the
    /// aligned query sequence by the length of the query sequence.
    #[clap(short = 'C', long, default_value = "0")]
    pub min_cov: f64,
    /// Minimum mapping quality of the alignment
    ///
    /// Filters (&) alignments by a minimum mapping quality.
    #[clap(short = 'M', long, default_value = "0")]
    pub min_mapq: u8,
    /// Verbose output statistics  
    ///
    /// Single flag (-v) adds whitespace separated tags in the last column,
    /// corresponding to the number of inferred alignment coverage regions.
    /// Tag fields are separated by vertical bars (|).
    ///
    /// When no grouping argument is given then each tag consists of:
    ///     start of coverage region
    ///     | end of coverage region
    ///     | number of alignments in the region
    ///
    /// When the output is grouped, then each tag consists of data from each
    /// output that has been grouped:
    ///     name of reference sequence
    ///     | number of regions
    ///     | number of reads in region
    ///     | number of alignments in region
    ///     | number of base pairs in region
    ///     | length of reference sequence
    ///     | coverage on reference sequence
    ///     | reference sequence description header
    #[clap(short, long, action = ArgAction::Count)]
    pub verbose: u64,
    /// Group alignments by a field in the reference sequence description
    ///
    /// Grouping can help to summarize regions and other output statistics
    /// aross members of the alignment group. Value specified as group should be
    /// in the target sequence headers separated by a delimiter, for example
    /// `taxid=1101 | species=BadVirus` where one can group by either the
    /// `taxid=` or `species=` fields. Alignment target sequence headers
    /// are split by the delimiter `|` and the field value is extracted
    /// with whitespace trimming for each alignment.
    /// Groups are summarized as follows:
    ///   - distinct regions are summed
    ///   - alignments are summed
    ///   - unique reads are recomputed across members of the group
    ///   - covered base pairs in the reference sequence lengths are set to 0
    ///   - coverage is selected to be the highest by a member of the group.
    #[clap(short = 'g', long)]
    pub group_by: Option<String>,
    /// Group field separator in the reference sequence description
    ///
    /// Often there are multiple fields in the header of the reference
    /// sequences (e.g. taxid=1101 | taxname=argh). A delimitor needs to
    /// be specified to signal the end of a specific grouping field
    /// (e.g. --group-sep "|" with --group-by "taxid=").
    #[clap(short = 's', long, default_value = ";")]
    pub group_sep: String,
    /// A file with a string per line to exclude alignments
    /// if the string occurs in the target sequence description
    ///
    /// This option can be used as a blacklist to filter out
    /// alignments of unwanted viruses, e.g. using taxonomy
    /// identifiers or species names in target sequence headers.
    /// Lines starting with "#" and empty lines are not considered.
    #[clap(short, long, value_parser = clap::builder::PathBufValueParser::new().map(|p| validate_file_path(&p)))]
    pub exclude: Option<PathBuf>,
    /// Prints pretty output table  
    ///
    /// Output the coverage statistics as a pretty table; may still get
    /// mangled on short terminals, as terminal width cannot currently
    /// be estimated in underlying library.
    #[clap(short = 't', long = "table")]
    pub table: bool,
    /// Outputs a header in the non-table output  
    ///
    /// Adds a machine-readable header to the standard output
    #[clap(short = 'H', long = "header")]
    pub header: bool,
    /// Prints coverage plots
    ///
    /// Output coverage plots below the coverage statistics.
    #[clap(short = 'k', long = "cov-plot")]
    pub covplot: bool,
    /// Width of coverage plots
    ///
    /// Adjusts the (approximate) width of the coverage plots by
    /// computing the bases covered by each coverage segment.
    #[clap(short = 'w', long = "width", default_value = "100")]
    pub width: u64,
    /// Minimum reference sequence length
    ///
    /// Filters results by minimum reference sequence length
    /// which can help remove alignments against small genes
    /// or genome fragments.
    #[clap(short = 'l', long = "length", default_value = "0")]
    pub seq_len: u64,
    /// Minimum number of coverage regions
    ///
    /// Filters results by a minimum number of coverage regions, the
    /// primary output to determine a positive hit.
    #[clap(short = 'r', long = "regions", default_value = "0")]
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
    #[clap(short = 't', long = "regions-coverage")]
    pub regions_coverage: Option<f64>,
    /// Minimum read threshold (unique reads in alignment)
    ///
    /// Filters results by a minimum reads in alignment; if results
    /// are grouped this is done before the grouping stage to weed
    /// out spurious alignments
    #[clap(short = 'u', long = "reads", default_value = "0")]
    pub reads: u64,
    /// Minimum coverage threshold (fraction)
    ///
    /// Filters results by a minimum coverage across the reference sequence
    #[clap(short = 'c', long = "coverage", default_value = "0")]
    pub coverage: f64,
    /// Minimum aligned reads
    ///
    /// Filters results by a minimum aligned reads across the reference sequence
    #[clap(short = 'a', long = "aligned", default_value = "0")]
    pub aligned: u64,
    /// Minimum number of grouped coverage regions
    ///
    /// Filters results by a minimum number of grouped coverage regions, can be
    /// used in addition to the pre-grouping coverage region filter
    #[clap(long, default_value = "0")]
    pub group_regions: u64,
    /// Minimum grouped coverage threshold (fraction)
    ///
    /// Filters results by a minimum average coverage
    /// across reference sequences if results are grouped
    #[clap(long, default_value = "0")]
    pub group_coverage: f64,
    /// Minimum grouped coverage alignments
    ///
    /// Filters results by a minimum number of alignments
    /// across reference sequences if results are grouped
    #[clap(long, default_value = "0")]
    pub group_aligned: u64,
    /// Minimum grouped coverage reads
    ///
    /// Filters results by a minimum number of unique reads
    /// across reference sequences if results are grouped
    #[clap(long, default_value = "0")]
    pub group_reads: u64,
    /// Output read identifiers of all alignments to file
    ///
    /// Creates a file (.txt) that contains the identifiers
    /// of all reads that passed the alignment and coverage
    /// filters (those involved in the outputs)
    #[clap(short = 'O', long = "read-ids")]
    pub read_ids: Option<PathBuf>,
    /// Output read identifiers per genome or grouped genomes to file
    ///
    /// Creates a directory with a file (.txt) for each genome or group
    /// that contains the identifiers of all reads that passed the
    /// alignment and coverage filters (those involved in the outputs)
    #[clap(short = 'S', long = "read-ids-split")]
    pub read_ids_split: Option<PathBuf>,
    /// Output a single reference for a grouped coverage assessment
    ///
    /// Creates a directory with a file (.fasta) for a selected genome
    /// from a grouped (e.g. species / taxid) coverage assessment using
    /// a selection criterion as outlined in --group-select-by  
    #[clap(short = 'R', long = "group-select-split")]
    pub group_select_split: Option<PathBuf>,
    /// Select a representative genome from the groups by reads or coverage
    #[clap(
        short = 'B',
        long = "group-select-by",
        value_name = "coverage",
        value_parser=clap::builder::PossibleValuesParser::new(["reads", "coverage"]),
        ignore_case=true,
        hide_possible_values=false
    )]
    pub group_select_by: Option<String>,
    /// Output selected sequences with a numeric prefix sorted by descending reads or coverage (--group-select-by)
    #[clap(short = 'G', long = "group-select-order")]
    pub group_select_order: bool,
    /// Output selected group summary data (coverage alignment data) to this file.
    ///
    /// Ungrouped summary is output to stdout. This is specifically to allow for combining
    /// scanning and remapping steps in the Cerebro pipeline.
    #[clap(short = 'D', long = "group-select-data")]
    pub group_select_data: Option<PathBuf>,
    /// Segment field identifier (e.g. "segment=")
    ///
    /// Use this value to identify segment fields in the referennce headers of grouped
    /// alignments to identify whether a group contains segmented genomes. If so
    /// the --group-select-split option is made "segment-aware". This means that
    /// the a representative segment is selected from all segments with the same
    /// identifier value (e.g. segment=L) using the metric in --group-select-by
    /// and all representative segments are output into a multi-fasta files into
    /// the directory specified by --group-select-split
    #[clap(long)]
    pub segment_field: Option<String>,
    /// Segment field identifier negative (e.g. "segment=N/A")
    ///
    /// This value identifies reference sequences that are not segmented. When --group-select-by
    /// is activated and --segment-field is specified, selection of segments DOES NOT occur for
    /// grroups in which all alignments are negative.
    #[clap(long)]
    pub segment_field_nan: Option<String>,
    /// Output zero-statistics reference alignments
    ///
    /// This option can be used to include reference sequences with and without alignments. Only
    /// for ungrouped alignments. Requires --fasta or will throw error. If filters are activated,
    /// this option will include all filtered sequences except statistics are now SET TO ZERO.
    #[clap(short = 'z', long = "zero")]
    pub zero: bool,
}

pub fn get_styles() -> clap::builder::Styles {
    clap::builder::Styles::styled()
        .header(
            anstyle::Style::new()
                .bold()
                .underline()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Yellow))),
        )
        .literal(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Green))),
        )
}

fn validate_file_path(path: &PathBuf) -> Result<(), String> {
    if !path.exists() {
        return Err(format!("File path does not exist: {}", path.display()));
    }
    Ok(())
}
