use clap::{ArgAction, Args, Parser, Subcommand};
use std::path::PathBuf;

use crate::alignment::{Aligner, AlignmentFormat, Preset, SelectHighest};

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
    /// Alignment the scan-remap-consensus pipeline
    Run(RunArgs),
    /// Alignment coverage from alignment
    Coverage(CoverageArgs),
    /// Subtyping of consensus assemblies against curated genotype collections
    Subtype(SubtypeArgs),
    /// Subtyping of consensus assemblies against curated genotype collections
    Abundance(SubtypeArgs),
    /// Process NCBI Virus meta data files to attempt genotype extraction
    FilterDatabase(FilterDatabaseArgs),
    /// Validate genotype table order with matching sequence names at the same index
    ValidateGenotypes(ValidateGenotypesArgs),
    /// Process NCBI Virus meta data files to attempt genotype extraction
    ProcessNcbi(ProcessNcbiArgs),
    /// Process GISAID and Nextstrain files to attempt genotype extraction
    ProcessGisaid(ProcessGisaidArgs),
}

#[derive(Debug, Args)]
pub struct RunArgs {
    /// Input read files (can be compressed with .gz)
    ///
    /// One or two input read files. These files can be in gzipped format.
    /// This parameter is required and multiple files can be specified (1 for long
    /// reads or 2 for paired-end short reads) either consecutively or using multiple
    /// input arguments, for example: '-i R1.fq.gz -i R2.fq.gz' or '-i R1.fq.gz R2.fq.gz'
    #[arg(short, long, num_args(0..))]
    pub input: Vec<PathBuf>,
    /// Reference index for aligner
    ///
    /// Depending on whether --aligner is chosen, the index is an alignment index 
    /// for 'bowtie2' (index), 'minimap2' and 'strobealign' (index or fasta).
    #[arg(long, short='I', default_value="reference.fasta")]
    pub index: PathBuf,
    /// Reference sequences used in aligner ()
    #[arg(long, short='R', default_value="reference.fasta")]
    pub reference: Option<PathBuf>,
    /// Aligner
    #[arg(long, short='a', default_value="minimap2")]
    pub aligner: Aligner,
    /// Aligner preset (minimap2)
    #[arg(long, short='P', default_value="sr")]
    pub preset: Option<Preset>,
    /// Output summary table of coverage and assembly metrics
    #[arg(short, long, default_value="results.tsv")]
    pub output: PathBuf,
    /// Group alignments by a field in the reference sequence description
    #[clap(long, short = 'g', long, default_value="taxid=")]
    pub group_by: Option<String>,
    /// Group field separator in the reference sequence description
    #[clap(long, short = 's', long, default_value = ";")]
    pub group_sep: String,
    /// Select a representative genome from the groups by reads or coverage
    #[clap(long, short='b', default_value="coverage")]
    pub select_by: SelectHighest,
    /// Parallel tasks for remapping alignment
    #[clap(long, short = 'p', default_value = "4")]
    pub parallel: usize,
    /// Working directory
    #[clap(long, short = 'w', default_value = ".")]
    pub workdir: Option<PathBuf>,
    /// Keep directory with working data
    #[clap(long, short = 'k')]
    pub keep: bool,
    // Segment field to identify segments in grouped alignments
    #[clap(long, default_value="segment=")]
    pub segment_field: Option<String>,
    /// Segment field identifier negative (e.g. "segment=N/A")
    #[clap(long, short='j', default_value="segment=N/A")]
    pub segment_field_nan: Option<String>,
    /// Threads for scanning alignment
    #[clap(long, default_value = "8")]
    pub scan_threads: usize,
    /// Threads for remapping alignment
    #[clap(long, default_value = "2")]
    pub remap_threads: usize,
    /// Create consensus genome from remapping 
    #[clap(long)]
    pub no_consensus: bool,
}




#[derive(Debug, Args)]
pub struct CoverageArgs {
}

#[derive(Debug, Args)]
pub struct SubtypeArgs {
    /// Consensus assemblies for subtyping
    #[clap(long, short = 'i', num_args(0..))]
    pub input: Vec<PathBuf>,
    /// Consensus assemblies for subtyping
    #[clap(long, short = 'o')]
    pub output: PathBuf,
    /// Databases for subtyping
    #[clap(long, short = 'd')]
    pub database: PathBuf,
    /// Threads for Diamond and BLAST
    #[clap(long, short = 't', default_value = "2")]
    pub threads: u32,
    /// Minimum percent nucleotide target coverage
    #[clap(long, short = 'c', default_value="20")]
    pub min_cov: f64,
    /// Minimum percent amino acid target coverage
    #[clap(long, short = 'a', default_value="0")]
    pub min_cov_aa: f64,
    /// Minimum percent shared protein coverage
    #[clap(long, short = 'p', default_value="0")]
    pub min_cov_prot: f64,
    /// Genomic neighbor typing metric to use for inference of genotypes
    #[clap(long, short = 'm', default_value="ani")]
    pub metric: String,
    /// Show the highest number of ranked matches of the selected metric
    #[clap(long, short = 'r', default_value="5")]
    pub ranks: usize,
    /// Show the highest number of ranked matches of the selected metric for reference isolates with genotype annotation
    #[clap(long, short = 'g',)]
    pub with_genotype: bool,
    /// Show the highest number of ranked matches of the selected metric
    #[clap(long, short = 'w')]
    pub workdir: Option<PathBuf>,
    /// Keep directory with working data
    #[clap(long, short = 'k')]
    pub keep: bool,
}


#[derive(Debug, Args)]
pub struct ProcessGisaidArgs {
    /// GISAID sequence file (.fasta)
    #[clap(long, short = 'f')]
    pub fasta: PathBuf,
    /// Nextstrain clade file (.tsv)
    #[clap(long, short = 'c')]
    pub clades: PathBuf,
    /// Optional output segment annotation
    #[clap(long, short = 's')]
    pub segment: Option<String>,
    /// Vircov database sequence (db_nuc.fasta) and meta-data (db.csv) output directory
    #[clap(long, short = 'o', default_value=".")]
    pub outdir: PathBuf,
}


#[derive(Debug, Args)]
pub struct ProcessNcbiArgs {
    /// NCBI Virus meta data file 
    #[clap(long, short = 'i')]
    pub input: PathBuf,
    /// Supported virus subtype processing scheme
    #[clap(long, short = 'v')]
    pub virus: String,
    /// Vircov database meta data file 
    #[clap(long, short = 'o')]
    pub output: PathBuf,
}


#[derive(Debug, Args)]
pub struct FilterDatabaseArgs {
    /// Vircov database sequence file (.fasta) 
    #[clap(long, short = 'f')]
    pub fasta: PathBuf,
    /// Vircov database genotypes file (.csv)
    #[clap(long, short = 'g')]
    pub genotypes: PathBuf,
    /// Filtered database sequence file (nucleotide or protein) 
    #[clap(long, short = 'o')]
    pub output_fasta: Option<PathBuf>,
    /// Filtered database meta data file 
    #[clap(long, short = 't')]
    pub output_genotypes: Option<PathBuf>,
    /// Minimum sequence length in database sequence file
    #[clap(long, short = 'm', default_value="0")]
    pub min_length: usize,
    /// Remove any duplicate accessions entirely from the database
    #[clap(long, short = 'r')]
    pub remove_duplicates: bool,
    /// Filtered database meta data file 
    #[clap(long, short = 'a', num_args(0..))]
    pub accessions: Option<Vec<String>>,
    /// Filtered database meta data file 
    #[clap(long, short = 'f')]
    pub accession_file: Option<PathBuf>,
}


#[derive(Debug, Args)]
pub struct ValidateGenotypesArgs {
    /// Vircov database sequence file (nucleotide or protein) 
    #[clap(long, short = 'f')]
    pub fasta: PathBuf,
    /// Vircov database genotypes file (.csv)
    #[clap(long, short = 'g')]
    pub genotypes: PathBuf,
}

#[derive(Debug, Args)]
pub struct AniArgs {
    /// Consensus assemblies for subtyping
    #[clap(long, short = 'f', num_args(0..))]
    pub fasta: Vec<PathBuf>,
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


// #[derive(Debug, Args)]
// pub struct CoverageArgs {
//     /// Alignment file (SAM/BAM/CRAM/PAF)
//     #[clap(
//         short='i', long, required = true
//     )]
//     pub alignment: PathBuf,
//     /// bam: SAM/BAM/CRAM alignment; paf: PAF alignment
//     ///
//     /// Default is to attempt to infer the input alignment format automatically from the filename
//     /// extension (bam|sam|cram|paf). This option is used to override that.
//     #[clap(
//         long,
//         ignore_case=true
//     )]
//     pub alignment_format: Option<AlignmentFormat>,
//     /// Reference sequences in FASTA format [default: None]
//     ///
//     /// If the input format is PAF, computation of the total coverage
//     /// against each target sequence requires sequence lengths. Must
//     /// be the same reference file as used for the alignment.
//     #[clap(
//         short = 'f',
//         long = "fasta"
//     )]
//     pub fasta: Option<PathBuf>,
//     /// Minimum length of the aligned query sequence
//     ///
//     /// Filters (&) alignments by minimum length of the aligned query sequence,
//     /// which corresponds to the difference between query alignment end and
//     /// start positions.
//     #[clap(short = 'L', long, default_value = "0")]
//     pub min_len: u64,
//     /// Minimum coverage of the aligned query sequence
//     ///
//     /// Filters (&) alignments by minimum proportion of the query sequence involved
//     /// in the alignment which corresponds to division of the length of the
//     /// aligned query sequence by the length of the query sequence.
//     #[clap(short = 'C', long, default_value = "0")]
//     pub min_cov: f64,
//     /// Minimum mapping quality of the alignment
//     ///
//     /// Filters (&) alignments by a minimum mapping quality.
//     #[clap(short = 'M', long, default_value = "0")]
//     pub min_mapq: u8,
//     /// Group alignments by a field in the reference sequence description
//     ///
//     /// Grouping can help to summarize regions and other output statistics
//     /// aross members of the alignment group. Value specified as group should be
//     /// in the target sequence headers separated by a delimiter, for example
//     /// `taxid=1101 | species=BadVirus` where one can group by either the
//     /// `taxid=` or `species=` fields. Alignment target sequence headers
//     /// are split by the delimiter `|` and the field value is extracted
//     /// with whitespace trimming for each alignment.
//     /// Groups are summarized as follows:
//     ///   - distinct regions are summed
//     ///   - alignments are summed
//     ///   - unique reads are recomputed across members of the group
//     ///   - covered base pairs in the reference sequence lengths are set to 0
//     ///   - coverage is selected to be the highest by a member of the group.
//     #[clap(short = 'g', long, default_value="taxid=")]
//     pub group_by: Option<String>,
//     /// Group field separator in the reference sequence description
//     ///
//     /// Often there are multiple fields in the header of the reference
//     /// sequences (e.g. taxid=1101 | taxname=argh). A delimitor needs to
//     /// be specified to signal the end of a specific grouping field
//     /// (e.g. --group-sep "|" with --group-by "taxid=").
//     #[clap(short = 's', long, default_value = ";")]
//     pub group_sep: String,
//     /// A file with a string per line to exclude alignments
//     /// if the string occurs in the target sequence description
//     ///
//     /// This option can be used as a blacklist to filter out
//     /// alignments of unwanted viruses, e.g. using taxonomy
//     /// identifiers or species names in target sequence headers.
//     /// Lines starting with "#" and empty lines are not considered.
//     #[clap(short, long)]
//     pub exclude: Option<PathBuf>,
//     /// Prints pretty output table  
//     ///
//     /// Output the coverage statistics as a pretty table; may still get
//     /// mangled on short terminals, as terminal width cannot currently
//     /// be estimated in underlying library.
//     #[clap(short = 't', long = "table")]
//     pub table: bool,
//     /// Outputs a header in the non-table output  
//     ///
//     /// Adds a machine-readable header to the standard output
//     #[clap(short = 'H', long = "header")]
//     pub header: bool,
//     /// Prints coverage plots
//     ///
//     /// Output coverage plots below the coverage statistics.
//     #[clap(short = 'k', long = "cov-plot")]
//     pub covplot: bool,
//     /// Width of coverage plots
//     ///
//     /// Adjusts the (approximate) width of the coverage plots by
//     /// computing the bases covered by each coverage segment.
//     #[clap(short = 'w', long = "width", default_value = "100")]
//     pub width: u64,
//     /// Minimum reference sequence length
//     ///
//     /// Filters results by minimum reference sequence length
//     /// which can help remove alignments against small genes
//     /// or genome fragments.
//     #[clap(short = 'l', long = "length", default_value = "0")]
//     pub seq_len: u64,
//     /// Minimum number of coverage regions
//     ///
//     /// Filters results by a minimum number of coverage regions, the
//     /// primary output to determine a positive hit.
//     #[clap(short = 'r', long = "regions", default_value = "0")]
//     pub regions: u64,
//     /// Conditional coverage threshold for regions or grouped-regions to apply
//     ///
//     /// Applies the regions or group-regions filter only if the sequence alignment
//     /// has <= regions-coverage threshold reference coverage. Setting this value to
//     /// e.g. 0.6  only applies the regions filter to alignments with at most 60%  
//     /// coverage against the reference sequence, thus allowing for high coverage hits
//     /// with fewer distinct alignment regions to pass the basic filter. For example,
//     /// a reference sequence with 100% coverage and 1 distinct alignment regions would
//     /// not be filtered from the report output.
//     #[clap(short = 'T', long = "regions-coverage")]
//     pub regions_coverage: Option<f64>,
//     /// Minimum read threshold (unique reads in alignment)
//     ///
//     /// Filters results by a minimum reads in alignment; if results
//     /// are grouped this is done before the grouping stage to weed
//     /// out spurious alignments
//     #[clap(short = 'u', long = "reads", default_value = "0")]
//     pub reads: u64,
//     /// Minimum coverage threshold (fraction)
//     ///
//     /// Filters results by a minimum coverage across the reference sequence
//     #[clap(short = 'c', long = "coverage", default_value = "0")]
//     pub coverage: f64,
//     /// Minimum aligned reads
//     ///
//     /// Filters results by a minimum aligned reads across the reference sequence
//     #[clap(short = 'a', long = "aligned", default_value = "0")]
//     pub aligned: u64,
//     /// Minimum number of grouped coverage regions
//     ///
//     /// Filters results by a minimum number of grouped coverage regions, can be
//     /// used in addition to the pre-grouping coverage region filter
//     #[clap(long, default_value = "0")]
//     pub group_regions: u64,
//     /// Minimum grouped coverage threshold (fraction)
//     ///
//     /// Filters results by a minimum average coverage
//     /// across reference sequences if results are grouped
//     #[clap(long, default_value = "0")]
//     pub group_coverage: f64,
//     /// Minimum grouped coverage alignments
//     ///
//     /// Filters results by a minimum number of alignments
//     /// across reference sequences if results are grouped
//     #[clap(long, default_value = "0")]
//     pub group_aligned: u64,
//     /// Minimum grouped coverage reads
//     ///
//     /// Filters results by a minimum number of unique reads
//     /// across reference sequences if results are grouped
//     #[clap(long, default_value = "0")]
//     pub group_reads: u64,
//     /// Output read identifiers of all alignments to file
//     ///
//     /// Creates a file (.txt) that contains the identifiers
//     /// of all reads that passed the alignment and coverage
//     /// filters (those involved in the outputs)
//     #[clap(short = 'O', long = "read-ids")]
//     pub read_ids: Option<PathBuf>,
//     /// Output read identifiers per genome or grouped genomes to file
//     ///
//     /// Creates a directory with a file (.txt) for each genome or group
//     /// that contains the identifiers of all reads that passed the
//     /// alignment and coverage filters (those involved in the outputs)
//     #[clap(short = 'S', long = "read-ids-split")]
//     pub read_ids_split: Option<PathBuf>,
//     /// Output a single reference for a grouped coverage assessment
//     ///
//     /// Creates a directory with a file (.fasta) for a selected genome
//     /// from a grouped (e.g. species / taxid) coverage assessment using
//     /// a selection criterion as outlined in --group-select-by  
//     #[clap(short = 'R', long = "group-select-split")]
//     pub group_select_split: Option<PathBuf>,
//     /// Select a representative genome from the groups by reads or coverage
//     #[clap(
//         short = 'B',
//         long = "group-select-by",
//         value_name = "reads|coverage", // reads / coverage
//         ignore_case=true,
//         hide_possible_values=false,
//         default_value="coverage"
//     )]
//     pub group_select_by: Option<String>,
//     /// Output selected sequences with a numeric prefix sorted by descending reads or coverage (--group-select-by)
//     #[clap(short = 'G', long = "group-select-order")]
//     pub group_select_order: bool,
//     /// Output selected group summary data (coverage alignment data) to this file.
//     ///
//     /// Ungrouped summary is output to stdout. This is specifically to allow for combining
//     /// scanning and remapping steps in the Cerebro pipeline.
//     #[clap(short = 'D', long = "group-select-data")]
//     pub group_select_data: Option<PathBuf>,
//     /// Segment field identifier (e.g. "segment=")
//     ///
//     /// Use this value to identify segment fields in the referennce headers of grouped
//     /// alignments to identify whether a group contains segmented genomes. If so
//     /// the --group-select-split option is made "segment-aware". This means that
//     /// the a representative segment is selected from all segments with the same
//     /// identifier value (e.g. segment=L) using the metric in --group-select-by
//     /// and all representative segments are output into a multi-fasta files into
//     /// the directory specified by --group-select-split
//     #[clap(long)]
//     pub segment_field: Option<String>,
//     /// Segment field identifier negative (e.g. "segment=N/A")
//     ///
//     /// This value identifies reference sequences that are not segmented. When --group-select-by
//     /// is activated and --segment-field is specified, selection of segments DOES NOT occur for
//     /// grroups in which all alignments are negative.
//     #[clap(long)]
//     pub segment_field_nan: Option<String>,
//     /// Output zero-statistics reference alignments
//     ///
//     /// This option can be used to include reference sequences with and without alignments. Only
//     /// for ungrouped alignments. Requires --fasta or will throw error. If filters are activated,
//     /// this option will include all filtered sequences except statistics are now SET TO ZERO.
//     #[clap(short = 'z', long = "zero")]
//     pub zero: bool,
// }

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
