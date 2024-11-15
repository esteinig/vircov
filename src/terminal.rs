use clap::{Args, Parser, Subcommand};
use std::path::PathBuf;

use crate::{alignment::{Aligner, Preset, SelectHighest}, annotation::AnnotationPreset, haplotype::Haplotyper};

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
    /// Alignment and consensus assembly pipeline
    Run(RunArgs),
    /// Coverage statistics from alignments
    Coverage(CoverageArgs),
    /// Abundance estimates from alignments
    Abundance(SubtypeArgs),
    /// Assembly subtyping using nearest neighbor graphs
    Subtype(SubtypeArgs),

    #[clap(subcommand)]
    /// Utilities for processing databases and outputs
    Tools(ToolsCommands)
}

#[derive(Debug, Args)]
pub struct RunArgs {
    /// Input read files (can be compressed with .gz)
    ///
    /// One or two input read files. These files can be in gzipped format.
    /// This parameter is required and multiple files can be specified (1 for long
    /// reads or 2 for paired-end short reads) either consecutively or using multiple
    /// input arguments, for example: '-i R1.fq.gz -i R2.fq.gz' or '-i R1.fq.gz R2.fq.gz'
    #[arg(short, long, num_args(1..), required=true, value_parser = validate_file)]
    pub input: Vec<PathBuf>,
    /// Output summary table of coverage and assembly metrics (.tsv)
    #[arg(short, long)]
    pub output: PathBuf,
    /// Reference sequences with header annotations for binning
    #[arg(long, short='r')]
    pub reference: PathBuf,
    /// Aligner used for scanning and binned realignment 
    #[arg(long, short='a', default_value="minimap2", help_heading="Alignment options")]
    pub aligner: Aligner,
    /// Reference index for aligner
    ///
    /// Alignment index for 'bowtie2' (index), 'minimap2' (index | fasta) 
    /// and 'strobealign' (index | fasta)
    #[arg(long, short='I', help_heading="Alignment options")]
    pub index: Option<PathBuf>,
    /// Aligner preset (minimap2)
    #[arg(long, short='P', default_value="sr", help_heading="Alignment options")]
    pub preset: Option<Preset>,
    /// Include secondary alignments in scanning alignment
    #[clap(long, help_heading="Alignment options")]
    pub secondary: bool,
    /// Annotation preset for reference headers (examples in long help)
    /// 
    /// 'default': bin=Species sp.|segment=NAN|name=Species|description=Species genome description
    /// 'ictv': species=Species sp.|segment=NAN|name=Species|description=Species genome description
    /// 'virosaurus': taxid=Species sp.;segment=N/A;name=Species;description=Species genome description
    #[arg(long, short='A', default_value="default", help_heading="Binning options")]
    pub annotation_preset: AnnotationPreset,
    /// Select a representative genome from the groups by reads or coverage
    #[clap(long, short='s', default_value="coverage", help_heading="Binning options")]
    pub select_by: SelectHighest,
    /// Minimum depth of bin-reference coverage for summary table (summary metric, does not affect other metrics)
    #[clap(long, default_value="1", help_heading="Filter options")]
    pub min_depth_coverage: Option<usize>,
    /// Minimum bin-reference coverage required for consensus assembly (fraction, 0.0 - 1.0) 
    #[clap(long, default_value="0.2", help_heading="Filter options")]
    pub min_remap_coverage: f64,
    /// Working directory
    #[clap(long, short = 'w', help_heading="Output options")]
    pub workdir: Option<PathBuf>,
    /// Keep directory with working data
    #[clap(long, short = 'k', help_heading="Output options")]
    pub keep: bool,
    /// Print formatted table to console
    #[clap(long, short = 'T', help_heading="Output options")]
    pub table: bool,
    /// Include all scanning records in output table
    #[clap(long, short='S', help_heading="Output options")]
    pub include_scans: bool,
    /// Threads for scanning alignment
    #[clap(long, default_value = "8", help_heading="Scanning stage")]
    pub scan_threads: usize,
    /// Additional arguments for scanning stage aligner
    #[clap(long, allow_hyphen_values=true, help_heading="Scanning stage")]
    pub scan_args: Option<String>,
    /// Additional arguments for scanning stage alignment filter (samtools)
    #[clap(long, allow_hyphen_values=true, default_value="-F 4", help_heading="Scanning stage")]
    pub scan_filter_args: Option<String>,
    /// Threads for remapping against bin reference
    #[clap(long, default_value = "2", help_heading="Remapping stage")]
    pub remap_threads: usize,
    /// Parallel tasks for remapping against bin reference
    #[clap(long, short = 'p', default_value = "4", help_heading="Remapping stage")]
    pub remap_parallel: usize,
    /// Additional arguments for remap stage aligner
    #[clap(long, allow_hyphen_values=true, help_heading="Remapping stage")]
    pub remap_args: Option<String>,
    /// Additional arguments for scanning stage alignment filter (samtools)
    #[clap(long, allow_hyphen_values=true, default_value="-F 4", help_heading="Remapping stage")]
    pub remap_filter_args: Option<String>,
    /// Remap all input reads instead of binned reads
    #[clap(long, help_heading="Remapping stage")]
    pub remap_all: bool,
    /// Exclude sequences from remapping according to their 'bin' (or equivalent annotation for scheme)
    /// 
    /// Can be used for example to exclude reads aligning to a host genome as a bait in the reference. 
    /// Use '--include-scans' to include the excluded scanning alignments in the output table, by default 
    /// these will not show due to not being taken into the remapping stage. 
    #[clap(long, num_args(1..), help_heading="Remapping stage")]
    pub remap_exclude_bins: Option<Vec<String>>,
    /// Do not create consensus genome from remapping stage
    #[clap(long, help_heading="Consensus stage")]
    pub consensus_disabled: bool,
    /// Maximum depth fo 'mpileup' to take into consensus assembly
    #[clap(long, default_value="10000", help_heading="Consensus stage")]
    pub consensus_max_depth: usize,
    /// Minimum consensus assembly read depth to call a site
    #[clap(long, default_value="3", help_heading="Consensus stage")]
    pub consensus_min_depth: usize,
    /// Minimum consensus frequency to call a variant site
    #[clap(long, default_value="0.75", help_heading="Consensus stage")]
    pub consensus_min_frequency: f64,
    /// Minimum base quality to consider a site
    #[clap(long, default_value="20", help_heading="Consensus stage")]
    pub consensus_min_quality: usize,
    /// Minimum consensus completeness to take consensus into remapping
    #[clap(long, default_value="60.0", help_heading="Consensus stage")]
    pub consensus_min_completeness: f64,
    /// Attempt haplotyping to recover strains from mixed samples
    #[clap(long, help_heading="Haplotype stage")]
    pub haplotype: bool,
    /// Haplotype calling pipelines available for strain haplotyping
    #[clap(long, default_value="floria", help_heading="Haplotype stage")]
    pub haplotyper: Haplotyper,
}

#[derive(Debug, Args)]
pub struct CoverageArgs {
    /// Input read files (can be compressed with .gz)
    ///
    /// One or two input read files. These files can be in gzipped format.
    /// This parameter is required and multiple files can be specified (1 for long
    /// reads or 2 for paired-end short reads) either consecutively or using multiple
    /// input arguments, for example: '-i R1.fq.gz -i R2.fq.gz' or '-i R1.fq.gz R2.fq.gz'
    #[arg(long, short, num_args(1..), value_parser = validate_file)]
    pub input: Vec<PathBuf>,
    /// Output summary table of coverage metrics (.tsv)
    #[arg(long, short)]
    pub output: PathBuf,
    /// Reference sequences used in alignment (required for --zero)
    #[arg(long, short='r')]
    pub reference: Option<PathBuf>,
    /// Aligner used for scanning coverage
    #[arg(long, short='a', default_value="minimap2")]
    pub aligner: Aligner,
    /// Reference index for aligner
    ///
    /// Depending on whether --aligner is chosen, the index is an alignment index 
    /// for 'bowtie2' (index), 'minimap2' and 'strobealign' (index or fasta).
    #[arg(long, short='I')]
    pub index: Option<PathBuf>,
    /// Aligner preset (minimap2)
    #[arg(long, short='P', default_value="sr")]
    pub preset: Option<Preset>,
    /// Annotation preset for reference headers
    #[arg(long, short='A', default_value="default")]
    pub annotation_preset: AnnotationPreset,
    /// Working directory
    #[clap(long, short = 'w', default_value = ".")]
    pub workdir: Option<PathBuf>,
    /// Keep directory with working data
    #[clap(long, short = 'k')]
    pub keep: bool,
    /// Print formatted table to console
    #[clap(long, short = 'T')]
    pub table: bool,
    /// Threads for alignment
    #[clap(long, short = 't', default_value = "4")]
    pub threads: usize,
    /// Additional arguments for aligner
    #[clap(long, allow_hyphen_values=true)]
    pub args: Option<String>,
    /// Output aligned read identifiers
    #[clap(long)]
    pub read_id: Option<PathBuf>,
    /// Include secondary alignments
    #[clap(long)]
    pub secondary: bool,
    /// Include interval alignments in output
    #[arg(long)]
    pub intervals: bool,
    /// Include references without alignment
    #[arg(long)]
    pub zero: bool,
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


#[derive(Debug, Subcommand)]
pub enum ToolsCommands {
    /// Process NCBI Virus meta data files to attempt genotype extraction
    AnnotateDatabase(AnnotateDatabaseArgs),
    /// Process NCBI Virus meta data files to attempt genotype extraction
    FilterDatabase(FilterDatabaseArgs),
    /// Validate genotype table order with matching sequence names at the same index
    ValidateGenotypes(ValidateGenotypesArgs),
    /// Process NCBI Virus meta data files to attempt genotype extraction
    ProcessNcbi(ProcessNcbiArgs),
    /// Process GISAID and Nextstrain files to attempt genotype extraction
    ProcessGisaid(ProcessGisaidArgs),
    /// Concatenate output tables
    ConcatOutput(ConcatArgs),
    /// Filter output tables 
    FilterOutput(FilterArgs)
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
pub struct ConcatArgs {
    /// Vircov run output table
    #[clap(long, short = 'i', num_args(0..))]
    pub input: Vec<PathBuf>,
    /// Concatenated output file
    #[clap(long, short = 'o')]
    pub output: PathBuf,
    /// Filter by minimum consensus completeness
    #[clap(long, short = 'm')]
    pub min_completeness: Option<f64>,
    /// Add the file parent directory name to column 'id'
    #[clap(long, short = 'd')]
    pub file_dir: bool,
    /// Add the file stem to column 'id'
    #[clap(long, short = 'f')]
    pub file_id: bool,
}

#[derive(Debug, Args)]
pub struct FilterArgs {
    /// Vircov run output table
    #[clap(long, short = 'i')]
    pub input: PathBuf,
    /// Filtered output file
    #[clap(long, short = 'o')]
    pub output: PathBuf,
    /// Filter by sample identifier (if concatenated)
    #[clap(long, num_args(0..))]
    pub id: Option<Vec<String>>,
    /// Filter by bin name
    #[clap(long, short = 'b')]
    pub bin: Option<String>,
    /// Filter by minimum consensus completeness percent
    #[clap(long, short = 'm')]
    pub min_consensus_completeness: Option<f64>,
    /// Filter by minimum consensus coverage percent (unique with MAPQ) 
    #[clap(long, short = 'n')]
    pub min_consensus_coverage_mapq: Option<f64>,
    /// Filter by minimum remap coverage
    #[clap(long, short = 'c')]
    pub min_remap_coverage: Option<f64>,
    /// Filter by minimum remap coverage
    #[clap(long, short = 'd')]
    pub min_remap_depth_coverage: Option<f64>,
    /// Filter by minimum scan alignments
    #[clap(long, short = 'a')]
    pub min_scan_alignments: Option<u64>,

}

#[derive(Debug, Args)]
pub struct AnnotateDatabaseArgs {
    /// Vircov database sequence file (.fasta) 
    #[clap(long, short = 'f')]
    pub fasta: PathBuf,
    /// Vircov database sequence annotation file (.tsv)
    /// 
    /// Required columns: id (string, sequence identifier, required value), bin (string, binning variable, required value), 
    /// segment (string, segment annotation, optional value), name (string, organism name, optional value), 
    /// description (string, sequence description, optional value)
    #[clap(long, short = 'a')]
    pub annotations: PathBuf,
    /// Annotated database sequence file (.fasta)
    #[clap(long, short = 'o')]
    pub output: PathBuf,
    /// Annotation configuration preset
    #[clap(long, short = 'p', default_value="default")]
    pub preset: AnnotationPreset,
    /// Skipped sequence identifiers for which provided annotations were missing (.tsv)
    #[clap(long, short = 's')]
    pub skipped: Option<PathBuf>,
    
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
pub struct DistArgs {
    /// Genomes for pairwise distance matrix in single file (.fasta)
    #[clap(long, short = 'f')]
    pub fasta: PathBuf,
    /// Output pairwise distance matrix as tab-delimited text file 
    #[clap(long, short = 'd')]
    pub dist: PathBuf,
    /// Output pairwise alignment fraction matrix as tab-delimited text file 
    #[clap(long, short = 'a')]
    pub afrac: Option<PathBuf>,
    /// Databases for subtyping
    #[clap(long, short = 'c', default_value="30")]
    pub compression_factor: usize,
    /// Output directory
    #[clap(long, short = 'm', default_value="200")]
    pub marker_compression_factor: usize,
    /// Minimum percent identity to include pairs 
    #[clap(long, short = 's', default_value="80.0")]
    pub min_percent_identity: f64,
    /// Minimum alignment fraction to include pair
    #[clap(long, short = 'n', default_value="15.0")]
    pub min_alignment_fraction: f64,
    /// Small genomes preset
    #[clap(long, short = 'g')]
    pub small_genomes: bool,
    /// Threads for distance matrix computation
    #[clap(long, short = 't', default_value = "8")]
    pub threads: usize,
}

/// Validator function to check if each file exists and is valid
fn validate_file(file: &str) -> Result<PathBuf, String> {
    let path = PathBuf::from(file);

    if !path.exists() {
        return Err(format!("File not found: {}", file));
    }

    if !path.is_file() {
        return Err(format!("Not a valid file: {}", file));
    }

    Ok(path)
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
