use clap::{Args, Parser, Subcommand};
use std::path::PathBuf;

use crate::{alignment::{Aligner, Preset, SelectHighest}, annotation::AnnotationPreset};

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
    #[arg(short, long, num_args(0..))]
    pub input: Vec<PathBuf>,
    /// Output summary table of coverage and assembly metrics
    #[arg(short, long)]
    pub output: PathBuf,
    /// Reference index for aligner
    ///
    /// Alignment index for 'bowtie2' (index), 'minimap2' (index | fasta) 
    /// and 'strobealign' (index | fasta)
    #[arg(long, short='I')]
    pub index: PathBuf,
    /// Reference sequences used in aligner (fasta)
    #[arg(long, short='R')]
    pub reference: PathBuf,
    /// Aligner
    #[arg(long, short='A', default_value="minimap2")]
    pub aligner: Aligner,
    /// Aligner preset (minimap2)
    #[arg(long, short='P', default_value="sr")]
    pub preset: Option<Preset>,
    /// Group alignments by a field in the reference sequence description
    #[clap(long, short = 'g', long, default_value="taxid=")]
    pub group_by: Option<String>,
    /// Group field separator in the reference sequence description
    #[clap(long, short = 's', long, default_value = ";")]
    pub field_sep: String,
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
    /// Print formatted table to console
    #[clap(long, short = 'T')]
    pub table: bool,
    // Segment field to identify segments in grouped alignments
    #[clap(long, default_value="segment=")]
    pub segment_field: Option<String>,
    /// Segment field identifier negative (e.g. "segment=N/A")
    #[clap(long, default_value="segment=N/A")]
    pub segment_field_nan: Option<String>,
    /// Threads for scanning alignment
    #[clap(long, default_value = "8")]
    pub scan_threads: usize,
    /// Threads for remapping alignment
    #[clap(long, default_value = "2")]
    pub remap_threads: usize,
    /// Include secondary alignments
    #[clap(long)]
    pub secondary: bool,
    /// Minimum remap coverage required for consensus assembly
    #[clap(long, default_value="0.2")]
    pub min_remap_coverage: f64,
    /// Minimum consensus assembly read depth to call a site
    #[clap(long, default_value="10")]
    pub min_consensus_depth: usize,
    /// Minimum consensus frequency to call a variant site
    #[clap(long, default_value="0.75")]
    pub min_consensus_frequency: f64,
    /// Minimum base quality to consider a site
    #[clap(long, default_value="20")]
    pub min_consensus_quality: usize,
    /// Do not create consensus genome from remapping stage
    #[clap(long)]
    pub no_consensus: bool,
    /// Remap all input reads instead of the grouped reads
    #[clap(long)]
    pub remap_all: bool,
}

#[derive(Debug, Args)]
pub struct CoverageArgs {
    /// Input read files (can be compressed with .gz)
    ///
    /// One or two input read files. These files can be in gzipped format.
    /// This parameter is required and multiple files can be specified (1 for long
    /// reads or 2 for paired-end short reads) either consecutively or using multiple
    /// input arguments, for example: '-i R1.fq.gz -i R2.fq.gz' or '-i R1.fq.gz R2.fq.gz'
    #[arg(short, long, num_args(0..))]
    pub input: Vec<PathBuf>,
    /// Output summary table of coverage metrics
    #[arg(short, long)]
    pub output: PathBuf,
    /// Reference index for aligner
    ///
    /// Depending on whether --aligner is chosen, the index is an alignment index 
    /// for 'bowtie2' (index), 'minimap2' and 'strobealign' (index or fasta).
    #[arg(long, short='I')]
    pub index: PathBuf,
    /// Reference sequences used in alignment (required for --zero)
    #[arg(long, short='R')]
    pub reference: Option<PathBuf>,
    /// Aligner
    #[arg(long, short='A', default_value="minimap2")]
    pub aligner: Aligner,
    /// Aligner preset (minimap2)
    #[arg(long, short='P', default_value="sr")]
    pub preset: Option<Preset>,
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
    ConcatOutput(ConcatArgs)
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
    #[clap(long, short = 'm', default_value="0")]
    pub min_completeness: f64,
    /// Add the file stem as identifier to column
    #[clap(long, short = 'f')]
    pub file_id: bool,
}

#[derive(Debug, Args)]
pub struct AnnotateDatabaseArgs {
    /// Vircov database sequence file (.fasta) 
    #[clap(long, short = 'f')]
    pub fasta: PathBuf,
    /// Vircov database sequence annotation file (.tsv)
    /// 
    /// Required columns: id (string, sequence identifier, required value), bin (string, binning variable, required value), 
    /// segment (string, segment annotation, optional value), description (string, sequence description, optional value)
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
