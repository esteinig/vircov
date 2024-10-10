use anyhow::Result;
use indexmap::IndexMap;
use noodles::fasta::record::{Record, Definition};
use regex::Regex;
use serde::{Deserialize, Serialize};
use tabled::settings::{Disable, Style};
use std::collections::{HashMap, HashSet};
use std::num::{ParseFloatError, ParseIntError};
use std::path::{Path, PathBuf};
use tabled::{settings::Width, settings::object::Columns, Table, Tabled};
use thiserror::Error;
use std::io::Write;
use std::fs::{create_dir_all, File};
use std::io::{BufRead, BufReader, BufWriter};
use std::process::Command;
use netview::prelude::*;
use tar::Archive;

/*
========================
Custom error definitions
========================
*/

#[derive(Error, Debug)]
pub enum SubtypeError {}

/*
====================
Output table display
====================
*/

fn display_coverage(cov: &f64) -> String {
    format!("{:.1}%", cov)
}

fn display_optional<T>(cov: &Option<T>) -> String
where
    T: IntoFormattedString,
{
    match cov {
        Some(value) => value.into_formatted_string(),
        None => "n/a".to_string(),
    }
}

trait IntoFormattedString {
    fn into_formatted_string(&self) -> String;
}

// Implement IntoFormattedString for f64 to format as percentage
impl IntoFormattedString for f64 {
    fn into_formatted_string(&self) -> String {
        format!("{:.1}%", self)
    }
}

// Implement IntoFormattedString for usize to format as an integer string
impl IntoFormattedString for usize {
    fn into_formatted_string(&self) -> String {
        format!("{}", self)
    }
}


// Implement IntoFormattedString for String (Option)
impl IntoFormattedString for String {
    fn into_formatted_string(&self) -> String {
        self.to_string()
    }
}

/*
================
Subtype Database
================
*/

pub struct SubtypeDatabaseFiles {
    csv: PathBuf,
    fasta_aa: PathBuf,
    fasta_nuc: PathBuf,
    db_blastn: PathBuf,
    db_blastx: PathBuf,
}


#[derive(Error, Debug)]
pub enum SubtypeDatabaseError {
    #[error("Decompression error")]
    Decompression(#[from] DecompressionError),
    #[error("Subtyping error")]
    Subtyping(#[from] SubtypeError),
    #[error("File name extraction error")]
    FileName(#[from] FileNameError),
    #[error("Database file missing ({0})")]
    DatabaseFileMissing(String),
    #[error(transparent)]
    IOError(#[from] std::io::Error),
    #[error(transparent)]
    NetviewError(#[from] netview::error::NetviewError),
    #[error("Failed to parse line: {0}")]
    LineParse(String),
    #[error("Parse integer error for record")]
    ParseInteger(#[from] ParseIntError),
    #[error("Parse value error for record")]
    ParseValue(String),
    #[error("Parse float error for record")]
    ParseFloat(#[from] ParseFloatError),
    #[error("Failed to execute '{0}' - is it installed?")]
    ProgramExecutionFailed(String),
    #[error("Failed to run command, output is:\n{0}")]
    CommandExecutionFailed(String),
    #[error("Failed to parse `skani` output matrix into symmetrical distance matrix")]
    ParseSkaniMatrix,
    #[error("Failed to serialize record to CSV")]
    SerializeRecord(#[source] csv::Error),
    #[error("Failed to serialize record to CSV")]
    DeserializeRecord(#[source] csv::Error),
    #[error("Regex error: {0}")]
    RegexError(#[from] regex::Error),
    #[error("Failed to delete working directory as it is the current user working directory: {0}")]
    DeleteWorkdir(String),
    #[error("error processing table data (CSV)")]
    CsvWriteError(#[from] csv::Error),
    #[error("failed to provide accession arguments")]
    AccessionArgument,
    #[error("failed to read accession file")]
    AccessionFileError,
    #[error("data mismatch error: {0}")]
    DataMismatchError(String),
}

pub struct SkaniConfig {
    pub compression_factor: usize,
    pub marker_compression_factor: usize,
    pub min_percent_identity: f64,
    pub min_alignment_fraction: f64,
    pub small_genomes: bool,
    pub threads: usize,
}
impl Default for SkaniConfig {
    fn default() -> Self {
        Self {
            compression_factor: 30,
            marker_compression_factor: 200,
            min_percent_identity: 80.0,
            min_alignment_fraction: 15.0,
            small_genomes: true,
            threads: 2
        }
    }
}

pub struct AlignmentConfig {
    blast: BlastConfig,
}

impl Default for AlignmentConfig {
    fn default() -> Self {
        Self {
            blast: BlastConfig::default(),
        }
    }
}

pub struct BlastConfig {
    min_percent_identity: f64,
    max_target_seqs: u64,
    min_evalue: f64,
}
impl Default for BlastConfig {
    fn default() -> Self {
        Self {
            min_percent_identity: 50.,
            max_target_seqs: 25000,
            min_evalue: 0.003
        }
    }
}
pub struct SubtypeDatabase {
    pub name: String,
    pub path: PathBuf,
    pub protein: bool,
    pub genotypes: Vec<Genotype>,
    pub files: SubtypeDatabaseFiles,
    pub skani: SkaniConfig
}

impl SubtypeDatabase {
    pub fn from(
        archive: &PathBuf,
        outdir: &Option<PathBuf>
    ) -> Result<Self, SubtypeDatabaseError> {
        
        log::info!("Unpacking and preparing reference databases...");
        
        let outdir = match outdir {
            Some(outdir) => outdir.to_owned(),
            None => PathBuf::from(chrono::Utc::now().format("%Y-%m-%d").to_string())
        };
        
        create_dir_all(&outdir)?;
        decompress_archive(archive, &outdir)?;

        let genome_fasta = outdir.join("db").join("db_nuc.fasta");
        let protein_fasta = outdir.join("db").join("db_aa.fasta");
        let db_csv = outdir.join("db").join("db.csv");

        if !genome_fasta.exists() {
            return Err(SubtypeDatabaseError::DatabaseFileMissing(
                genome_fasta.display().to_string(),
            ));
        }

        if !db_csv.exists() {
            return Err(SubtypeDatabaseError::DatabaseFileMissing(db_csv.display().to_string()))
        }


        log::info!("Reading genotype database file...");
        let genotypes = read_genotypes_from_file(&db_csv, false)?;
        let name = get_base_file_name(&archive)?;

        let fasta_aa = outdir.join("db").join("db_aa_matched.fasta");
        let db_blastx = outdir.join("db_blastx");
        let db_blastn = outdir.join("db_blastn");

        let db = Self {
            name,
            path: outdir.to_path_buf(),
            protein: protein_fasta.exists(),
            genotypes,
            files: SubtypeDatabaseFiles {
                csv: db_csv,
                fasta_aa: fasta_aa.clone(),
                fasta_nuc: genome_fasta.clone(),
                db_blastn: db_blastn.clone(),
                db_blastx: db_blastx.clone(),
            },
            skani: SkaniConfig::default()
        };

        match protein_fasta.exists() {
            true => {
                db.match_fasta_dbs(&protein_fasta, &genome_fasta, &fasta_aa)?;
                db.makedb_blast_prot(&fasta_aa, &db_blastx)?;
                true
            },
            false => {
                log::warn!("Database does not contain a protein reference");
                false
            }
        };
        db.makedb_blast_nuc(&genome_fasta, &db_blastn)?;
                
        Ok(db)

    }

    pub fn subtype(
        &self,
        fasta: &PathBuf,
        min_cov: f64,
        min_cov_aa: f64,
        min_cov_prot: f64,
        alignment_config: Option<AlignmentConfig>,
        threads: u32
    ) -> Result<Vec<SubtypeSummary>, SubtypeDatabaseError> {

        let input_name = get_base_file_name(&fasta)?;
        let config = alignment_config.unwrap_or(
            AlignmentConfig::default()
        );

        let blastn_output = self.path.join(format!("{}.{}.blastn.tsv", input_name, self.name));
        let blastx_output = self.path.join(format!("{}.{}.blastx.tsv", input_name, self.name));

        log::info!("Computing average nucleotide identity and coverage...");
        self.run_blastn(
            fasta,
            &self.files.db_blastn,
            &blastn_output,
            config.blast.min_percent_identity,
            config.blast.max_target_seqs,
            threads,
        )?;
        
        let ani = self.compute_blastn_ani_cov(&blastn_output)?;
        let aai = match self.protein {
            true => {
                log::info!("Computing average amino acid identity and coverage...");
                self.run_blastx(
                    fasta,
                    &self.files.db_blastx,
                    &blastx_output,
                    config.blast.min_evalue,
                    config.blast.max_target_seqs,
                    threads,
                )?;
                Some(self.compute_blastx_aai_cov(&blastx_output)?)
            },
            false => None
        };
        
        log::info!("Computing mutual nearest neighbors graph...");

        let combined_fasta = self.path.join(format!("{}.genomes.fasta", input_name));

        self.compute_mknn_graph(fasta, &combined_fasta)?;

       
        let subtype_data = SubtypeSummary::from(
            &input_name, ani, aai, Some(self.genotypes.clone())
        );

        log::debug!("Subtype summaries main: {:#?}", subtype_data.len());

        log::info!("Coverage filters: min_cov={} min_cov_aa={} min_cov_prot={}", min_cov, min_cov_aa, min_cov_prot);
        let subtype_data = filter_subtype_summaries(
            &subtype_data, min_cov, min_cov_aa, min_cov_prot
        )?;
        
        log::debug!("Subtype summaries filtered: {}", subtype_data.len());
        
        let table_output = self.path.join(format!("{input_name}.filtered.tsv"));
        write_subtype_summaries_to_csv(&subtype_data, &table_output, true)?;

        Ok(subtype_data)
    }
    pub fn compute_mknn_graph(&self, fasta: &PathBuf, combined_fasta: &PathBuf) -> Result<(), SubtypeDatabaseError> {

        // // Concatentate the input genome with the reference nucleotide data 
        // log::info!("Concatenating query sequences to database sequences....");
        // concatenate_fasta_files(&self.files.fasta_nuc, &[fasta], &combined_fasta)?;

        // let netview = Netview::new();

        // let (distance, af, ids) = netview.skani_distance(
        //     fasta, 
        //     self.skani.marker_compression_factor, 
        //     self.skani.compression_factor, 
        //     self.skani.threads, 
        //     self.skani.min_percent_identity, 
        //     self.skani.min_alignment_fraction, 
        //     self.skani.small_genomes
        // )?;

        // let graph = netview.graph_from_vecs(distance, 30, Some(af), Some(ids))?;

        


        Ok(())
    }
    pub fn create_ranked_tables(&self, output: &PathBuf, subtype_summaries: &Vec<SubtypeSummary>, metric: &str, ranks: usize, print_ranks: bool, with_genotype: bool, protein: bool) -> Result<(), SubtypeDatabaseError> {

        let subtype_summaries = match with_genotype {
            true => subtype_summaries.into_iter().filter(|s| {
                s.target_genotype.is_some()
            }).cloned().collect(),
            false => subtype_summaries.clone()
        };

        let file_groups = group_and_sort_summaries_by_file(&subtype_summaries, metric)?;

        let mut output_rows = Vec::new();
        for file_summary in file_groups {
            let query_groups = group_and_sort_summaries_by_query(&file_summary.summaries, metric)?;
            

            let mut query_subsets_ranks = Vec::new();
            for query_summary in query_groups {
                
                let mut ranked_subset: Vec<SubtypeSummary> = query_summary.summaries.iter().take(ranks).map(|s| s.clone()).collect();

                query_subsets_ranks.append(&mut ranked_subset);
            }

            let ranked_file_subset_path = self.path.join(format!("{}.ranks.tsv", file_summary.file));
            write_subtype_summaries_to_csv(&query_subsets_ranks, &ranked_file_subset_path, true)?;
            

            if print_ranks {
                let mut subtable = Table::new(&query_subsets_ranks);
                subtable.modify(Columns::new(1..3), Width::truncate(match protein { true => 16, false => 32}).suffix("...")).with(Style::modern());

                if !protein {
                    subtable.with(Disable::column(Columns::new(3..9)));
                }

                println!("{}", subtable);
            }
            output_rows.push(query_subsets_ranks);

        }

        let rows: Vec<SubtypeSummary> = output_rows.iter().flatten().cloned().collect();
        write_subtype_summaries_to_csv(&rows, &output, true)?;

        Ok(())
    }
    pub fn remove_workdir(&self) -> Result<(), SubtypeDatabaseError> {
        if self.path == std::env::current_dir()? {
            return Err(SubtypeDatabaseError::DeleteWorkdir(self.path.display().to_string()))
        };
        std::fs::remove_dir_all(&self.path)?;
        Ok(())
    }
    pub fn compute_blastn_ani_cov(
        &self,
        blastn_results: &PathBuf,
    ) -> Result<Vec<BlastAniSummary>, SubtypeDatabaseError> {
        let file = File::open(blastn_results)?;
        let reader = BufReader::new(file);
        let parsed_lines = parse_blast(reader.lines().filter_map(Result::ok));
        let alignment_blocks = collect_alignment_blocks(parsed_lines)?;
        let summaries = process_alignments(alignment_blocks)?;

        Ok(summaries)
    }
    pub fn compute_blastx_aai_cov(
        &self,
        blastx_results: &PathBuf,
    ) -> Result<Vec<BlastAaiSummary>, SubtypeDatabaseError> {
        let genome_info = parse_aa_fasta(&self.files.fasta_aa)?;
        let genome_alns = parse_and_organize_alignments(blastx_results, &genome_info)?;
        let summaries = compute_summaries(&genome_info, &genome_alns);

        Ok(summaries)
    }
    pub fn match_fasta_dbs(&self, protein_fasta: &PathBuf, genome_fasta: &PathBuf, protein_fasta_renamed: &PathBuf) -> Result<(), SubtypeDatabaseError> {
        
        let proteins = parse_fasta_annotations_file(&protein_fasta)?;
        let genomes = parse_fasta_annotations_file(genome_fasta)?;
        let matches = match_accessions(&proteins, &genomes);

        log::debug!("{:#?}", matches);

        let output_file = File::create(protein_fasta_renamed)?;
        let mut writer = noodles::fasta::Writer::new(BufWriter::new(output_file));
        
        let file = File::open(&protein_fasta)?;
        let reader = BufReader::new(file);
        let mut fasta_reader = noodles::fasta::Reader::new(reader);


        let mut not_matched = 0;
        let mut genome_counts = HashMap::new();
        for (_, result) in fasta_reader.records().enumerate() { 
            let record = result?;

            if let Some(genome_accession) = matches.get(&record.name().to_string()) {

                let protein_in_genome_count = genome_counts.get(&genome_accession).unwrap_or(&0);
                let new_id = format!("{}__{}", genome_accession, protein_in_genome_count); // Replace with matching genome accession

                genome_counts.entry(genome_accession)
                    .and_modify(|e| { *e += 1 })
                    .or_insert(1);

                let header = Definition::new(new_id, None);
                let new_record = noodles::fasta::Record::new(header, record.sequence().clone());

                writer.write_record(&new_record)?;

            } else {
                // Unmachted proteins are not written to the matching database file!
                not_matched += 1;
            }
            
        }

        log::debug!("Matched {} proteins to {} genomes with {} proteins unmatched to any genome", proteins.len(), genomes.len(), not_matched);

        Ok(())
    }
    pub fn makedb_blast_nuc(
        &self,
        db_fasta: &PathBuf,
        output: &PathBuf,
    ) -> Result<(), SubtypeDatabaseError> {
        let args = vec![
            "-in".to_string(),
            db_fasta.display().to_string(),
            "-out".to_string(),
            output.display().to_string(),
            "-dbtype".to_string(),
            "nucl".to_string(),
        ];
        log::debug!("Running command: makeblastdb {}", &args.join(" "));
        run_command(args, "makeblastdb")
    }
    pub fn makedb_blast_prot(
        &self,
        db_fasta: &PathBuf,
        output: &PathBuf,
    ) -> Result<(), SubtypeDatabaseError> {
        let args = vec![
            "-in".to_string(),
            db_fasta.display().to_string(),
            "-out".to_string(),
            output.display().to_string(),
            "-dbtype".to_string(),
            "prot".to_string(),
        ];
        log::debug!("Running command: makeblastdb {}", &args.join(" "));
        run_command(args, "makeblastdb")
    }
    pub fn run_blastn(
        &self,
        input: &PathBuf,
        database: &PathBuf,
        output: &PathBuf,
        min_percent_identity: f64,
        max_target_seqs: u64,
        threads: u32,
    ) -> Result<(), SubtypeDatabaseError> {
        let args = vec![
            "-query".to_string(),
            input.display().to_string(),
            "-db".to_string(),
            database.display().to_string(),
            "-out".to_string(),
            output.display().to_string(),
            "-outfmt".to_string(),
            "6 std qlen slen".to_string(), // No need for single quotes around the format
            "-max_target_seqs".to_string(),
            max_target_seqs.to_string(),
            "-num_threads".to_string(),
            threads.to_string(),
            "-perc_identity".to_string(),
            min_percent_identity.to_string(),
        ];

        log::debug!("Running command: blastn with args {}", &args.join(" "));
        run_command(args, "blastn")
    }
    pub fn run_blastx(
        &self,
        input: &PathBuf,
        database: &PathBuf,
        output: &PathBuf,
        max_evalue: f64,
        max_target_seqs: u64,
        threads: u32,
    ) -> Result<(), SubtypeDatabaseError> {
        let args = vec![
            "-query".to_string(),
            input.display().to_string(),
            "-db".to_string(),
            database.display().to_string(),
            "-out".to_string(),
            output.display().to_string(),
            "-outfmt".to_string(),
            "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore"
                .to_string(),
            "-evalue".to_string(),
            max_evalue.to_string(),
            "-max_target_seqs".to_string(),
            max_target_seqs.to_string(),
            "-num_threads".to_string(),
            threads.to_string(),
        ];
        log::debug!("Running command: blastx {}", &args.join(" "));
        run_command(args, "blastx")
    }
}


/// Represents a selected group of subtype summaries for a specific query.
#[derive(Debug, Serialize)]
pub struct SubtypeFileSelection {
    file: String,
    summaries: Vec<SubtypeSummary>,
}

pub fn group_and_sort_summaries_by_file(
    summaries: &Vec<SubtypeSummary>,
    metric: &str
) -> Result<Vec<SubtypeFileSelection>, SubtypeDatabaseError> {
    let mut map: HashMap<String, Vec<SubtypeSummary>> = HashMap::new();

    // Grouping summaries by 'file'
    for summary in summaries {
        map.entry(summary.file.clone())
            .or_default()
            .push(summary.clone());
    }

    // Sorting groups by metric descending and wrapping in 'SubtypeSelection'
    let mut results = Vec::new();
    for (file, mut summaries) in map {
        summaries.sort_by(|a, b| {
            match metric {
                "ani" => b.ani.partial_cmp(&a.ani).unwrap_or(std::cmp::Ordering::Equal),
                "wani" => b.ani_weighted.partial_cmp(&a.ani_weighted).unwrap_or(std::cmp::Ordering::Equal),
                "aai" => b.aai.partial_cmp(&a.aai).unwrap_or(std::cmp::Ordering::Equal),
                "waai" => b.aai_weighted.partial_cmp(&a.aai_weighted).unwrap_or(std::cmp::Ordering::Equal),
                _ => b.ani.partial_cmp(&a.ani).unwrap_or(std::cmp::Ordering::Equal) // defualt to ANI
            }
            
        });
        results.push(SubtypeFileSelection { file, summaries });
    }

    Ok(results)
}

/// Represents a selected group of subtype summaries for a specific query.
#[derive(Debug, Serialize)]
pub struct SubtypeQuerySelection {
    query: String,
    summaries: Vec<SubtypeSummary>,
}

pub fn group_and_sort_summaries_by_query(
    summaries: &Vec<SubtypeSummary>,
    metric: &str
) -> Result<Vec<SubtypeQuerySelection>, SubtypeDatabaseError> {
    let mut map: HashMap<String, Vec<SubtypeSummary>> = HashMap::new();

    // Grouping summaries by 'file'
    for summary in summaries {
        map.entry(summary.query.clone())
            .or_default()
            .push(summary.clone());
    }

    // Sorting groups by metric descending and wrapping in 'SubtypeQuerySelection'
    let mut results = Vec::new();
    for (query, mut summaries) in map {
        summaries.sort_by(|a, b| {
            match metric {
                "ani" => b.ani.partial_cmp(&a.ani).unwrap_or(std::cmp::Ordering::Equal),
                "wani" => b.ani_weighted.partial_cmp(&a.ani_weighted).unwrap_or(std::cmp::Ordering::Equal),
                "aai" => b.aai.partial_cmp(&a.aai).unwrap_or(std::cmp::Ordering::Equal),
                "waai" => b.aai_weighted.partial_cmp(&a.aai_weighted).unwrap_or(std::cmp::Ordering::Equal),
                _ => b.ani.partial_cmp(&a.ani).unwrap_or(std::cmp::Ordering::Equal) // defualt to ANI
            }
            
        });
        results.push(SubtypeQuerySelection { query, summaries });
    }

    Ok(results)
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Genotype {
    accession: String,
    genotype: Option<String>,
    segment: Option<String>,
}

/// Parses genotype data from a CSV or TSV file into a vector of `Genotype`.
///
/// # Arguments
///
/// * `file_path` - Path to the CSV or TSV file containing genotype data.
/// * `tsv` - Whether the file is formatted as TSV instead of CSV.
///
/// # Example
///
/// ```no_run
/// let genotypes = read_genotypes_from_file("path/to/genotypes.csv", false);
/// match genotypes {
///     Ok(genotypes) => println!("Read {} genotype entries", genotypes.len()),
///     Err(e) => eprintln!("Failed to read genotype data: {}", e),
/// }
/// ```
///
/// # Errors
///
/// Returns `SubtypeDatabaseError` if there is an error opening the file or parsing its contents.
pub fn read_genotypes_from_file<P: AsRef<Path>>(
    file_path: P,
    tsv: bool,
) -> Result<Vec<Genotype>, SubtypeDatabaseError> {
    let file = File::open(file_path)?;
    let mut rdr = if tsv {
        csv::ReaderBuilder::new().delimiter(b'\t').trim(csv::Trim::All).from_reader(file)
    } else {
        csv::ReaderBuilder::new().trim(csv::Trim::All).from_reader(file)
    };

    rdr.deserialize()
        .collect::<Result<Vec<Genotype>, csv::Error>>()
        .map_err(SubtypeDatabaseError::DeserializeRecord)
}


/// Updates a vector of `SubtypeSummary` with corresponding `Genotype` information based on matching accession fields.
///
/// # Arguments
///
/// * `summaries` - A mutable reference to a vector of `SubtypeSummary`.
/// * `genotypes` - A reference to a vector of `Genotype`.
///
/// # Example
///
/// ```no_run
/// let mut summaries = vec![SubtypeSummary {
///     query: "Q1".to_string(),
///     target: "T1".to_string(),
///     target_proteins: Some(100),
///     shared_proteins: Some(50),
///     aai_tcov_protein: Some(99.9),
///     aai_tcov_aa: Some(98.5),
///     aai: Some(99.1),
///     aai_weighted: Some(99.0),
///     ani_qcov: Some(98.8),
///     ani_tcov: Some(97.5),
///     ani: Some(98.0),
///     ani_weighted: Some(98.2),
///     target_segment: None,
///     target_genotype: None,
/// }];
///
/// let genotypes = vec![Genotype {
///     accession: "T1".to_string(),
///     genotype: Some("G1".to_string()),
///     segment: Some("S1".to_string()),
/// }];
///
/// update_subtype_summaries_with_genotypes(&mut summaries, &genotypes);
/// ```
///
/// # Note
///
/// This function does not return any value but updates the `summaries` vector in place.
pub fn update_subtype_summaries_with_genotypes(
    summaries: &mut Vec<SubtypeSummary>,
    genotypes: &[Genotype],
) {
    let genotype_map: HashMap<String, &Genotype> = genotypes
        .iter()
        .map(|g| (g.accession.clone(), g))
        .collect();

    for summary in summaries.iter_mut() {
        if let Some(genotype) = genotype_map.get(&summary.target) {
            summary.target_genotype = genotype.genotype.clone();
            summary.target_segment = genotype.segment.clone();
        }
    }
}

#[derive(Debug, Clone, Tabled, Serialize)]
pub struct SubtypeSummary {
    #[tabled(rename="File")]
    file: String,
    #[tabled(rename="Query")]
    query: String,
    #[tabled(rename="Target")]
    target: String,
    #[tabled(rename="Target proteins")]
    #[tabled(display_with = "display_optional")]
    target_proteins: Option<usize>,
    #[tabled(rename="Shared proteins")]
    #[tabled(display_with = "display_optional")]
    shared_proteins: Option<usize>,
    #[tabled(rename="Protein coverage")]
    #[tabled(display_with = "display_optional")]
    aai_tcov_protein: Option<f64>,
    #[tabled(rename="AA query coverage")]
    #[tabled(display_with = "display_optional")]
    aai_qcov_aa: Option<f64>,
    #[tabled(rename="AA target coverage")]
    #[tabled(display_with = "display_optional")]
    aai_tcov_aa: Option<f64>,
    #[tabled(rename="AAI")]
    #[tabled(display_with = "display_optional")]
    aai: Option<f64>,
    #[tabled(rename="wAAI")]
    #[tabled(display_with = "display_optional")]
    aai_weighted: Option<f64>,
    #[tabled(rename="Query coverage")]
    #[tabled(display_with = "display_optional")]
    ani_qcov: Option<f64>,
    #[tabled(rename="Target coverage")]
    #[tabled(display_with = "display_optional")]
    ani_tcov: Option<f64>,
    #[tabled(rename="ANI")]
    #[tabled(display_with = "display_optional")]
    ani: Option<f64>,
    #[tabled(rename="wANI")]
    #[tabled(display_with = "display_optional")]
    ani_weighted: Option<f64>,    
    #[tabled(rename="Segment")]
    #[tabled(display_with = "display_optional")]
    target_segment: Option<String>,
    #[tabled(rename="Genotype")]
    #[tabled(display_with = "display_optional")]
    target_genotype: Option<String>,
}

impl SubtypeSummary {
    pub fn from(input_name: &str, ani_summary: Vec<BlastAniSummary>, aai_summary: Option<Vec<BlastAaiSummary>>, genotypes: Option<Vec<Genotype>>) -> Vec<Self> {

        let mut summaries = Vec::new();
        let mut unmatched_aais = if let Some(aai) = &aai_summary { aai.clone() } else { Vec::new() };

        for ani in ani_summary.iter() {
            let matched_aai = aai_summary.as_ref()
                .and_then(|aais| aais.iter().find(|aai| aai.query == ani.qname && aai.target == ani.tname));

            if let Some(aai) = matched_aai {
                summaries.push(SubtypeSummary {
                    file: input_name.to_string(),
                    query: ani.qname.clone(),
                    target: ani.tname.clone(),
                    target_proteins: Some(aai.target_proteins),
                    shared_proteins: Some(aai.shared_proteins),
                    ani: Some(ani.ani),
                    ani_weighted: Some(ani.ani_weighted),
                    ani_qcov: Some(ani.qcov),
                    ani_tcov: Some(ani.tcov),
                    aai: Some(aai.aai),
                    aai_tcov_protein: Some(aai.protein_tcov),
                    aai_qcov_aa: Some(aai.aa_qcov),
                    aai_tcov_aa: Some(aai.aa_tcov),
                    aai_weighted: Some(aai.weighted_aai),
                    target_genotype: None,
                    target_segment: None
                });
                unmatched_aais.retain(|u| !(u.query == ani.qname && u.target == ani.tname));
            } else {
                summaries.push(SubtypeSummary {
                    file: input_name.to_string(),
                    query: ani.qname.clone(),
                    target: ani.tname.clone(),
                    target_proteins: None,
                    shared_proteins: None,
                    ani: Some(ani.ani),
                    ani_weighted: Some(ani.ani_weighted),
                    ani_qcov: Some(ani.qcov),
                    ani_tcov: Some(ani.tcov),
                    aai: None,
                    aai_tcov_protein: None,
                    aai_qcov_aa: None,
                    aai_tcov_aa: None,
                    aai_weighted: None,
                    target_genotype: None,
                    target_segment: None
                });
            }
        }

        // Add entries for AAI summaries not matched with any ANI summary
        for aai in unmatched_aais {
            summaries.push(SubtypeSummary {
                file: input_name.to_string(),
                query: aai.query.clone(),
                target: aai.target.clone(),
                target_proteins: Some(aai.target_proteins),
                shared_proteins: Some(aai.shared_proteins),
                ani: None,
                ani_weighted: None,
                ani_qcov: None,
                ani_tcov: None,
                aai: Some(aai.aai),
                aai_tcov_protein: Some(aai.protein_tcov),
                aai_qcov_aa: Some(aai.aa_qcov),
                aai_tcov_aa: Some(aai.aa_tcov),
                aai_weighted: Some(aai.weighted_aai),
                target_genotype: None,
                target_segment: None
            });
        }

        if let Some(genotypes) = genotypes {
            update_subtype_summaries_with_genotypes(&mut summaries, &genotypes);
        }

        summaries
    }
}

/// Represents a summary of BLAST ANI and coverage results.
#[derive(Debug, Clone, Tabled)]
pub struct BlastAniSummary {
    #[tabled(rename="Query")]
    pub qname: String,
    #[tabled(rename="Target")]
    pub tname: String,
    #[tabled(rename="Alignments")]
    pub aln_count: usize,
    #[tabled(rename="Query coverage")]
    pub qcov: f64,
    #[tabled(rename="Target coverage")]
    pub tcov: f64,
    #[tabled(rename="ANI")]
    pub ani: f64,
    #[tabled(rename="ANIw")]
    pub ani_weighted: f64,
}

#[derive(Debug, Clone, Tabled)]
pub struct BlastAaiSummary {
    #[tabled(rename="Query")]
    query: String,
    #[tabled(rename="Target")]
    target: String,
    #[tabled(rename="Target proteins")]
    target_proteins: usize,
    #[tabled(rename="Shared proteins")]
    shared_proteins: usize,
    // #[tabled(rename="Query protein coverage")]
    // #[tabled(display_with = "display_optional")]
    // protein_qcov: Option<f64>,
    #[tabled(rename="Protein coverage")]
    #[tabled(display_with = "display_coverage")]
    protein_tcov: f64,
    #[tabled(rename="AA query coverage")]
    #[tabled(display_with = "display_coverage")]
    aa_qcov: f64,
    #[tabled(rename="AA target coverage")]
    #[tabled(display_with = "display_coverage")]
    aa_tcov: f64,
    #[tabled(rename="AAI")]
    #[tabled(display_with = "display_coverage")]
    aai: f64,
    #[tabled(rename="AAIw")]
    #[tabled(display_with = "display_coverage")]
    weighted_aai: f64,
}


#[derive(Debug)]
pub struct BlastRecord {
    qname: String,
    tname: String,
    pid: f64,
    len: f64,
    qcoords: (i64, i64),
    tcoords: (i64, i64),
    qlen: f64,
    tlen: f64,
    evalue: f64,
}

#[derive(Debug)]
pub struct Alignment {
    qseqid: String,
    sseqid: String,
    pident: f64,
    _length: usize,
    qlen: usize,
    slen: usize,
    qstart: usize,
    qend: usize,
    sstart: usize,
    send: usize,
    _evalue: f64,
    _bitscore: f64,
    qcov: f64,
    scov: f64,
}

#[derive(Debug, Clone)]
pub struct GenomeInfo {
    protein_to_genome: HashMap<String, String>,
    num_proteins: HashMap<String, usize>,
}


/// Processes blocks of alignments to compute ANI and coverage,
/// filtering alignments based on length and e-value.
///
/// # Arguments
///
/// * `alignment_blocks` - A vec of alignment blocks (`BlastRecord`) to be processed.
///
/// # Returns
///
/// This function returns a `Result` containing a vector of `BlastAniSummary`,
/// or a `SubtypeDatabaseError` in case of failure.
pub fn process_alignments(
    alignment_blocks: Vec<Vec<BlastRecord>>,
) -> Result<Vec<BlastAniSummary>, SubtypeDatabaseError> {
    let mut summaries = Vec::new();

    for alns in alignment_blocks {

        let pruned_alns = prune_alns(
            alns, 0.0, 1e-3
        );
        if pruned_alns.is_empty() {
            continue;
        }

        let ani = compute_ani(&pruned_alns, false);
        let (qcov, tcov) = compute_cov(&pruned_alns, false);

        summaries.push(BlastAniSummary {
            qname: pruned_alns[0].qname.clone(),
            tname: pruned_alns[0].tname.clone(),
            aln_count: pruned_alns.len(),
            ani,
            qcov: qcov * 100.0,
            tcov: tcov * 100.0,
            ani_weighted: ani*tcov*qcov
        });
    }

    log::debug!("{:#?}", summaries);

    Ok(summaries)
}

/// Parses a BLASTX output file to extract alignments meeting specific coverage criteria.
///
/// This function reads a BLASTX output file and filters alignments based on minimum query
/// and target coverage percentages. Each line in the BLASTX output is expected to be a 
/// tab-delimited alignment record. Records not meeting the specified coverage thresholds 
/// or lacking sufficient data fields are excluded.
///
/// # Arguments
///
/// * `file_path` - A `&PathBuf` reference to the BLASTX output file.
/// * `min_qcov` - The minimum query coverage percentage required for an alignment to be included.
/// * `min_tcov` - The minimum target coverage percentage required for an alignment to be included.
///
/// # Returns
///
/// A `Result<Vec<Alignment>, SubtypeDatabaseError>`, where `Vec<Alignment>` contains the 
/// alignments meeting the specified coverage criteria, or `SubtypeDatabaseError` on failure.
///
/// # Errors
///
/// Returns `SubtypeDatabaseError::IOError` if the file cannot be opened or read.
/// Returns `SubtypeDatabaseError::ParseFloat` or `SubtypeDatabaseError::ParseInteger` 
/// if parsing of numeric fields fails.
///
/// # Example
///
/// ```
/// use std::path::PathBuf;
/// use vircov::{parse_blastx_output, Alignment, SubtypeDatabaseError};
///
/// let file_path = PathBuf::from("blastx_output.txt");
/// let min_qcov = 70.0;
/// let min_tcov = 70.0;
///
/// match parse_blastx_output(&file_path, min_qcov, min_tcov) {
///     Ok(alignments) => println!("Alignments found: {}", alignments.len()),
///     Err(error) => eprintln!("Failed to parse BLASTX output: {:?}", error),
/// }
/// ```
pub fn parse_blastx_output(
    file_path: &PathBuf,
    min_qcov: f64,
    min_tcov: f64,
) -> Result<Vec<Alignment>, SubtypeDatabaseError> {
    let file = File::open(file_path).map_err(SubtypeDatabaseError::IOError)?;
    let reader = BufReader::new(file);
    let mut alignments: Vec<Alignment> = Vec::new();

    for line in reader.lines() {
        let line = line.map_err(SubtypeDatabaseError::IOError)?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 12 {
            continue; // Skip lines that do not have enough fields.
        }

        let mut alignment = Alignment {
            qseqid: parts[0].to_string(),
            sseqid: parts[1].to_string(),
            pident: parts[2].parse().map_err(SubtypeDatabaseError::ParseFloat)?,
            _length: parts[3]
                .parse()
                .map_err(SubtypeDatabaseError::ParseInteger)?,
            qlen: parts[4]
                .parse()
                .map_err(SubtypeDatabaseError::ParseInteger)?,
            slen: parts[5]
                .parse()
                .map_err(SubtypeDatabaseError::ParseInteger)?,
            qstart: parts[6]
                .parse()
                .map_err(SubtypeDatabaseError::ParseInteger)?,
            qend: parts[7]
                .parse()
                .map_err(SubtypeDatabaseError::ParseInteger)?,
            sstart: parts[8]
                .parse()
                .map_err(SubtypeDatabaseError::ParseInteger)?,
            send: parts[9]
                .parse()
                .map_err(SubtypeDatabaseError::ParseInteger)?,
            _evalue: parts[10]
                .parse()
                .map_err(SubtypeDatabaseError::ParseFloat)?,
            _bitscore: parts[11]
                .parse()
                .map_err(SubtypeDatabaseError::ParseFloat)?,
            qcov: 0.0, // To be calculated
            scov: 0.0, // To be calculated
        };

        // Calculate coverages
        alignment.qcov =
            (alignment.qend - alignment.qstart + 1) as f64 / alignment.qlen as f64 * 100.0;
        alignment.scov =
            (alignment.send - alignment.sstart + 1) as f64 / alignment.slen as f64 * 100.0;

        if alignment.qcov >= min_qcov && alignment.scov >= min_tcov {
            alignments.push(alignment);
        }
    }

    Ok(alignments)
}

/// Parses a FASTA file of amino acid sequences to extract protein-to-genome mappings and count the number of proteins per genome.
///
/// This function reads a FASTA file where each entry's header starts with '>' followed by the protein ID,
/// which is expected to contain the genome ID as a prefix, separated from the protein-specific part by "__".
/// It collects information about which proteins belong to which genome and counts the total number of proteins per genome.
///
/// # Arguments
///
/// * `file_path` - A `&PathBuf` reference to the FASTA file containing the amino acid sequences.
///
/// # Returns
///
/// A `Result<GenomeInfo, SubtypeDatabaseError>`, where `GenomeInfo` is a struct containing two fields:
/// `protein_to_genome`, a mapping from protein IDs to genome IDs, and `num_proteins`, a count of proteins per genome.
/// Returns `SubtypeDatabaseError::IOError` on file reading errors.
///
/// # Example
///
/// ```
/// use std::path::PathBuf;
/// use vircov::{parse_aa_fasta, GenomeInfo, SubtypeDatabaseError};
///
/// let file_path = PathBuf::from("path/to/your/aa_sequences.fasta");
///
/// match parse_aa_fasta(&file_path) {
///     Ok(info) => {
///         println!("Number of proteins per genome: {:?}", info.num_proteins);
///         println!("Protein to genome mapping: {:?}", info.protein_to_genome);
///     },
///     Err(error) => eprintln!("Error parsing FASTA file: {:?}", error),
/// }
/// ```
pub fn parse_aa_fasta(file_path: &PathBuf) -> Result<GenomeInfo, SubtypeDatabaseError> {
    let file = match File::open(file_path) {
        Ok(file) => file,
        Err(e) => return Err(SubtypeDatabaseError::IOError(e)),
    };
    let reader = BufReader::new(file);
    let mut info = GenomeInfo {
        protein_to_genome: HashMap::new(),
        num_proteins: HashMap::new(),
    };

    for line in reader.lines() {
        let line = line.map_err(SubtypeDatabaseError::IOError)?;

        if line.starts_with('>') {
            let protein_id = line.trim_start_matches('>').to_owned();
            let ids: Vec<&str> = protein_id.split("__").collect();

            let genome_id = ids[0].to_string();

            *info.num_proteins.entry(genome_id.clone()).or_insert(0) += 1;

            info.protein_to_genome.insert(protein_id, genome_id);
        }
    }

    Ok(info)
}

/// Parses BLASTX output and organizes alignments by genome and target protein.
///
/// This function first parses the BLASTX output file to extract alignments, then organizes these alignments
/// into a nested `HashMap` structure. The top-level `HashMap` keys are genome IDs from the query sequences. 
/// Each genome ID maps to another `HashMap`, where the keys are target genome IDs derived from the protein IDs
/// in the `GenomeInfo` mapping, and the values are vectors of `Alignment` objects representing the alignments
/// of the query genome to proteins of the target genome.
///
/// # Arguments
///
/// * `file_path` - A reference to a `PathBuf` specifying the path to the BLASTX output file.
/// * `info` - A reference to a `GenomeInfo` struct containing mappings from protein IDs to genome IDs, used to organize the alignments.
///
/// # Returns
///
/// A `Result<HashMap<String, HashMap<String, Vec<Alignment>>>, SubtypeDatabaseError>`, where the outer `HashMap` maps
/// genome IDs to another `HashMap` that maps target genome IDs to vectors of `Alignment` objects. Returns `SubtypeDatabaseError`
/// in case of parsing failures or if the target protein is not found in the `GenomeInfo`.
///
/// # Example
///
/// ```no_run
/// use std::path::PathBuf;
/// use vircov::{parse_and_organize_alignments, GenomeInfo, SubtypeDatabaseError, Alignment};
///
/// let file_path = PathBuf::from("path/to/blastx/output.txt");
/// let info = GenomeInfo {
///     // Assuming these fields are populated accordingly
///     protein_to_genome: HashMap::new(),
///     num_proteins: HashMap::new(),
/// };
///
/// match parse_and_organize_alignments(&file_path, &info) {
///     Ok(genome_alns) => println!("Organized alignments: {:?}", genome_alns),
///     Err(error) => eprintln!("Failed to parse and organize alignments: {:?}", error),
/// }
/// ```
pub fn parse_and_organize_alignments(
    file_path: &PathBuf,
    info: &GenomeInfo,
) -> Result<HashMap<String, HashMap<String, Vec<Alignment>>>, SubtypeDatabaseError> {
    let mut genome_alns: HashMap<String, HashMap<String, Vec<Alignment>>> = HashMap::new();
    let alignments = parse_blastx_output(file_path, 0., 0.)?;

    for aln in alignments {
        let genome = aln.qseqid.clone();
        let protein = aln.sseqid.clone();

        let target = info
            .protein_to_genome
            .get(&protein)
            .ok_or(SubtypeDatabaseError::ParseValue(format!(
                "Target protein not found: {}",
                &protein
            )))?
            .to_string();
        genome_alns
            .entry(genome)
            .or_insert_with(HashMap::new)
            .entry(target)
            .or_insert_with(Vec::new)
            .push(aln);
    }

    Ok(genome_alns)
}

/// Computes summary statistics for genome alignments against target genomes.
///
/// This function aggregates alignment data to compute summary statistics for each query genome
/// against each target genome. Summary statistics include the total number of proteins in the target,
/// the number of shared proteins (based on gene IDs), the percentage of shared proteins relative to the
/// total target proteins (protein target coverage), average amino acid identity (AAI), and weighted AAI,
/// which factors in both the protein target coverage and the amino acid target coverage.
///
/// # Arguments
///
/// * `genome_info` - A reference to a `GenomeInfo` struct containing mappings of proteins to genomes
///   and the number of proteins per genome.
/// * `genome_alns` - A reference to a nested `HashMap` where the outer key is the query genome ID,
///   the inner key is the target genome ID, and the value is a vector of `Alignment` objects
///   representing the alignments between the query genome and proteins of the target genome.
///
/// # Returns
///
/// Returns a `Vec<BlastAaiSummary>`, where each `BlastAaiSummary` contains AAI summary statistics for
/// a query genome aligned against a target genome.
///
/// # Example
///
/// ```no_run
/// use std::collections::HashMap;
/// use vircov::{compute_summaries, GenomeInfo, Alignment, BlastAaiSummary};
///
/// let genome_info = GenomeInfo {
///     // Assuming these fields are populated accordingly
///     protein_to_genome: HashMap::new(),
///     num_proteins: HashMap::new(),
/// };
///
/// let genome_alns: HashMap<String, HashMap<String, Vec<Alignment>>> = HashMap::new(); // Assuming populated
///
/// let summaries = compute_summaries(&genome_info, &genome_alns);
/// ```
pub fn compute_summaries(
    genome_info: &GenomeInfo,
    genome_alns: &HashMap<String, HashMap<String, Vec<Alignment>>>,
) -> Vec<BlastAaiSummary> {
    let mut summaries = Vec::new();

    for (query, targets) in genome_alns {
        for (target, alignments) in targets {
            let total_target_proteins = *genome_info.num_proteins.get(target).unwrap_or(&0);

            let mut aln_by_gene = HashMap::new();
            for alignment in alignments {
                
                let ids: Vec<&str> = alignment.sseqid.split("__").collect();
                let gene_id = ids[1];

                aln_by_gene
                    .entry(gene_id)
                    .or_insert_with(Vec::new)
                    .push(alignment);
            }

            let mut shared_proteins = HashSet::new();
            for (gene_id, alignments) in aln_by_gene.into_iter() {
                let (_, _) = compute_coverage_single_target(alignments);
                shared_proteins.insert(gene_id);
            }

            let shared_proteins = shared_proteins.len();
            let protein_tcov = 100.0 * shared_proteins as f64 / total_target_proteins as f64;
            let (aa_qcov, aa_tcov) = compute_coverage(alignments);

            let aai =
                alignments.iter().map(|aln| aln.pident).sum::<f64>() / alignments.len() as f64;

            summaries.push(BlastAaiSummary {
                query: query.clone(),
                target: target.clone(),
                target_proteins: total_target_proteins,
                shared_proteins,
                // protein_qcov: None, // full genome input for now
                protein_tcov,
                aa_qcov,
                aa_tcov,
                aai,
                weighted_aai: aai * (aa_qcov / 100.) * (aa_tcov / 100.),
            });
        }
    }

    summaries
}

/// Merges overlapping or consecutive intervals for single target.
///
/// Given a vector of intervals `(usize, usize)`, where each tuple represents a start and end, this function merges
/// all overlapping or consecutive intervals into the smallest number of non-overlapping intervals.
///
/// # Arguments
///
/// * `intervals` - A vector of tuples where each tuple contains the start and end of an interval.
///
/// # Returns
///
/// Returns a vector of tuples representing the merged intervals.
///
/// # Example
///
/// ```
/// use vircov::merge_intervals_single_target;
///
/// let intervals = vec![(1, 3), (2, 4), (5, 7), (6, 8)];
/// let merged_intervals = merge_intervals_single_target(intervals);
/// assert_eq!(merged_intervals, vec![(1, 4), (5, 8)]);
/// ```
pub fn merge_intervals_single_target(intervals: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
    if intervals.is_empty() {
        return vec![];
    }
    let mut merged: Vec<(usize, usize)> = Vec::new();
    let mut sorted_intervals = intervals.clone();
    sorted_intervals.sort_by(|a, b| a.0.cmp(&b.0));

    merged.push(sorted_intervals[0].clone());
    for current in sorted_intervals.iter().skip(1) {
        let prev = merged.last_mut().unwrap();
        if current.0 <= prev.1 {
            prev.1 = prev.1.max(current.1);
        } else {
            merged.push(current.clone());
        }
    }
    merged
}

/// Computes the query and subject coverage for a set of alignments to a single target.
///
/// This function calculates the total coverage of query and subject sequences by merging overlapping
/// or consecutive alignments and then calculating the coverage as a percentage of the total sequence lengths.
///
/// # Arguments
///
/// * `alignments` - A vector of references to `Alignment` objects representing alignments between a query
///   and a target sequence.
///
/// # Returns
///
/// Returns a tuple containing the query coverage (`qcov`) and subject coverage (`scov`) as percentages.
///
/// # Example
///
/// ```no_run
/// use vircov::{compute_coverage_single_target, Alignment};
///
/// // Assuming alignments is a Vec<&Alignment> populated with relevant data
/// let alignments: Vec<&Alignment> = vec![]; // Example placeholder
///
/// let (qcov, scov) = compute_coverage_single_target(alignments);
/// println!("Query coverage: {:.2}%, Subject coverage: {:.2}%", qcov, scov);
/// ```
pub fn compute_coverage_single_target(alignments: Vec<&Alignment>) -> (f64, f64) {
    let q_intervals: Vec<(usize, usize)> = alignments
        .iter()
        .map(|aln| (aln.qstart, aln.qend))
        .collect();
    let s_intervals: Vec<(usize, usize)> = alignments
        .iter()
        .map(|aln| (aln.sstart, aln.send))
        .collect();

    let q_merged = merge_intervals_single_target(q_intervals);
    let s_merged = merge_intervals_single_target(s_intervals);

    let q_coverage: usize = q_merged.iter().map(|(start, end)| end - start + 1).sum();
    let s_coverage: usize = s_merged.iter().map(|(start, end)| end - start + 1).sum();

    let qlen: usize = alignments.first().unwrap().qlen;
    let slen: usize = alignments.first().unwrap().slen;

    let qcov = q_coverage as f64 / qlen as f64 * 100.0;
    let scov = s_coverage as f64 / slen as f64 * 100.0;

    (qcov, scov)
}

/// Merges overlapping or consecutive intervals into the smallest number of non-overlapping intervals.
///
/// This function sorts the intervals by their start positions and then merges any intervals
/// that overlap or are consecutive into a single interval. This is useful for simplifying
/// the representation of ranges, such as genomic alignments or time intervals.
///
/// # Arguments
///
/// * `intervals` - A mutable vector of tuples `(usize, usize)` representing the intervals to be merged,
///   where the first element of each tuple is the start and the second is the end.
///
/// # Returns
///
/// Returns a vector of tuples representing the merged intervals.
///
/// # Example
///
/// ```
/// use vircov::merge_intervals;
///
/// let intervals = vec![(1, 3), (2, 4), (6, 8), (7, 10)];
/// let merged_intervals = merge_intervals(intervals);
/// assert_eq!(merged_intervals, vec![(1, 4), (6, 10)]);
/// ```
pub fn merge_intervals(mut intervals: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
    if intervals.is_empty() {
        return vec![];
    }

    intervals.sort_by(|a, b| a.0.cmp(&b.0));
    let mut merged: Vec<(usize, usize)> = vec![intervals[0]];

    for interval in intervals.iter().skip(1) {
        let last = merged.last_mut().unwrap();
        if interval.0 <= last.1 {
            last.1 = std::cmp::max(last.1, interval.1);
        } else {
            merged.push(*interval);
        }
    }
    merged
}

/// Computes the coverage of query and subject sequences from a set of alignments.
///
/// This function calculates the coverage for query and subject sequences based on their alignments.
/// Coverage is computed by first organizing the start and end positions of alignments by sequence ID,
/// then merging overlapping or consecutive intervals to calculate the total coverage. Coverage is
/// expressed as a percentage of the total length of all query or subject sequences.
///
/// # Arguments
///
/// * `alignments` - A slice of `Alignment` structs representing the alignments from which to compute coverage.
///
/// # Returns
///
/// Returns a tuple containing the query coverage (qcov) and subject coverage (scov) as percentages.
///
/// # Example
///
/// ```
/// use vircov::{compute_coverage, Alignment};
///
/// // Assuming alignments is an array or Vec<Alignment> populated with relevant data
/// let alignments: Vec<Alignment> = vec![]; // Example placeholder
///
/// let (qcov, scov) = compute_coverage(&alignments);
/// println!("Query coverage: {:.2}%, Subject coverage: {:.2}%", qcov, scov);
/// ```
pub fn compute_coverage(alignments: &[Alignment]) -> (f64, f64) {
    let mut q_intervals_by_seqid: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
    let mut s_intervals_by_seqid: HashMap<String, Vec<(usize, usize)>> = HashMap::new();

    // Organize intervals by qseqid and sseqid
    for aln in alignments {
        q_intervals_by_seqid
            .entry(aln.qseqid.clone())
            .or_default()
            .push((aln.qstart, aln.qend));
        s_intervals_by_seqid
            .entry(aln.sseqid.clone())
            .or_default()
            .push((aln.sstart, aln.send));
    }

    // Merge intervals and calculate total coverage
    let mut total_q_coverage: usize = 0;
    let mut total_s_coverage: usize = 0;

    for (_, intervals) in q_intervals_by_seqid {
        let merged = merge_intervals(intervals);
        total_q_coverage += merged
            .iter()
            .map(|(start, end)| end - start + 1)
            .sum::<usize>();
    }

    for (_, intervals) in s_intervals_by_seqid {
        let merged = merge_intervals(intervals);
        total_s_coverage += merged
            .iter()
            .map(|(start, end)| end - start + 1)
            .sum::<usize>();
    }

    let total_qlen: usize = alignments.iter().map(|aln| aln.qlen).sum();
    let total_slen: usize = alignments.iter().map(|aln| aln.slen).sum();

    let qcov = total_q_coverage as f64 / total_qlen as f64 * 100.0;
    let scov = total_s_coverage as f64 / total_slen as f64 * 100.0;

    (qcov, scov)
}

/// Parses BLASTN output lines into an iterator of `BlastRecord` results.
///
/// This function takes an iterator over `String`s, each representing a line of BLASTN output,
/// and attempts to parse each line into a `BlastRecord`. It handles parsing errors by returning
/// `Result<BlastRecord, SubtypeDatabaseError>`.
///
/// # Arguments
///
/// * `handle` - An iterator over `String` instances, each a line from a BLASTN output file.
///
/// # Returns
///
/// An iterator over `Result<BlastRecord, SubtypeDatabaseError>`, where each `Ok(BlastRecord)` represents
/// a successfully parsed record, and each `Err(SubtypeDatabaseError)` represents a failed parsing attempt.
///
/// # Examples
///
/// ```
/// let blast_output = vec![
///     "query1 target1 99.9 1000 1 1000 1 1000 1000 1000 1e-5".to_string(),
/// ];
/// let parsed = parse_blast(blast_output.into_iter()).collect::<Result<Vec<_>, _>>();
/// assert!(parsed.is_ok());
/// ```

pub fn parse_blast<I>(handle: I) -> impl Iterator<Item = Result<BlastRecord, SubtypeDatabaseError>>
where
    I: Iterator<Item = String>,
{
    handle.map(|line| {
        let r: Vec<&str> = line.split_whitespace().collect();
        if r.len() < 10 {
            return Err(SubtypeDatabaseError::LineParse(line));
        }
        let qcoords = {
            let mut coords = [r[6].parse::<i64>()?, r[7].parse::<i64>()?];
            coords.sort_unstable();
            (coords[0], coords[1])
        };
        let tcoords = {
            let mut coords = [r[8].parse::<i64>()?, r[9].parse::<i64>()?];
            coords.sort_unstable();
            (coords[0], coords[1])
        };
        Ok(BlastRecord {
            qname: r[0].to_string(),
            tname: r[1].to_string(),
            pid: r[2]
                .parse()
                .map_err(|_| SubtypeDatabaseError::LineParse(line.clone()))?,
            len: r[3]
                .parse()
                .map_err(|_| SubtypeDatabaseError::LineParse(line.clone()))?,
            qcoords,
            tcoords,
            qlen: r[r.len() - 2]
                .parse()
                .map_err(|_| SubtypeDatabaseError::LineParse(line.clone()))?,
            tlen: r[r.len() - 1]
                .parse()
                .map_err(|_| SubtypeDatabaseError::LineParse(line.clone()))?,
            evalue: r[r.len() - 4]
                .parse()
                .map_err(|_| SubtypeDatabaseError::LineParse(line.clone()))?,
        })
    })
}

/// Groups alignments by query and target name, excluding self-hits.
///
/// This function processes an iterator over `Result<BlastRecord, SubtypeDatabaseError>` and groups consecutive
/// alignments by their query and target names. It skips over self-hits (where query and target names are the same)
/// and accumulates the alignments into blocks, each block containing alignments between the same query and target.
///
/// # Arguments
///
/// * `handle` - An iterator over `Result<BlastRecord, SubtypeDatabaseError>`, typically produced by `parse_blast`.
///
/// # Returns
///
/// A `Result` containing either a `Vec<Vec<BlastRecord>>` representing groups of alignments between the same query and target,
/// or a `SubtypeDatabaseError` if an error occurs during processing.
///
/// # Examples
///
pub fn collect_alignment_blocks(
    mut handle: impl Iterator<Item = Result<BlastRecord, SubtypeDatabaseError>>,
) -> Result<Vec<Vec<BlastRecord>>, SubtypeDatabaseError> {
    let mut blocks: Vec<Vec<BlastRecord>> = Vec::new();
    let mut current_block: Vec<BlastRecord> = Vec::new();

    while let Some(result) = handle.next() {
        let aln = result?;
        if aln.qname == aln.tname {
            continue; // Skip self-hits as before
        }

        // If we're starting a new block or still within the same block
        if current_block.is_empty()
            || (current_block[0].qname == aln.qname && current_block[0].tname == aln.tname)
        {
            current_block.push(aln);
        } else {
            // We've encountered the start of a new block, so save the old one and start fresh
            blocks.push(current_block);
            current_block = vec![aln]; // Start new block with current alignment
        }
    }

    // Add the last block if it's not empty
    if !current_block.is_empty() {
        blocks.push(current_block);
    }

    Ok(blocks)
}

/// Filters alignments based on length and e-value thresholds.
///
/// Given a vector of `BlastRecord`s, this function returns a new vector containing
/// only those records that meet the specified minimum length and maximum e-value criteria.
///
/// # Arguments
///
/// * `alns` - A vector of `BlastRecord` instances to filter.
/// * `min_length` - The minimum alignment length to include.
/// * `min_evalue` - The maximum e-value to include.
///
/// # Returns
///
/// A vector of `BlastRecord` instances that meet the specified criteria.
///
/// # Examples
///
/// ```
/// let alns = vec![
///     BlastRecord {
///         qname: "query1".to_string(),
///         tname: "target1".to_string(),
///         pid: 99.9,
///         len: 1000.0,
///         qcoords: (1, 1000),
///         tcoords: (1, 1000),
///         qlen: 1000.0,
///         tlen: 1000.0,
///         evalue: 1e-5,
///     },
///     BlastRecord {
///         qname: "query1".to_string(),
///         tname: "target2".to_string(),
///         pid: 95.0,
///         len: 800.0,
///         qcoords: (1, 800),
///         tcoords: (1, 800),
///         qlen: 1000.0,
///         tlen: 850.0,
///         evalue: 2e-4,
///     },
/// ];
/// let filtered_alns = prune_alns(alns, 900.0, 1e-3);
/// assert_eq!(filtered_alns.len(), 1);
/// assert_eq!(filtered_alns[0].tname, "target1");
/// ```

pub fn prune_alns(alns: Vec<BlastRecord>, min_length: f64, max_evalue: f64) -> Vec<BlastRecord> {
    alns.into_iter()
        .filter(|aln| aln.len >= min_length && aln.evalue <= max_evalue)
        .collect()
}

/// Computes the Average Nucleotide Identity (ANI) from a slice of alignments.
///
/// The ANI is calculated as the weighted average of the percent identity (`pid`)
/// of each alignment, weighted by the alignment length.
///
/// # Arguments
///
/// * `alns` - A slice of `BlastRecord` instances for which to compute the ANI.
/// * `round` - A boolean indicating whether to round the resuling coverage values.
///
/// # Returns
///
/// The ANI as a floating-point number.
///
/// # Examples
///
/// ```
/// let alns = vec![
///     BlastRecord {
///         qname: "query1".to_string(),
///         tname: "target1".to_string(),
///         pid: 99.9,
///         len: 1000.0,
///         qcoords: (1, 1000),
///         tcoords: (1, 1000),
///         qlen: 1000.0,
///         tlen: 1000.0,
///         evalue: 1e-5,
///     },
///     BlastRecord {
///         qname: "query1".to_string(),
///         tname: "target2".to_string(),
///         pid: 98.0,
///         len: 800.0,
///         qcoords: (1, 800),
///         tcoords: (1, 800),
///         qlen: 1000.0,
///         tlen: 850.0,
///         evalue: 2e-4,
///     },
/// ];
/// let ani = compute_ani(&alns, true);
/// assert_eq!(ani, 99.24);
/// ```
pub fn compute_ani(alns: &[BlastRecord], round: bool) -> f64 {
    let total_pid_len: f64 = alns.iter().map(|a| a.len * a.pid).sum();
    let total_len: f64 = alns.iter().map(|a| a.len).sum();
    if round {
        (total_pid_len / total_len * 100.0).round() / 100.0
    } else {
        (total_pid_len / total_len * 100.0) / 100.0
    }
}

/// Calculates the query and target coverage percentages from a slice of alignments.
///
/// Coverage is computed by merging overlapping alignment coordinates and calculating
/// the total covered length as a percentage of the query and target lengths, respectively.
///
/// # Arguments
///
/// * `alns` - A slice of `BlastRecord` instances for which to compute coverage.
/// * `round` - A boolean indicating whether to round the resuling coverage values.
///
/// # Returns
///
/// A tuple containing the query coverage and target coverage percentages, respectively.
///
/// # Examples
///
/// ```
/// let alns = vec![
///     BlastRecord {
///         qname: "query1".to_string(),
///         tname: "target1".to_string(),
///         pid: 99.9,
///         len: 1000.0,
///         qcoords: (1, 1000),
///         tcoords: (1, 1000),
///         qlen: 1200.0,
///         tlen: 1000.0,
///         evalue: 1e-5,
///     },
///     BlastRecord {
///         qname: "query1".to_string(),
///         tname: "target2".to_string(),
///         pid: 98.0,
///         len: 200.0,
///         qcoords: (801, 1000),
///         tcoords: (601, 800),
///         qlen: 1200.0,
///         tlen: 850.0,
///         evalue: 2e-4,
///     },
/// ];
/// let (qcov, tcov) = compute_cov(&alns, true);
/// assert_eq!(qcov, 100.0);
/// assert_eq!(tcov, 84.0);
/// ```
pub fn compute_cov(alns: &[BlastRecord], round: bool) -> (f64, f64) {
    let merge_coords = |coords: &[(i64, i64)]| -> Vec<(i64, i64)> {
        let mut nr_coords: Vec<(i64, i64)> = vec![];
        for &(start, stop) in coords {
            if let Some(last) = nr_coords.last_mut() {
                if start <= last.1 + 1 {
                    last.1 = last.1.max(stop);
                    continue;
                }
            }
            nr_coords.push((start, stop));
        }
        nr_coords
    };

    let qcov = {
        let coords: Vec<_> = alns.iter().map(|a| a.qcoords).collect();
        let nr_coords = merge_coords(&coords);
        100.0 * nr_coords.iter().map(|&(s, e)| e - s + 1).sum::<i64>() as f64 / alns[0].qlen
    };

    let tcov = {
        let coords: Vec<_> = alns.iter().map(|a| a.tcoords).collect();
        let nr_coords = merge_coords(&coords);
        100.0 * nr_coords.iter().map(|&(s, e)| e - s + 1).sum::<i64>() as f64 / alns[0].tlen
    };

    if round {
        (qcov.round() / 100.0, tcov.round() / 100.0)
    } else {
        (qcov / 100.0, tcov / 100.0)
    }
}

/*
================
Helper functions
================
*/

// Run a command that has an output arg
pub fn run_command(args: Vec<String>, program: &str) -> Result<(), SubtypeDatabaseError> {
    let output = Command::new(program)
        .args(args)
        .output()
        .map_err(|_| SubtypeDatabaseError::ProgramExecutionFailed(program.to_string()))?;

    // Ensure command ran successfully
    if !output.status.success() {
        return Err(SubtypeDatabaseError::CommandExecutionFailed(
            String::from_utf8_lossy(&output.stderr).to_string(),
        ));
    }
    Ok(())
}

#[derive(Error, Debug)]
pub enum DecompressionError {
    #[error("I/O error")]
    Io(#[from] std::io::Error),
    #[error("Decompression error")]
    Decompression(#[from] niffler::Error),
}

/// Decompresses a `.tar.gz` or `.tar.xz` archive into a specified output directory.
///
/// The function automatically detects the compression format of the archive (either gzip or xz)
/// using the `niffler` crate and extracts its contents to the given output directory.
///
/// # Arguments
///
/// * `archive_path` - A `PathBuf` pointing to the archive file to be decompressed.
/// * `output_dir` - A `PathBuf` indicating the directory where the archive's contents will be extracted.
///
/// # Errors
///
/// Returns a `Result<(), DecompressionError>` indicating the operation's success or failure.
/// Errors can occur due to issues with file I/O or during the decompression process.
///
/// # Examples
///
/// Assuming you have a `.tar.gz` or `.tar.xz` archive at `"/path/to/archive.tar.gz"` and
/// you want to extract it to `"/path/to/output_dir"`, you can use the function as follows:
///
/// ```
/// use std::path::PathBuf;
///
/// fn main() -> Result<(), Box<dyn std::error::Error>> {
///     let archive_path = PathBuf::from("/path/to/archive.tar.gz");
///     let output_dir = PathBuf::from("/path/to/output_dir");
///
///     // Attempt to decompress the archive
///     match decompress_archive(&archive_path, &output_dir) {
///         Ok(_) => println!("Decompression successful."),
///         Err(e) => println!("Failed to decompress: {}", e),
///     }
///
///     Ok(())
/// }
/// ```
pub fn decompress_archive(
    archive_path: &PathBuf,
    output_dir: &PathBuf,
) -> Result<(), DecompressionError> {
    let archive_file = File::open(archive_path)?;
    let buffered_reader = BufReader::new(archive_file);
    let (decompressor, _compression_format) = niffler::get_reader(Box::new(buffered_reader))?;
    let mut archive = Archive::new(decompressor);
    archive.unpack(output_dir)?;
    Ok(())
}


#[derive(Error, Debug)]
pub enum FileNameError {
    #[error("The path does not have a file name")]
    NoFileName,
    #[error("The file name is not valid Unicode")]
    InvalidUnicode,
}

/// Extracts the file name without any extensions from a given `PathBuf`.
///
/// This function is designed to handle cases with multiple extensions (e.g., "archive.tar.gz")
/// and will strip all extensions to return just the base file name.
///
/// # Arguments
///
/// * `path` - A reference to the `PathBuf` from which to extract the file name.
///
/// # Returns
///
/// This function returns a `Result<String, FileNameError>`, where the `Ok` variant
/// contains the base file name as a `String`, and the `Err` variant contains an error
/// of type `FileNameError`.
///
/// # Examples
///
/// ```
/// use std::path::PathBuf;
/// use vircov::subtype::get_base_file_name;
///
/// fn main() -> Result<(), Box<dyn std::error::Error>> {
///     let path = PathBuf::from("/path/to/archive.tar.gz");
///     let file_name = get_base_file_name(&path)?;
///     assert_eq!(file_name, "archive");
///     Ok(())
/// }
/// ```
pub fn get_base_file_name(path: &PathBuf) -> Result<String, FileNameError> {
    let file_name = path.file_name().ok_or(FileNameError::NoFileName)?;
    let file_name_str = file_name.to_str().ok_or(FileNameError::InvalidUnicode)?;
    let base_name = file_name_str
        .split('.')
        .next()
        .ok_or(FileNameError::NoFileName)?
        .to_string();
    Ok(base_name)
}



/// Filters a list of `SubtypeSummary` based on specified minimum coverages.
///
/// The function applies filters for:
/// - `ani_tcov`: Minimum nucleotide target coverage.
/// - `aai_tcov_aa`: Minimum amino acid coverage.
/// - `aai_tcov_protein`: Minimum target protein coverage.
///
/// # Examples
///
/// ```
/// use vircov::filter_subtype_summaries;
/// let summaries = vec![vircov::SubtypeSummary {
///     query: "query1".to_string(),
///     target: "target1".to_string(),
///     target_proteins: Some(10),
///     shared_proteins: Some(5),
///     aai_tcov_protein: Some(0.9),
///     aai_tcov_aa: Some(0.8),
///     aai: Some(0.85),
///     aai_weighted: Some(0.90),
///     ani_qcov: Some(0.95),
///     ani_tcov: Some(0.97),
///     ani: Some(0.99),
///     ani_weighted: Some(0.98),
/// }];
/// let filtered = filter_subtype_summaries(&summaries, 0.90, 0.85, 0.80).unwrap();
/// assert_eq!(filtered.len(), 1);
/// ```
///
/// # Errors
///
/// This function returns `Err(FilterError::MissingData)` if any required coverage data is missing.
pub fn filter_subtype_summaries(
    summaries: &Vec<SubtypeSummary>,
    min_cov: f64,
    min_cov_aa: f64,
    min_cov_prot: f64,
) -> Result<Vec<SubtypeSummary>, SubtypeDatabaseError> {
    let filtered_summaries = summaries.iter().filter(|s| {
        s.ani_tcov.map_or(true, |tcov| tcov >= min_cov)
            && s.aai_tcov_aa.map_or(true, |tcov_aa| tcov_aa >= min_cov_aa)
            && s.aai_tcov_protein.map_or(true, |tcov_prot| tcov_prot >= min_cov_prot)
    })
    .cloned()
    .collect::<Vec<_>>();
    Ok(filtered_summaries)
}

#[derive(Debug, PartialEq, Eq, Hash)]
struct AnnotationRecord {
    accession: String,
    assembly: String
}


trait IntoPathBuf {
    fn into_path_buf(&self) -> PathBuf;
}

// Implement IntoFormattedString for f64 to format as percentage
impl IntoPathBuf for &str {
    fn into_path_buf(&self) -> PathBuf {
        PathBuf::from(self)
    }
}

// Implement IntoFormattedString for PathBuf to format as percentage
impl IntoPathBuf for PathBuf {
    fn into_path_buf(&self) -> PathBuf {
        self.to_path_buf()
    }
}

/// Parses a FASTA file and returns a vector of `Record`s.
///
/// # Arguments
///
/// * `file_path` - A string slice that holds the path to the file
///
/// # Examples
///
/// ```
/// let proteins = parse_fasta_file("proteins.fasta").unwrap();
/// let genomes = parse_fasta_file("genomes.fasta").unwrap();
/// ```
fn parse_fasta_annotations_file<T>(file_path: T) -> Result<Vec<AnnotationRecord>, SubtypeDatabaseError>
where T: AsRef<Path>
{
    let file = File::open(&file_path)?;
    let reader = BufReader::new(file);
    let mut fasta_reader = noodles::fasta::Reader::new(reader);

    let mut records = Vec::new();

    for result in fasta_reader.records() {
        let record = result?;
        let header = record.description();
        
        match header {
            Some(header) => {
                let record = AnnotationRecord {
                    accession: record.name().to_string(),
                    assembly: header.trim_start_matches("|").to_string(),
                };
                records.push(record);
            },
            None => {
                continue
            }
        }
       
    }
    Ok(records)
}

/// Matches protein accessions to nucleotide accessions from record descriptions
///
/// # Arguments
///
/// * `proteins` - A reference to a vector of protein `AnnotationRecord`s.
/// * `genomes` - A reference to a vector of genome `AnnotationRecord`s.
///
/// # Returns
///
/// A HashMap where each tuple contains a protein accession and a matching genome accession.
fn match_accessions(proteins: &[AnnotationRecord], genomes: &[AnnotationRecord]) -> HashMap<String, String> {
    let mut matches = HashMap::new();

    for protein in proteins {
        for genome in genomes {
            if protein.assembly == genome.assembly
            {
                matches.insert(protein.accession.clone(), genome.accession.clone());
            }
        }
    }

    matches
}

#[derive(Debug, Deserialize)]
struct NcbiGenotype {
    #[serde(rename(deserialize = "Accession"))]
    accession: String,
    #[serde(rename(deserialize = "Genotype"))]
    genotype: Option<String>,
    #[serde(rename(deserialize = "Segment"))]
    segment: Option<String>,
    #[serde(rename(deserialize = "Isolate"))]
    isolate: Option<String>,
    #[serde(rename(deserialize = "Organism_Name"))]
    organism_name: Option<String>,
    #[serde(rename(deserialize = "GenBank_Title"))]
    genbank_title: Option<String>,
}


pub fn process_ncbi_genotypes(
    ncbi_genotypes: &PathBuf,
    vircov_genotypes: &PathBuf,
    subtypes: &IndexMap<&str, Vec<&str>>
) -> Result<(), SubtypeDatabaseError> {

    let tsv = false;

    let file = File::open(ncbi_genotypes)?;
    let mut rdr = csv::ReaderBuilder::new()
        .flexible(true)
        .from_reader(file);

    let records = rdr.deserialize::<NcbiGenotype>();
    let mut genotypes = Vec::new();
    for record in records {
        let genotype = record.map_err(SubtypeDatabaseError::DeserializeRecord)?;
        
        let mut genotype_extracted = None;
        for (subtype, terms) in subtypes {
            let mut subtype_found = false;
            if *subtype == "regex" {
                for regex_pattern in terms {
                    let re = Regex::new(regex_pattern)?;
                    if let Some(caps) = re.captures(genotype.organism_name.as_deref().unwrap_or_default()) {
                        if let Some(mat) = caps.get(1) {
                            genotype_extracted = Some(mat.as_str().to_string());
                            subtype_found = true;
                            break; // Term found
                        }
                    }
                }
            } else {
                for term in terms {
                    if genotype.isolate.as_ref().map_or(false, |title| title.contains(term)) || 
                    genotype.organism_name.as_ref().map_or(false, |name| name.contains(term)) ||
                    genotype.genbank_title.as_ref().map_or(false, |title| title.contains(term))  {
                        log::info!("{} Genotype: {:?}, Subtype: {}", genotype.accession, genotype.genotype, subtype);
                        genotype_extracted = Some(subtype.to_string());
                        subtype_found = true;
                        break; // Term found
                    }
                }
            }
            if subtype_found {
                break
            }
        }
        genotypes.push(Genotype {
            accession: genotype.accession.clone(),
            genotype: genotype_extracted,
            segment: genotype.segment.clone()
        });
    }

    let file = File::create(vircov_genotypes)?;
    let mut writer = if tsv {
        csv::WriterBuilder::new().delimiter(b'\t').has_headers(true).from_writer(file)
    } else {
        csv::WriterBuilder::new().has_headers(true).from_writer(file)
    };


    for genotype in genotypes {
        writer.serialize(genotype).map_err(SubtypeDatabaseError::SerializeRecord)?;
    }
    writer.flush()?;


    Ok(())
}




/// Concatenates multiple Fasta files into a single file.
///
/// The function takes a base Fasta file and a list of Fasta files to append to the base file.
/// It writes the output to a new file specified by `output_path`.
///
/// # Arguments
///
/// * `base_file` - A `PathBuf` to the base Fasta file.
/// * `files_to_append` - A vector of `PathBuf` references to the Fasta files to append.
/// * `output_path` - A `PathBuf` to the output Fasta file.
///
/// # Returns
///
/// This function returns a `Result<(), ConcatError>`, which is `Ok` if the files were
/// successfully concatenated, or an `Err` with a `ConcatError` detailing what went wrong.
///
/// # Examples
///
/// ```
/// use std::path::PathBuf;
/// use vircov::subtype::concatenate_fasta_files;
///
/// let base_file = PathBuf::from("base.fasta");
/// let files_to_append = vec![PathBuf::from("append1.fasta"), PathBuf::from("append2.fasta")];
/// let output_path = PathBuf::from("output.fasta");
///
/// if let Err(e) = concatenate_fasta_files(base_file, &files_to_append, output_path) {
///     println!("An error occurred: {}", e);
/// }
/// ```
pub fn concatenate_fasta_files(base_file: &PathBuf, files_to_append: &[&PathBuf], output_path: &PathBuf) -> Result<(), SubtypeDatabaseError> {
    let mut output_file = File::create(&output_path)?;

    // Append base file content to the output file.
    let base_content = std::fs::read(&base_file)?;
    output_file.write_all(&base_content)?;

    // Iterate over files to append and write their content to the output file.
    for file_path in files_to_append {
        let content = std::fs::read(file_path)?;
        output_file.write_all(&content)?;
    }

    Ok(())
}


/// Writes a vector of `SubtypeSummary` to a CSV file with specified field names as column headers.
///
/// # Arguments
///
/// * `summaries` - A vector of `SubtypeSummary` structs to be written to the CSV.
/// * `file_path` - The path to the CSV file where the data will be written.
/// * `use_tsv` - Optional argument; if set to true, writes to a TSV file instead.
///
/// # Example
///
/// ```no_run
/// let summaries = vec![
///     SubtypeSummary {
///         query: "Query1".to_string(),
///         target: "Target1".to_string(),
///         target_proteins: Some(100),
///         shared_proteins: None,
///         aai_tcov_protein: Some(98.6),
///         aai_tcov_aa: None,
///         aai: Some(99.2),
///         aai_weighted: None,
///         ani_qcov: Some(95.5),
///         ani_tcov: None,
///         ani: Some(96.4),
///         ani_weighted: Some(97.0),
///     },
/// ];
/// write_subtype_summaries_to_csv(&summaries, "path/to/output.csv", false).unwrap();
/// ```
///
/// # Errors
///
/// Returns an error if writing to the file fails.
pub fn write_subtype_summaries_to_csv(
    summaries: &[SubtypeSummary],
    file_path: &PathBuf,
    tsv: bool,
) -> Result<(), SubtypeDatabaseError> {


    let file = File::create(file_path)?;
    let mut writer = if tsv {
        csv::WriterBuilder::new().delimiter(b'\t').has_headers(true).from_writer(file)
    } else {
        csv::WriterBuilder::new().has_headers(true).from_writer(file)
    };


    for summary in summaries {
        writer.serialize(summary).map_err(SubtypeDatabaseError::SerializeRecord)?;
    }
    writer.flush()?;
    Ok(())
}

#[derive(Debug, Deserialize)]
struct NextstrainClade {
    name: String,
    clade: String,
}


/// Parses a TSV file with NextstrainClade data and a FASTA file to associate clades with sequence IDs,
/// and outputs modified sequences and genotype records.
///
/// # Arguments
///
/// * `tsv_path` - Path to the TSV file containing the Nextstrain clade information.
/// * `fasta_path` - Path to the FASTA file containing sequence data.
/// * `output_fasta` - Path to output modified FASTA sequences.
/// * `output_csv` - Path to output matched genotype records.
/// * `segment` - Optional parameter specifying the segment of the sequence.
///
/// # Errors
///
/// This function can return errors related to file reading and data processing as described by `SubtypeDatabaseError`.
///
/// # Examples
///
/// ```
/// use std::path::PathBuf;
///
/// let result = process_gisaid_genotypes(
///     Path::new("clades.tsv"),
///     Path::new("sequences.fasta"),
///     PathBuf::from("output_sequences.fasta"),
///     PathBuf::from("output_genotypes.csv"),
///     Some("Segment 1")
/// );
/// match result {
///     Ok(()) => println!("Processing completed successfully."),
///     Err(e) => println!("Error: {}", e),
/// }
/// ```
pub fn process_gisaid_genotypes(
    tsv_path: &Path,
    fasta_path: &Path,
    output_fasta: PathBuf,
    output_csv: PathBuf,
    segment: Option<&str>,
) -> Result<(), SubtypeDatabaseError> {

    // Parse the TSV file to get clade information
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(tsv_path)?;

    let clades: Vec<NextstrainClade> = reader.deserialize().collect::<Result<_, _>>()?;

    // Prepare to read the FASTA file and write matched sequences
    let mut fasta_reader = noodles::fasta::Reader::new(BufReader::new(File::open(fasta_path)?));
    let fasta_writer = BufWriter::new(File::create(output_fasta)?);
    let mut fasta_writer = noodles::fasta::Writer::new(fasta_writer);

    // Prepare to write matched Genotypes to CSV
    let csv_writer = BufWriter::new(File::create(output_csv)?);
    let mut csv_writer = csv::WriterBuilder::new().from_writer(csv_writer);

    // Process each FASTA record
    for result in fasta_reader.records() {
        let record = result?;
        let id = record.name().to_string();

        // Check if any clade name is contained within the sequence ID
        for clade in &clades {
            if id.contains(&clade.name) {
                let modified_record = Record::new(
                    Definition::new(clade.name.clone(), None), record.sequence().clone()
                );

                let genotype = Genotype {
                    accession: clade.name.clone(),
                    genotype: Some(clade.clade.clone()),
                    segment: segment.map(|s| s.to_string()),
                };

                // Write modified sequence to FASTA file
                fasta_writer.write_record(&modified_record)?;

                // Write matched genotype to CSV file
                csv_writer.serialize(&genotype)?;
                break; // Stop searching after the first match
            }
        }
    }

    Ok(())
}

/// Synchronizes entries in genotype CSV and sequence FASTA, applying filters and ensuring matched outputs.
///
/// # Arguments
///
/// * `genotype_csv_path` - Path to the CSV file containing genotype records.
/// * `fasta_path` - Path to the FASTA file containing sequence data.
/// * `output_csv_path` - Path for the output filtered CSV file.
/// * `output_fasta_path` - Path for the output filtered FASTA file.
/// * `min_length` - Minimum sequence length required to retain a record.
/// * `remove_duplicates` - Whether to remove duplicate accessions.
/// * `entries_to_remove` - List of accession entries to remove from consideration.
///
/// # Errors
///
/// This function can return various errors related to file I/O and data parsing as encapsulated by `SubtypeDatabaseError`.
///
/// # Examples
///
/// ```
/// use std::path::PathBuf;
///
/// let result = filter_database(
///     PathBuf::from("genotypes.csv"),
///     PathBuf::from("sequences.fasta"),
///     PathBuf::from("synchronized_genotypes.csv"),
///     PathBuf::from("synchronized_sequences.fasta"),
///     100, // Minimum sequence length
///     true, // Remove duplicates
///     vec!["accession1".to_string(), "accession2".to_string()] // Entries to remove
/// );
///
/// match result {
///     Ok(()) => println!("Entries synchronized successfully."),
///     Err(e) => println!("Error occurred: {}", e),
/// }
/// ```
pub fn filter_database(
    genotype_csv_path: PathBuf,
    fasta_path: PathBuf,
    output_csv_path: PathBuf,
    output_fasta_path: PathBuf,
    min_length: usize,
    remove_duplicates: bool,
    entries_to_remove: Vec<String>
) -> Result<(), SubtypeDatabaseError> {
    // Set up CSV reading
    let mut rdr = csv::ReaderBuilder::new().from_path(genotype_csv_path)?;
    let mut genotype_map: HashMap<String, Genotype> = HashMap::new();
    let mut seen_accessions = HashSet::new();

    // Filter genotype entries based on entries to remove and duplicates
    for result in rdr.deserialize() {
        let genotype: Genotype = result?;
        if !entries_to_remove.contains(&genotype.accession) &&
            (!remove_duplicates || seen_accessions.insert(genotype.accession.clone())) {
            genotype_map.insert(genotype.accession.clone(), genotype);
        }
    }

    // Set up FASTA reading and writing
    let mut fasta_reader = noodles::fasta::Reader::new(BufReader::new(File::open(&fasta_path)?));
    let fasta_writer = BufWriter::new(File::create(&output_fasta_path)?);
    let mut fasta_writer = noodles::fasta::Writer::new(fasta_writer);

    // Set up CSV writing
    let mut wtr = csv::WriterBuilder::new().from_path(output_csv_path)?;

    let mut matched_genotypes = Vec::new();

    // Process FASTA entries and synchronize with genotypes
    for result in fasta_reader.records() {
        let record = result?;
        if record.sequence().len() >= min_length &&
           !entries_to_remove.contains(&record.name().to_string()) &&
           genotype_map.contains_key(record.name()) {
            if let Some(genotype) = genotype_map.get(record.name()) {
                // Write synchronized sequence to FASTA
                fasta_writer.write_record(&record)?;

                // Collect matched genotype for later output to maintain order
                matched_genotypes.push(genotype.clone());
            }
        }
    }

    // Write synchronized genotypes to CSV
    for genotype in matched_genotypes {
        wtr.serialize(&genotype)?;
    }

    Ok(())
}



/// Validates that the sequence names in a FASTA file match the accessions in a corresponding CSV file.
///
/// # Arguments
///
/// * `csv_path` - Path to the CSV file containing genotype records.
/// * `fasta_path` - Path to the FASTA file containing sequence data.
///
/// # Errors
///
/// This function will return an error if there are mismatches between the accessions and sequence names,
/// or if any I/O or parsing errors occur.
///
/// # Examples
///
/// ```
/// use std::path::Path;
///
/// let result = validate_genotypes(
///     Path::new("genotypes.csv"),
///     Path::new("sequences.fasta")
/// );
///
/// match result {
///     Ok(_) => println!("Validation successful."),
///     Err(e) => println!("Validation error: {}", e),
/// }
/// ```
pub fn validate_genotypes(csv_path: &Path, fasta_path: &Path) -> Result<(), SubtypeDatabaseError> {
    // Set up CSV reader
    let mut rdr = csv::ReaderBuilder::new().has_headers(true).from_path(csv_path)?;
    
    // Set up FASTA reader
    let mut fasta_reader = noodles::fasta::Reader::new(BufReader::new(File::open(fasta_path)?));

    // Collect genotype data
    let genotypes: Vec<Genotype> = rdr.deserialize().collect::<Result<Vec<_>, _>>()?;
    
    // Initialize mismatch collection
    let mut mismatches = Vec::new();
    
    // Check for mismatches
    for (index, record) in fasta_reader.records().enumerate() {
        let record = record?;
        if index < genotypes.len() {
            let genotype = &genotypes[index];
            if record.name() != genotype.accession {
                let message = format!("Mismatch at index {}: CSV '{}', FASTA '{}'", index, genotype.accession, record.name());
                log::info!("{}", message);
                mismatches.push(message);
            }
        }
    }

    // If there are mismatches, raise an error
    if !mismatches.is_empty() {
        return Err(SubtypeDatabaseError::DataMismatchError("Accessions do not match sequence names.".to_string()));
    } else {
        log::info!("Sequence and genotype file have matching order and accessions 🥳")
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::{self, File};
    use std::io::{BufWriter, Write};
    use tempfile::tempdir;

    #[test]
    fn test_decompress_archive() -> Result<()> {
        let temp_dir = tempdir()?;
        let output_dir = temp_dir.path().join("output");
        fs::create_dir(&output_dir)?;

        let archive_path = temp_dir.path().join("archive.tar.gz");
        let extracted_file_path = output_dir.join("extracted_file.txt");
        let extracted_file_contents = "Hello, world!";

        // Create a dummy .tar.gz archive for testing using niffler for compression
        {
            let file = File::create(&archive_path)?;
            let writer = BufWriter::new(file);
            let compressor = niffler::get_writer(
                Box::new(writer),
                niffler::Format::Gzip,
                niffler::compression::Level::Nine,
            )?;
            let mut tar = tar::Builder::new(compressor);
            let mut file = File::create(&extracted_file_path)?;
            writeln!(file, "{}", extracted_file_contents)?;
            tar.append_path_with_name(&extracted_file_path, "extracted_file.txt")?;
            // Ensure all data is written to the tar builder
            tar.into_inner()?;
            // Going out of scope drops the compressor and flushes the write
        }

        decompress_archive(&archive_path, &output_dir).expect("Decompression failed");

        assert!(
            PathBuf::from(&extracted_file_path).exists(),
            "The extracted file does not exist"
        );
        let contents = fs::read_to_string(extracted_file_path)?;
        assert_eq!(
            contents, extracted_file_contents,
            "The contents of the extracted file do not match"
        );

        Ok(())
    }
}
