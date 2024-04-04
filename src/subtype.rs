use anyhow::Result;
use std::num::{ParseFloatError, ParseIntError};
use std::path::PathBuf;
use tabled::{Table, Tabled};
use thiserror::Error;
use std::collections::{HashMap, HashSet};

use std::process::{Command, Stdio};
use std::fs::{create_dir_all, File};
use std::io::{BufRead, BufReader};

use tar::Archive;


/*
========================
Custom error definitions
========================
*/

#[derive(Error, Debug)]
pub enum SubtypeError {
}

/*
====================
Output table display
====================
*/

fn display_coverage(cov: &f64) -> String {
    format!("{:.1}%", cov)
}


fn display_optional_coverage(cov: &Option<f64>) -> String {
    match cov {
        Some(cov) => format!("{:.1}%", cov),
        None => String::from("n/a")
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
    db_blastx: PathBuf
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
    #[error("Failed to parse line: {0}")]
    LineParse(String),
    #[error("Parse integer error for record")]
    ParseInteger(#[from] ParseIntError),
    #[error("Parse value error for record")]
    ParseValue(String),
    #[error("Parse float error for record")]
    ParseFloat(#[from] ParseFloatError),
    /// Represents a failure to excecute a program in the simulation pipeline
    #[error("Failed to execute '{0}' - is it installed?")]
    ProgramExecutionFailed(String),
    /// Represents a failure to excecute a process command in the simulation pipeline
    #[error("Failed to run command, output is:\n{0}")]
    CommandExecutionFailed(String),
}

pub struct AlignmentConfig {
    blast: BlastConfig,
    diamond: DiamondConfig
}

impl Default for AlignmentConfig {
    fn default() -> Self {
        Self {
            blast: BlastConfig::default(),
            diamond: DiamondConfig::default()
        }
    }
}

pub struct BlastConfig {
    min_percent_identity: f64,
    max_target_seqs: u64
}
impl Default for BlastConfig {
    fn default() -> Self {
        Self {
            min_percent_identity: 50.,
            max_target_seqs: 25000,
        }
    }
}
pub struct DiamondConfig {
    max_evalue: f64,
    max_target_seqs: u64,
    min_query_cover: f64,
    min_target_cover: f64
}
impl Default for DiamondConfig {
    fn default() -> Self {
        Self {
            max_evalue: 0.00001,
            max_target_seqs: 10000,
            min_query_cover: 50.,
            min_target_cover: 50.
        }
    }
}


pub struct SubtypeDatabase {
    pub name: String,
    pub path: PathBuf,
    pub files: SubtypeDatabaseFiles
}

impl SubtypeDatabase {
    pub fn from(archive: &PathBuf, outdir: &PathBuf, threads: u32) -> Result<Self, SubtypeDatabaseError>  {

        create_dir_all(&outdir)?;

        // Decompress database archive into subtype database output directory
        decompress_archive(archive, outdir)?;

        let db_nuc = outdir.join("db").join("db_nuc.fasta");
        let db_aa = outdir.join("db").join("db_aa.fasta");
        let db_csv = outdir.join("db").join("db.csv");

        if !db_nuc.exists() {
            return Err(SubtypeDatabaseError::DatabaseFileMissing(db_nuc.display().to_string()))
        }
        if !db_aa.exists() {
            return Err(SubtypeDatabaseError::DatabaseFileMissing(db_aa.display().to_string()))
        }
        // if !db_csv.exists() {
        //     return Err(SubtypeDatabaseError::DatabaseFileMissing(db_csv.display().to_string()))
        // }
        
        let name = get_base_file_name(&archive)?;

        let db_blastx = outdir.join("blastx");
        let db_blastn = outdir.join("blastn");

        SubtypeDatabase::makedb_blast_nuc(&db_nuc, &db_blastn)?;
        SubtypeDatabase::makedb_blast_prot(&db_aa, &db_blastx)?;
        
        // SubtypeDatabase::makedb_diamond(&db_aa, &db_dmd, threads)?;

        Ok(Self { 
            name,
            path: outdir.to_path_buf(),
            files: SubtypeDatabaseFiles { 
                csv: db_csv, 
                fasta_aa: db_aa, 
                fasta_nuc: db_nuc,
                db_blastn,
                db_blastx
            }
        })
    }
    pub fn subtype(&self, fasta: &PathBuf, outdir: &PathBuf, alignment_config: Option<AlignmentConfig>, threads: u32) -> Result<(), SubtypeDatabaseError> {

        let config = alignment_config.unwrap_or(AlignmentConfig::default());
        let input_name = get_base_file_name(&fasta)?;

        let blastn_output = outdir.join(
            format!("{}.{}.blastn.tsv", input_name, self.name)
        );
        let blastx_output = outdir.join(
            format!("{}.{}.blastx.tsv", input_name, self.name)
        );

        log::info!("Computing nucleotide and protein alignments with BLAST...");

        self.run_blastn(
            fasta, 
            &self.files.db_blastn, 
            &blastn_output, 
            config.blast.min_percent_identity, 
            config.blast.max_target_seqs, 
            threads
        )?;
        // self.run_blast(
        //     fasta, 
        //     &self.files.db_blastx, 
        //     &blastx_output, 
        //     config.blast.min_percent_identity, 
        //     config.blast.max_target_seqs, 
        //     threads,
        //     true
        // )?;

        self.run_blastx(
            fasta, 
            &self.files.db_blastx, 
            &blastx_output, 
            config.diamond.max_evalue, 
            config.diamond.max_target_seqs, 
            threads
        )?;

        // log::info!("Computing amino acid alignments with Diamond...");

        // self.run_diamond(
        //     fasta, 
        //     &self.files.db_dmd, 
        //     &diamond_output, 
        //     config.diamond.max_evalue, 
        //     config.diamond.max_target_seqs, 
        //     config.diamond.min_query_cover, 
        //     config.diamond.min_target_cover, 
        //     threads
        // )?;

        log::info!("Computing average nucleotide identity and coverage from alignments with BLAST...");
        self.compute_blastn_ani_cov(&blastn_output)?;

        log::info!("Computing average amino acid identity and coverage from alignments with BLAST...");
        self.compute_blastx_aai_cov(&blastx_output)?;

        // self.compute_diamond_aai_cov(&diamond_output)?;


        Ok(())
    }

    // ANI

    pub fn compute_blastn_ani_cov(&self, blastn_results: &PathBuf) -> Result<(), SubtypeDatabaseError> {
        let file = File::open(blastn_results)?;
        let reader = BufReader::new(file);

        // Convert lines of the file into an iterator of BlastRecord results
        let parsed_lines = parse_blast(reader.lines().filter_map(Result::ok));

        let alignment_blocks = collect_alignment_blocks(parsed_lines)?;

        // Process blocks of alignments
        for alns in alignment_blocks {
            // Filter alignments based on length and e-value
            let pruned_alns = prune_alns(alns, 100.0, 1e-3);
            if pruned_alns.is_empty() {
                continue;
            }

            // Compute ANI and coverage
            let ani = compute_ani(&pruned_alns, false);
            let (qcov, tcov) = compute_cov(&pruned_alns, false);

            // Output results
            println!(
                "{qname} {tname} {aln_count} {ani:.2} {qcov:.2}% {tcov:.2}%",
                qname = pruned_alns[0].qname,
                tname = pruned_alns[0].tname,
                aln_count = pruned_alns.len(),
                ani = ani,
                qcov = qcov * 100.0,
                tcov = tcov * 100.0,
            );
        }

        Ok(())
    }

    // AAI

    pub fn compute_blastx_aai_cov(&self, blastx_results: &PathBuf) -> Result<(), SubtypeDatabaseError> {

        let genome_info = parse_aa_fasta(&self.files.fasta_aa)?;
        let genome_alns = parse_and_organize_alignments(blastx_results, &genome_info)?;
        let mut summaries = compute_summaries(&genome_info, &genome_alns);

        summaries.sort_by(|a, b| b.aai.partial_cmp(&a.aai).expect("Should not be an optional field") );

        let table = Table::new(summaries);

        println!("{}",  &table);

        Ok(())
    }

    // Helper methods

    pub fn makedb_blast_nuc(db_fasta: &PathBuf, output: &PathBuf) -> Result<(), SubtypeDatabaseError> {
        let args = vec![
            "-in".to_string(),
            db_fasta.display().to_string(),
            "-out".to_string(),
            output.display().to_string(),
            "-dbtype".to_string(),
            "nucl".to_string(),
        ];
        log::info!("Running command: makeblastdb {}", &args.join(" "));
        run_command(args, "makeblastdb")
    }

    pub fn makedb_blast_prot(db_fasta: &PathBuf, output: &PathBuf) -> Result<(), SubtypeDatabaseError> {
        let args = vec![
            "-in".to_string(),
            db_fasta.display().to_string(),
            "-out".to_string(),
            output.display().to_string(),
            "-dbtype".to_string(),
            "prot".to_string(),
        ];
        log::info!("Running command: makeblastdb {}", &args.join(" "));
        run_command(args, "makeblastdb")
    }
    
    pub fn run_blastn(
        &self,
        input: &PathBuf,
        database: &PathBuf,
        output: &PathBuf,
        min_percent_identity: f64,
        max_target_seqs: u64,
        threads: u32
    ) -> Result<(), SubtypeDatabaseError> {
        let mut args = vec![
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
            min_percent_identity.to_string()
        ];

        log::info!("Running command: blastn with args {}", &args.join(" "));
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
            "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore".to_string(),
            "-evalue".to_string(),
            max_evalue.to_string(),
            "-max_target_seqs".to_string(),
            max_target_seqs.to_string(),
            "-num_threads".to_string(),
            threads.to_string(),
        ];
        log::info!("Running command: blastx {}", &args.join(" "));
        run_command(args, "blastx")
    }
}

// AAI

#[derive(Debug, Clone)]
struct GenomeInfo {
    protein_to_genome: HashMap<String, String>,
    num_proteins: HashMap<String, usize>,
}

#[derive(Debug, Clone, Tabled)]
struct AlignmentSummary {
    #[header("Query")]
    query: String,
    #[header("Target")]
    target: String,
    #[header("Target proteins")]
    target_proteins: usize,
    #[header("Shared proteins")]
    shared_proteins: usize,
    // #[header("Query protein coverage")]
    // #[field(display_with = "display_optional_coverage")]
    // protein_qcov: Option<f64>,
    #[header("Protein coverage")]
    #[field(display_with = "display_coverage")]
    protein_tcov: f64,
    #[header("AA coverage")]
    #[field(display_with = "display_coverage")]
    aa_tcov: f64,
    #[header("AAI")]
    #[field(display_with = "display_coverage")]
    aai: f64,
    #[header("AAIw")]
    #[field(display_with = "display_coverage")]
    weighted_aai: f64
}


#[derive(Debug)]
struct Alignment {
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

fn parse_blastx_output(file_path: &PathBuf, min_qcov: f64, min_tcov: f64) -> Result<Vec<Alignment>, SubtypeDatabaseError> {
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
            _length: parts[3].parse().map_err(SubtypeDatabaseError::ParseInteger)?,
            qlen: parts[4].parse().map_err(SubtypeDatabaseError::ParseInteger)?,
            slen: parts[5].parse().map_err(SubtypeDatabaseError::ParseInteger)?,
            qstart: parts[6].parse().map_err(SubtypeDatabaseError::ParseInteger)?,
            qend: parts[7].parse().map_err(SubtypeDatabaseError::ParseInteger)?,
            sstart: parts[8].parse().map_err(SubtypeDatabaseError::ParseInteger)?,
            send: parts[9].parse().map_err(SubtypeDatabaseError::ParseInteger)?,
            _evalue: parts[10].parse().map_err(SubtypeDatabaseError::ParseFloat)?,
            _bitscore: parts[11].parse().map_err(SubtypeDatabaseError::ParseFloat)?,
            qcov: 0.0, // To be calculated
            scov: 0.0, // To be calculated
        };

        // Calculate coverages
        alignment.qcov = (alignment.qend - alignment.qstart + 1) as f64 / alignment.qlen as f64 * 100.0;
        alignment.scov = (alignment.send - alignment.sstart + 1) as f64 / alignment.slen as f64 * 100.0;

        if alignment.qcov >= min_qcov && alignment.scov >= min_tcov {
            alignments.push(alignment);
        }
    }

    Ok(alignments)
}


fn parse_aa_fasta(file_path: &PathBuf) -> Result<GenomeInfo, SubtypeDatabaseError> {
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

fn parse_and_organize_alignments(
    file_path: &PathBuf, 
    info: &GenomeInfo
) -> Result<HashMap<String, HashMap<String, Vec<Alignment>>>, SubtypeDatabaseError> {

    let mut genome_alns: HashMap<String, HashMap<String, Vec<Alignment>>> = HashMap::new();
    let alignments = parse_blastx_output(file_path, 0., 0.)?;

    for aln in alignments {

        let genome = aln.qseqid.clone();
        let protein = aln.sseqid.clone();

        let target = info.protein_to_genome.get(&protein).ok_or(SubtypeDatabaseError::ParseValue(format!("Target protein not found: {}", &protein)))?.to_string();
        genome_alns.entry(genome)
            .or_insert_with(HashMap::new)
            .entry(target)
            .or_insert_with(Vec::new)
            .push(aln);
    }

    Ok(genome_alns)
}

fn compute_summaries(
    genome_info: &GenomeInfo,
    genome_alns: &HashMap<String, HashMap<String, Vec<Alignment>>>,
) -> Vec<AlignmentSummary> {
    let mut summaries = Vec::new();

    for (query, targets) in genome_alns {

        for (target, alignments) in targets {

            let total_target_proteins = *genome_info.num_proteins.get(target).unwrap_or(&0);

            let mut aln_by_gene = HashMap::new();
            for alignment in alignments {
                let ids: Vec<&str> = alignment.sseqid.split("__").collect();
                let gene_id = ids[1];

                aln_by_gene.entry(gene_id)
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
            let (_, aa_tcov) = compute_coverage(alignments);

            let aai = alignments.iter().map(|aln| aln.pident).sum::<f64>() / alignments.len() as f64;

            summaries.push(AlignmentSummary {
                query: query.clone(),
                target: target.split("::").last().unwrap().to_string(),
                target_proteins: total_target_proteins,
                shared_proteins,
                // protein_qcov: None, // full genome input for now
                protein_tcov,
                aa_tcov,
                aai,
                weighted_aai: aai*(protein_tcov/100.)*(aa_tcov/100.)
            });
        }
    }

    summaries
}


fn merge_intervals_single_target(intervals: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
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

fn compute_coverage_single_target(alignments: Vec<&Alignment>) -> (f64, f64) {
    let q_intervals: Vec<(usize, usize)> = alignments.iter().map(|aln| (aln.qstart, aln.qend)).collect();
    let s_intervals: Vec<(usize, usize)> = alignments.iter().map(|aln| (aln.sstart, aln.send)).collect();

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

fn merge_intervals(mut intervals: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
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

fn compute_coverage(alignments: &[Alignment]) -> (f64, f64) {
    let mut q_intervals_by_seqid: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
    let mut s_intervals_by_seqid: HashMap<String, Vec<(usize, usize)>> = HashMap::new();

    // Organize intervals by qseqid and sseqid
    for aln in alignments {
        q_intervals_by_seqid.entry(aln.qseqid.clone())
            .or_default()
            .push((aln.qstart, aln.qend));
        s_intervals_by_seqid.entry(aln.sseqid.clone())
            .or_default()
            .push((aln.sstart, aln.send));
    }

    // Merge intervals and calculate total coverage
    let mut total_q_coverage: usize = 0;
    let mut total_s_coverage: usize = 0;

    for (_, intervals) in q_intervals_by_seqid {
        let merged = merge_intervals(intervals);
        total_q_coverage += merged.iter().map(|(start, end)| end - start + 1).sum::<usize>();
    }

    for (_, intervals) in s_intervals_by_seqid {
        let merged = merge_intervals(intervals);
        total_s_coverage += merged.iter().map(|(start, end)| end - start + 1).sum::<usize>();
    }

    let total_qlen: usize = alignments.iter().map(|aln| aln.qlen).sum();
    let total_slen: usize = alignments.iter().map(|aln| aln.slen).sum();

    let qcov = total_q_coverage as f64 / total_qlen as f64 * 100.0;
    let scov = total_s_coverage as f64 / total_slen as f64 * 100.0;

    (qcov, scov)
}


#[derive(Debug)]
struct BlastRecord {
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

fn parse_blast<I>(handle: I) -> impl Iterator<Item = Result<BlastRecord, SubtypeDatabaseError>>
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
            pid: r[2].parse().map_err(|_| SubtypeDatabaseError::LineParse(line.clone()))?,
            len: r[3].parse().map_err(|_| SubtypeDatabaseError::LineParse(line.clone()))?,
            qcoords,
            tcoords,
            qlen: r[r.len() - 2].parse().map_err(|_| SubtypeDatabaseError::LineParse(line.clone()))?,
            tlen: r[r.len() - 1].parse().map_err(|_| SubtypeDatabaseError::LineParse(line.clone()))?,
            evalue: r[r.len() - 4].parse().map_err(|_| SubtypeDatabaseError::LineParse(line.clone()))?,
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
fn collect_alignment_blocks(
    mut handle: impl Iterator<Item = Result<BlastRecord, SubtypeDatabaseError>>,
) -> Result<Vec<Vec<BlastRecord>>, SubtypeDatabaseError> {
    let mut blocks: Vec<Vec<BlastRecord>> = Vec::new();
    let mut current_block: Vec<BlastRecord> = Vec::new();

    while let Some(result) = handle.next() {
        let aln = result?; // This propagates errors upwards
        if aln.qname == aln.tname {
            continue; // Skip self-hits as before
        }

        // If we're starting a new block or still within the same block
        if current_block.is_empty() || (current_block[0].qname == aln.qname && current_block[0].tname == aln.tname) {
            current_block.push(aln);
        } else {
            // We've encountered the start of a new block, so save the old one and start fresh
            blocks.push(current_block);
            current_block = vec![aln]; // Start new block with current alignment
        }
    }

    // Don't forget to add the last block if it's not empty
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

fn prune_alns(alns: Vec<BlastRecord>, min_length: f64, min_evalue: f64) -> Vec<BlastRecord> {
    alns.into_iter()
        .filter(|aln| aln.len >= min_length && aln.evalue <= min_evalue)
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
fn compute_ani(alns: &[BlastRecord], round: bool) -> f64 {
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
fn compute_cov(alns: &[BlastRecord], round: bool) -> (f64, f64) {
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
        return Err(SubtypeDatabaseError::CommandExecutionFailed(String::from_utf8_lossy(&output.stderr).to_string()));
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
pub fn decompress_archive(archive_path: &PathBuf, output_dir: &PathBuf) -> Result<(), DecompressionError> {
    let archive_file = File::open(archive_path)?;
    let buffered_reader = BufReader::new(archive_file);
    let (decompressor, _compression_format) = niffler::get_reader(Box::new(buffered_reader))?;
    let mut archive = Archive::new(decompressor);
    archive.unpack(output_dir)?;
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
                niffler::compression::Level::Nine)
            ?;
            let mut tar = tar::Builder::new(compressor);
            let mut file = File::create(&extracted_file_path)?;
            writeln!(file, "{}", extracted_file_contents)?;
            tar.append_path_with_name(&extracted_file_path, "extracted_file.txt")?;
            // Ensure all data is written to the tar builder
            tar.into_inner()?;
            // Going out of scope drops the compressor and flushes the write
            
        }

        decompress_archive(&archive_path, &output_dir).expect("Decompression failed");

        assert!(PathBuf::from(&extracted_file_path).exists(), "The extracted file does not exist");
        let contents = fs::read_to_string(extracted_file_path)?;
        assert_eq!(contents, extracted_file_contents, "The contents of the extracted file do not match");

        Ok(())
    }
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
