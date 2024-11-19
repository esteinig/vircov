use crate::consensus::VircovConsensus;
use crate::utils::FileComponent;
use crate::vircov::ConsensusConfig;
use crate::{utils::get_file_component, vircov::HaplotypeConfig};
use crate::error::VircovError;

use regex::Regex;
use serde::{Deserialize, Serialize};
use std::fs::{create_dir_all, remove_dir_all};
use std::{path::{Path, PathBuf}, process::{Command, Output}};

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, clap::ValueEnum)]
pub enum Haplotyper {
    Floria,
    Viloca
}

pub struct HaplotypeRecord {
    id: String,
    snps: usize,
    
}


#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct VircovHaplotype {
    outdir: PathBuf,
    config: HaplotypeConfig
}
impl VircovHaplotype {
    pub fn new(outdir: &PathBuf, config: HaplotypeConfig) -> Result<Self, VircovError> {
        Self::check_haplotype_dependencies(&config)?;

        Ok(Self { outdir: outdir.to_path_buf(), config })
    }
    pub fn haplotype(&self) -> Result<Vec<HaplotypeRecord>, VircovError> {
        let records = match self.config.haplotyper {
            Haplotyper::Floria => self.run_floria()?,
            Haplotyper::Viloca => self.run_viloca()?
        };
        Ok(records)
    }
    pub fn check_haplotype_dependencies(config: &HaplotypeConfig) -> Result<(), VircovError> {
        match config.haplotyper {
            Haplotyper::Floria => {
                run_command("freebayes --version").map_err(|_| VircovError::HaplotypeDependencyMissing(String::from("freebayes")))?;
                run_command("floria --version").map_err(|_| VircovError::HaplotypeDependencyMissing(String::from("floria")))?;
                run_command("floria-strainer --version").map_err(|_| VircovError::HaplotypeDependencyMissing(String::from("floria-strainer")))?;
            },
            Haplotyper::Viloca => {
                run_command("viloca --version").map_err(|_| VircovError::HaplotypeDependencyMissing(String::from("floria-strainer")))?;
            }
        };
        Ok(())
    }
    pub fn run_floria(&self) -> Result<Vec<HaplotypeRecord>, VircovError> {

        let vcf = self.config.alignment.with_extension("vcf");

        log::info!("File name");
        let bin_name = get_file_component(&self.config.alignment, FileComponent::FileStem)?;
        let bin_outdir = self.outdir.join(format!("{bin_name}"));

        log::info!("Freebayes");
        self.run_freebayes(&self.config.alignment, &vcf)?;
        log::info!("Floria");
        let strain_bam = self.phase_floria(&vcf, &bin_outdir, &bin_name)?;
        log::info!("Consensus");
        self.assemble_consensus_haplotypes(&strain_bam, &bin_outdir, &bin_name)?;

        Ok(Vec::new())
    }
    pub fn run_viloca(&self) -> Result<Vec<HaplotypeRecord>, VircovError> {

        let bin_name = get_file_component(&self.config.alignment, FileComponent::FileStem)?;
        let bin_outdir = self.outdir.join(format!("{bin_name}"));

        self.phase_viloca(&bin_outdir, &bin_name)?;
        
        Ok(Vec::new())
    }

    pub fn phase_viloca(&self, outdir: &Path, name: &str) -> Result<Vec<PathBuf>, VircovError> {

        let cmd = format!(
            "viloca run -b {bam} -f {fasta} --threads {threads} --mode use_quality_scores",
            bam=self.config.alignment.display(),
            fasta=self.config.fasta.display(),
            threads=self.config.threads,
        );

        match run_command_dir(&cmd, outdir) {
            Ok(_) => {
                log::warn!("Viloca processed for bin '{name}'");
            },
            Err(_) => {

                log::error!("No strain haplotypes could be extracted for bin '{name}'");
                remove_dir_all(outdir)?;
                return Ok(Vec::new())
            }
        }

        Ok(Vec::new())

    }
    pub fn phase_floria(&self, vcf: &Path, outdir: &PathBuf, name: &str) -> Result<Vec<PathBuf>, VircovError> {
        
        let cmd = format!(
            "floria -b {bam} -r {fasta} -v {vcf} -o {outdir} --snp-count-filter {snp_count_filter} --snp-density {snp_density} --ploidy-sensitivity {ploidy_sensitivity} --threads {threads} --overwrite",
            bam=self.config.alignment.display(),
            fasta=self.config.fasta.display(),
            vcf=vcf.display(),
            outdir=outdir.display(),
            snp_count_filter=self.config.snp_count_filter,
            snp_density=self.config.snp_density,
            ploidy_sensitivity=self.config.ploidy_sensitivity,
            threads=self.config.threads
        );

        log::info!("Calling haplotypes for bin '{}'", name);

        match run_command(&cmd) {
            Ok(_) => {
                let basename = outdir.join(name);
                let cmd = format!(
                    "floria-strainer -b {bam} --mode split --basename {basename} {outdir}",
                    bam=self.config.alignment.display(),
                    basename=basename.display(),
                    outdir=outdir.display()
                );
                match run_command(&cmd) {
                    Ok(_) => {
                        let strain_bams = self.get_bam_files_sorted(outdir, name)?;
                        log::info!("Recovered {} strain haplotypes for bin '{}'", strain_bams.len(), name);
                        Ok(strain_bams)
                    },
                    Err(_) => {
                        log::info!("No strain haplotypes could be extracted for bin '{name}'");
                        remove_dir_all(outdir)?;
                        return Ok(Vec::new())
                    }
                }
            },
            Err(_) => {
                log::info!("No haplotypes could be called for bin '{name}'");
                log::warn!("{}", cmd);
                if outdir.exists() { remove_dir_all(outdir)?; }
                return Ok(Vec::new())
            }
        }
    }
    pub fn assemble_consensus_haplotypes(&self, bam: &Vec<PathBuf>, outdir: &PathBuf, name: &str) -> Result<(), VircovError> {

        for strain in bam {
            let strain_id = get_file_component(strain, FileComponent::FileStem)?;
            let strain_assembly = format!("{strain_id}.consensus.fa");

            log::info!("Assembling strain consensus genome for '{strain_id}'");
            let strain_config = ConsensusConfig::with_default(
                strain, 
                &outdir.join(strain_assembly), 
                self.config.reference.clone(), 
                Some(format!("{} bin={}", strain_id, name)), 
                self.config.min_consensus_quality, 
                self.config.min_consensus_freq,
                self.config.min_consensus_depth, 
                1000000
            );

            let consensus = VircovConsensus::new(strain_config)?;
            consensus.assemble()?;
        }
        Ok(())
    } 
    pub fn run_freebayes(&self, bam: &Path, vcf: &Path) -> Result<(), VircovError> {

        let cmd = format!(
            "freebayes -f {fasta} --pooled-continuous {bam} > {vcf}",
            fasta=self.config.fasta.display(),
            bam=bam.display(),
            vcf=vcf.display()
        );

        run_command(&cmd)?;

        Ok(())
    }
    fn get_bam_files_sorted(&self, outdir: &Path, prefix: &str) -> Result<Vec<PathBuf>, VircovError> {
        let re = Regex::new(&format!(r"^{}.(\d+).bam$", regex::escape(prefix)))?;
        
        let mut bam_files: Vec<(u32, PathBuf)> = std::fs::read_dir(outdir)
            .expect("Failed to read directory")
            .filter_map(|entry| entry.ok())
            .filter_map(|entry| {
                let path = entry.path();
                if let Some(filename) = path.file_name().and_then(|f| f.to_str()) {
                    if let Some(captures) = re.captures(filename) {
                        if let Ok(n) = captures[1].parse::<u32>() {
                            return Some((n, path));
                        }
                    }
                }
                None
            })
            .collect();
        
        bam_files.sort_by_key(|&(n, _)| n);

        Ok(bam_files.into_iter().map(|(_, path)| path).collect())
    }
}


fn run_command(command: &str) -> Result<Output, VircovError> {
    let output = Command::new("sh")
        .arg("-c")
        .arg(command)
        .output()
        .map_err(|e| VircovError::CommandExecutionFailed(command.to_string(), e.to_string()))?;

    if !output.status.success() {
        return Err(VircovError::CommandFailed(command.to_string(), output.status.code().unwrap_or(-1)));
    }

    Ok(output)
}


fn run_command_dir(command: &str, dir: &Path) -> Result<Output, VircovError> {

    if !dir.exists() {
        create_dir_all(&dir)?;
    }

    let output = Command::new("sh")
        .arg("-c")
        .arg(command)
        .current_dir(dir)
        .output()
        .map_err(|e| VircovError::CommandExecutionFailed(command.to_string(), e.to_string()))?;

    if !output.status.success() {
        return Err(VircovError::CommandFailed(command.to_string(), output.status.code().unwrap_or(-1)));
    }

    Ok(output)
}