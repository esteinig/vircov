use crate::vircov::HaplotypeConfig;
use crate::error::VircovError;

use serde::{Deserialize, Serialize};
use std::{path::{Path, PathBuf}, process::{Command, Output}};

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, clap::ValueEnum)]
pub enum HaplotypeVariantCaller {
    Freebayes
}

pub struct HaplotypeRecord {
    id: String,
    snps: usize,
    
}


#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct VircovHaplotype {
    config: HaplotypeConfig
}
impl VircovHaplotype {
    pub fn new(config: HaplotypeConfig) -> Result<Self, VircovError> {
        Self::check_haplotype_dependencies()?;

        Ok(Self { config })
    }
    pub fn haplotype(&self) -> Result<Vec<HaplotypeRecord>, VircovError> {
        let records = match self.config.variant_caller {
            HaplotypeVariantCaller::Freebayes => self.run_haplotype_recovery()?,
        };
        Ok(records)
    }
    pub fn check_haplotype_dependencies() -> Result<(), VircovError> {
        run_command("freebayes --version").map_err(|_| VircovError::HaplotypeDependencyMissing(String::from("freebayes")))?;
        run_command("floria --version").map_err(|_| VircovError::HaplotypeDependencyMissing(String::from("floria")))?;
        Ok(())
    }
    fn run_haplotype_recovery(&self) -> Result<Vec<HaplotypeRecord>, VircovError> {

        // Call variants

        let vcf = self.config.alignment.with_extension("vcf");
        let tsv = self.config.alignment.with_extension("tsv");

        self.call_variants(&self.config.alignment, &vcf)?;
        let haplotype_fastq = self.phase_strain_haplotypes(&vcf, &tsv)?;
        self.assemble_consensus_haplotypes(&vcf, &haplotype_fastq)?;

        Ok(Vec::new())
    }

    pub fn phase_strain_haplotypes(&self, vcf: &Path, tsv: &Path) -> Result<Vec<PathBuf>, VircovError> {

        let cmd = match self.config.variant_caller {
            HaplotypeVariantCaller::Freebayes => {
                format!(
                    "floria {vcf}",
                    vcf=vcf.display()
                )
            }
        };

        run_command(&cmd)?;

        Ok(Vec::new())
    }
    pub fn assemble_consensus_haplotypes(&self, vcf: &Path, fastq: &Vec<PathBuf>) -> Result<(), VircovError> {

        let cmd = match self.config.variant_caller {
            HaplotypeVariantCaller::Freebayes => {
                format!(
                    " {vcf}",
                    vcf=vcf.display()
                )
            }
        };

        run_command(&cmd)?;

        Ok(())
    }
    // Can ad the '-i -u -X' flags that may influence results from Bayesian variant caller 
    pub fn call_variants(&self, bam: &Path, vcf: &Path) -> Result<(), VircovError> {


        let cmd = match self.config.variant_caller {
            HaplotypeVariantCaller::Freebayes => {
                format!(
                    "freebayes -f {fasta} -F {min_var_freq} --pooled-continuous {bam} > {vcf}",
                    fasta=self.config.fasta.display(),
                    min_var_freq=self.config.min_var_freq,
                    bam=bam.display(),
                    vcf=vcf.display()
                )
            }
        };

        run_command(&cmd)?;

        Ok(())
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