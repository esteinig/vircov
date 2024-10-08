use crate::vircov::HaplotypeConfig;
use crate::error::VircovError;

use serde::{Deserialize, Serialize};
use std::process::{Command, Output};

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, clap::ValueEnum)]
pub enum HaplotypeVariantCaller {
    Freebayes
}

pub struct HaplotypeRecord {
    id: String,
    snp: usize,

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
            HaplotypeVariantCaller::Freebayes => self.run_haplotyping()?,
        };
        Ok(records)
    }
    pub fn check_haplotype_dependencies() -> Result<(), VircovError> {
        Self::run_version_command("freebayes --version").map_err(|_| VircovError::HaplotypeDependencyMissing(String::from("freebayes")))?;
        Self::run_version_command("floria --version").map_err(|_| VircovError::HaplotypeDependencyMissing(String::from("floria")))?;
        Ok(())
    }
    fn run_version_command(command: &str) -> Result<Output, VircovError> {
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
    fn run_haplotyping(&self) -> Result<Vec<HaplotypeRecord>, VircovError> {

        

        Ok(Vec::new())
    }
}