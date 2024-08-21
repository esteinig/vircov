use std::{path::PathBuf, process::{Command, Output, Stdio}, str::from_utf8};

use needletail::{parse_fastx_file, parser::SequenceRecord, Sequence};
use serde::{Deserialize, Serialize};

use crate::{error::VircovError, vircov::ConsensusConfig};

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, clap::ValueEnum)]
pub enum ConsensusAssembler {
    Ivar
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct ConsensusRecord {
    pub id: String,
    pub length: u64,
    pub missing: u64,
    pub completeness: f64
}

impl std::fmt::Display for ConsensusAssembler {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ConsensusAssembler::Ivar => write!(f, "ivar"),
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct VircovConsensus {
    config: ConsensusConfig
}
impl VircovConsensus {
    pub fn new(config: ConsensusConfig) -> Result<Self, VircovError> {
        Self::check_consensus_dependency(&config.assembler)?;

        Ok(Self { config })
    }
    pub fn assemble(&self) -> Result<Vec<ConsensusRecord>, VircovError> {
        let records = match self.config.assembler {
            ConsensusAssembler::Ivar => self.run_ivar()?,
        };
        Ok(records)
    }
    pub fn check_consensus_dependency(assembler: &ConsensusAssembler) -> Result<(), VircovError> {
        let command = match assembler {
            ConsensusAssembler::Ivar => "ivar --version",
        };
        Self::run_version_command(command).map_err(|_| VircovError::ConsensusDependencyMissing(assembler.clone()))?;
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
    fn run_ivar(&self) -> Result<Vec<ConsensusRecord>, VircovError> {
        let assembler_args = self.config.args.as_deref().unwrap_or("");
        let mpileup_args = self.config.mpileup.as_deref().unwrap_or("");


        let samtools_ref = match &self.config.reference {
            Some(reference) => format!("samtools index '{}' && samtools view -h '{}' '{reference}' |", self.config.alignment.display(), self.config.alignment.display()),
            None => format!("")
        };

        let ivar_header = match &self.config.header {
            Some(header) => format!("-i '{header}'"),
            None => "".to_string(),
        };

        let cmd = format!(
            "{} samtools mpileup -a -A -d 0 -Q 0 {} - | ivar consensus -p '{}' -q {} -t {} -m {} -n {} {} {}", // Do not use -aa as segmented genomes have multiple ref-segment headers (@SQ)
            samtools_ref,
            mpileup_args,
            self.config.output.display(),
            self.config.min_quality,
            self.config.min_frequency,
            self.config.min_depth,
            self.config.missing,
            ivar_header,
            assembler_args,
        );

        self.run_command(&cmd)?;

        Ok(self.parse_consensus_sequence(&self.config.output)?)
    }
    pub fn parse_consensus_sequence(&self, fasta: &PathBuf) -> Result<Vec<ConsensusRecord>, VircovError>{
        
        let mut reader = match parse_fastx_file(fasta) {
            Ok(reader) => reader, Err(err) => {
                log::error!("{}", err.to_string());
                log::warn!("Possibly no consensus sequence generated in: {}", fasta.display());
                return Ok(Vec::new())
            }
        };

        let mut consensus_records = Vec::new();

        while let Some(record) = reader.next() {
            
            let rec = record?;
            let count = rec.sequence().iter().filter(|x| *x == &b'N').count() as u64;
            let rec_id = get_seq_record_identifier(&rec)?;
            let rec_len = rec.num_bases() as u64;
            let completeness  = 100.0 - (count as f64 / rec_len as f64)*100.0;

            consensus_records.push(ConsensusRecord {
                id: rec_id,
                length: rec_len,
                missing: count,
                completeness
            })
        }
        Ok(consensus_records)
    }
    fn run_command(&self, cmd: &str) -> Result<(), VircovError> {
        log::debug!("Running command: {}", cmd);

        let status = Command::new("sh")
            .arg("-c")
            .arg(cmd)
            .stderr(Stdio::null())
            .stdout(Stdio::null())
            .status()
            .map_err(|e| VircovError::CommandExecutionFailed(cmd.to_string(), e.to_string()))?;

        if !status.success() {
            return Err(VircovError::CommandFailed(cmd.to_string(), status.code().unwrap_or(-1)));
        }

        Ok(())
    }
}



pub fn get_seq_record_identifier(rec: &SequenceRecord) -> Result<String, VircovError> {
    match from_utf8(rec.id())?.split(' ').next() {
        Some(rec_id) => Ok(rec_id.to_string()),
        None => return Err(VircovError::NeedltailRecordIdentifierNotParsed)
    }
}