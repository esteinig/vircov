use std::{collections::HashMap, path::{Path, PathBuf}};

use needletail::parser::write_fasta;
use serde::{Deserialize, Serialize};

use crate::{error::VircovError, utils::{get_niffler_fastx_writer, get_seq_id, parse_fastx_file_with_check, read_tsv, write_tsv}};

#[derive(Debug, Clone, Serialize, Deserialize, clap::ValueEnum)]
pub enum AnnotationPreset {
    Virosaurus,
    Ictv,
    Default
}

/// Database annotation for Vircov

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AnnotationRecord {
    id: String,
    bin: String,
    segment: Option<String>,
    description: Option<String>
}

pub struct AnnotationConfig {
    pub field_delimiter: String, 
    pub value_delimiter: String, 
    pub bin: String, 
    pub segment: String, 
    pub segment_nan: String, 
    pub description: String
}
impl AnnotationConfig {
    pub fn from_preset(preset: AnnotationPreset) -> Self {
        match preset {
            AnnotationPreset::Default => Self::default(),
            AnnotationPreset::Ictv => Self::ictv(),
            AnnotationPreset::Virosaurus => Self::virosaurus()
        }
    }
    pub fn virosaurus() -> Self {
        Self {
            field_delimiter: String::from(";"),
            value_delimiter: String::from("="),
            bin: String::from("taxid"),
            segment: String::from("segment"),
            segment_nan: String::from("N/A"),
            description: String::from("description")
        }
    }
    pub fn ictv() -> Self {
        Self {
            field_delimiter: String::from("|"),
            value_delimiter: String::from("="),
            bin: String::from("species"),
            segment: String::from("segment"),
            segment_nan: String::from("NAN"),
            description: String::from("description")
        }
    }
    pub fn str_from_record(&self, record: &AnnotationRecord) -> String {

        format!(
            "{bin}{vsep}{vbin}{fsep}{segment}{vsep}{vsegment}{fsep}{description}{vsep}{vdescription}",
            bin=self.bin,
            vsep=self.value_delimiter,
            vbin=record.bin,
            fsep=self.field_delimiter,
            segment=self.segment,
            vsegment=record.segment.clone().unwrap_or(self.segment_nan.clone()),
            description=self.description,
            vdescription=record.description.clone().unwrap_or(String::from(""))
        )

    }
}
impl Default for AnnotationConfig {
    fn default() -> Self {
        Self {
            field_delimiter: String::from("|"),
            value_delimiter: String::from("="),
            bin: String::from("bin"),
            segment: String::from("segment"),
            segment_nan: String::from("NAN"),
            description: String::from("description")
        }
    }
}

pub struct DatabaseAnnotation {
    config: AnnotationConfig,
    annotations: Vec<AnnotationRecord>
}
impl DatabaseAnnotation {
    pub fn new(tsv: &Path, config: AnnotationConfig) -> Result<Self, VircovError> {

        Ok(Self {
            config,
            annotations: read_tsv(&tsv, false)?,
        })
    }
    pub fn annotate(&self, fasta_in: &PathBuf, fasta_out: &PathBuf, skipped_out: Option<PathBuf>) -> Result<(), VircovError> {
 
        let annotation_map: HashMap<String, AnnotationRecord> = HashMap::from_iter(
            self.annotations.iter().map(|a| (a.id.clone(), a.clone())).collect::<Vec<_>>()
        );

        let mut reader = match parse_fastx_file_with_check(&fasta_in)? {
            Some(reader) => reader, 
            None => {
                log::error!("Failed to read database sequence file - is it empty?");
                return Err(VircovError::FastaFileIsEmpty(fasta_in.to_path_buf()))
            }
        };

        let mut writer = get_niffler_fastx_writer(&fasta_out)?;
        
        let mut skipped = Vec::new();
        while let Some(rec) = reader.next() {
            let record = rec?;
            let id = get_seq_id(record.id())?;

            match annotation_map.get(&id) {
                Some(annotation) => {
                    let new_id = format!("{id} {}", self.config.str_from_record(annotation));
                    write_fasta(new_id.as_bytes(), record.raw_seq(), &mut writer, record.line_ending())?;
                },
                None => {
                    log::warn!("Failed to find annotation for sequence {id} (skipped)");
                    skipped.push(id)
                }
            }
        }
        if let Some(path) = skipped_out {

            log::info!("Skipped {} sequences for which annotations were not provided", skipped.len());
            write_tsv(&skipped, &path)?;
            log::info!("Skipped sequence identifiers written to: {}", path.display());
        }

        Ok(())

    }
}
