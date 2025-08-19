use std::{collections::HashMap, path::{Path, PathBuf}};

use needletail::parser::write_fasta;
use serde::{Deserialize, Serialize};

use crate::{error::VircovError, utils::{get_niffler_fastx_writer, get_seq_id, parse_fastx_file_with_check, read_tsv, write_tsv}};

#[derive(Debug, Clone, Serialize, Deserialize, clap::ValueEnum)]
pub enum AnnotationPreset {
    Virosaurus,
    Ictv,
    Default,
    Cipher
}

/// Database annotation for Vircov

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AnnotationRecord {
    pub id: String,
    pub bin: String,
    pub segment: Option<String>,
    pub name: Option<String>,
    pub description: Option<String>,
    pub taxid: Option<String>,
    pub genome: Option<String>
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AnnotationConfig {
    pub field_delimiter: String, 
    pub value_delimiter: String, 
    pub bin: String, 
    pub segment: String, 
    pub segment_nan: String, 
    pub genome: String,
    pub genome_nan: String,
    pub name: String,
    pub description: String,
    pub taxid: String
}
impl AnnotationConfig {
    pub fn from_preset(preset: AnnotationPreset) -> Self {
        match preset {
            AnnotationPreset::Default => Self::default(),
            AnnotationPreset::Ictv => Self::ictv(),
            AnnotationPreset::Virosaurus => Self::virosaurus(),
            AnnotationPreset::Cipher => Self::cipher()
        }
    }
    pub fn virosaurus() -> Self {
        Self {
            field_delimiter: String::from(";"),
            value_delimiter: String::from("="),
            bin: String::from("taxid"),
            segment: String::from("segment"),
            name: String::from("usual_name"),
            segment_nan: String::from("N/A"),
            description: String::from("description"),
            taxid: String::from("taxid"),
            genome: String::from("genome"),
            genome_nan: String::from("N/A")
        }
    }
    pub fn ictv() -> Self {
        Self {
            field_delimiter: String::from("|"),
            value_delimiter: String::from("="),
            bin: String::from("species"),
            segment: String::from("segment"),
            segment_nan: String::from("NAN"),
            name: String::from("name"),
            description: String::from("description"),
            taxid: String::from("taxid"),
            genome: String::from("genome"),
            genome_nan: String::from("NAN")
        }
    }
    pub fn cipher() -> Self {
        Self {
            field_delimiter: String::from("|"),
            value_delimiter: String::from("="),
            bin: String::from("species"),
            segment: String::from("segment"),
            segment_nan: String::from("NAN"),
            name: String::from("name"),
            description: String::from("description"),
            taxid: String::from("taxid"),
            genome: String::from("genome"),
            genome_nan: String::from("NAN")
        }
    }
    pub fn str_from_record(&self, record: &AnnotationRecord) -> String {
        format!(
            "{bin}{vsep}{vbin}{fsep}{segment}{vsep}{vsegment}{fsep}{name}{vsep}{vname}{fsep}{description}{vsep}{vdescription}{fsep}{taxid}{vsep}{vtaxid}{fsep}{genome}{vsep}{vgenome}",
            bin=self.bin,
            vsep=self.value_delimiter,
            vbin=record.bin,
            fsep=self.field_delimiter,
            segment=self.segment,
            vsegment=record.segment.clone().unwrap_or(self.segment_nan.clone()),
            name=self.name,
            vname=record.name.clone().unwrap_or(String::from("")),
            description=self.description,
            vdescription=record.description.clone().unwrap_or(String::from("")),
            taxid=self.taxid,
            vtaxid=record.taxid.clone().unwrap_or(String::from("")),
            genome=self.genome,
            vgenome=record.genome.clone().unwrap_or(self.genome_nan.clone())
        )
    }
    pub fn sanitize(&self, s: &str) -> String {
        s.replace(" ", "_").replace("|", "_").replace(";", "_").replace("'", "_")
    }
    pub fn segment_is_nan(&self, segment: &str) -> bool {
        segment == self.segment_nan
    }
    pub fn genome_is_nan(&self, genome: &str) -> bool {
        genome == self.genome_nan
    }
    pub fn bin_field(&self) -> String {
        format!("{}{}", self.bin, self.value_delimiter)
    }
    pub fn segment_field(&self) -> String {
        format!("{}{}", self.segment, self.value_delimiter)
    }
    pub fn genome_field(&self) -> String {
        format!("{}{}", self.genome, self.value_delimiter)
    }
    pub fn name_field(&self) -> String {
        format!("{}{}", self.name, self.value_delimiter)
    }
    pub fn description_field(&self) -> String {
        format!("{}{}", self.description, self.value_delimiter)
    }
    pub fn taxid_field(&self) -> String {
        format!("{}{}", self.taxid, self.value_delimiter)
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
            name: String::from("name"),
            description: String::from("description"),
            taxid: String::from("taxid"),
            genome: String::from("genome"),
            genome_nan: String::from("NAN"),
        }
    }
}

// Used in processing
pub struct Annotation {
    pub bin: Option<String>,
    pub name: Option<String>,
    pub segment: Option<String>,
    pub genome: Option<String>,
    pub taxid: Option<String>
}
impl Annotation {
    pub fn from(description: &str, options: &AnnotationConfig) -> Self {

        let fields: Vec<&str> = description.split(&options.field_delimiter).map(|field| field.trim()).collect();

        let bin_fields: Vec<&str> = fields.iter().filter(|field| {
            field.starts_with(&options.bin_field())
        }).map(|x| *x).collect();
        
        let bin = match bin_fields.first() {
            Some(f) => Some(f.trim().trim_start_matches(&options.bin_field()).to_string()),
            _ => None
        };

        let name_fields: Vec<&str> = fields.iter().filter(|field| {
            field.starts_with(&options.name_field())
        }).map(|x| *x).collect();
        
        let name = match name_fields.first() {
            Some(f) => Some(f.trim().trim_start_matches(&options.name_field()).to_string()),
            _ => None
        };

        let segment_fields: Vec<&str> = fields.iter().filter(|field| {
            field.starts_with(&options.segment_field())
        }).map(|x| *x).collect();
        
        let segment = match segment_fields.first() {
            Some(f) => Some(f.trim().trim_start_matches(&options.segment_field()).to_string()),
            _ => None
        };

        let genome_fields: Vec<&str> = fields.iter().filter(|field| {
            field.starts_with(&options.genome_field())
        }).map(|x| *x).collect();
        
        let genome = match genome_fields.first() {
            Some(f) => Some(f.trim().trim_start_matches(&options.genome_field()).to_string()),
            _ => None
        };

        let taxid_fields: Vec<&str> = fields.iter().filter(|field| {
            field.starts_with(&options.taxid_field())
        }).map(|x| *x).collect();
        
        let taxid = match taxid_fields.first() {
            Some(f) => Some(f.trim().trim_start_matches(&options.taxid_field()).to_string()),
            _ => None
        };

        Self {
            bin, name, segment, genome, taxid
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
            annotations: read_tsv(&tsv, false, true)?,
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
