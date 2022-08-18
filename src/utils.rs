use std::fs::File;
use anyhow::{Result, Error};
use std::collections::BTreeMap;
use crate::align::CoverageFields;
use crate::align::ReadAlignmentError;
use std::collections::btree_map::Entry;
use ordered_float::OrderedFloat;

pub fn get_sanitized_fasta_writer(name: &str, path: &std::path::PathBuf) -> Result<noodles::fasta::Writer<File>, Error> {

    let sanitized_name = name.replace(" ", "_");
    let sanitized_name = sanitized_name.trim_matches(';'); // Virosaurus specific: sanitize remaining header separator on sequence identifier (weird format)
    let file_path = path.join(&sanitized_name).with_extension("fasta");
    let file_handle = std::fs::File::create(file_path.as_path())
        .expect(&format!("Could not create fasta file"));

    Ok(noodles::fasta::Writer::new(file_handle))

}

/*
========================
ReadAlignment utilities
========================
*/ 


/// Get a dictionary / map of segments grouped by their unique fields:
/// 
/// If `segment_field` = "segment=" the dictionary groupings will be for 
/// example: "segment=L": [CoverageFields, ...]
pub fn get_grouped_segments(tags: Vec<CoverageFields>, segment_field: Option<String>, group_sep: String) -> Result<BTreeMap<String, Vec<CoverageFields>>, ReadAlignmentError>{

    match segment_field {
        Some(seg_field) => {
            let mut grouped_segments: BTreeMap<String, Vec<CoverageFields>> = BTreeMap::new();
            for tag_cov_field in tags {
                let header_fields = tag_cov_field.description.split(&group_sep);

                for field in header_fields {
                    if field.contains(&seg_field) {
                        match grouped_segments.entry(field.to_string()) {
                            Entry::Occupied(mut entry) => {
                                entry.get_mut().push(tag_cov_field.clone());
                            }
                            Entry::Vacant(entry) => {
                                entry.insert(vec![tag_cov_field.clone()]);
                            }
                        }
                    }
                }
            };
            Ok(grouped_segments)
        }, 
        None => return Err(ReadAlignmentError::GroupSelectSplitError())
    }

}


// Select the best segments from grouped references of the same segment identifier by reads or coverage
pub fn get_segment_selections(grouped_segments: BTreeMap<String, Vec<CoverageFields>>, group_select_by: Option<String>) -> Result<BTreeMap<String, CoverageFields>, ReadAlignmentError> {

    let mut selected_segments: BTreeMap<String, CoverageFields> = BTreeMap::new();

    for (segment, cov_fields) in grouped_segments {
        
        let selected = match group_select_by.clone() {
            Some(value) => match value.as_str() {
                "reads" => cov_fields.iter().max_by_key(|x| x.reads),
                "coverage" => cov_fields.iter().max_by_key(|x| OrderedFloat(x.coverage)),
                _ => return Err(ReadAlignmentError::GroupSelectByError())
            }, 
            None => return Err(ReadAlignmentError::GroupSelectByError())
        };

        match selected {
            Some(selected_cov_field) => {
                let sanitized_name = selected_cov_field.name.replace(" ", "_");
                let sanitized_name = sanitized_name.trim_matches(';');
                let segment_sanitized = segment.trim_matches(' ').to_string();  // trim whitespaces fro msegment field
                let seg_name = format!("{}_{}", sanitized_name, segment_sanitized);
                selected_segments.entry(seg_name).or_insert(selected_cov_field.clone());
            },
            None => return Err(ReadAlignmentError::GroupSelectReference())
        }

    };

    Ok(selected_segments)

}