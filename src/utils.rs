use crate::alignment::Coverage;
use crate::error::VircovError;
use anyhow::Result;
use env_logger::{fmt::Color, Builder};
use log::{Level, LevelFilter};
use ordered_float::OrderedFloat;
use std::collections::btree_map::Entry;
use std::collections::BTreeMap;
use std::ffi::OsStr;
use std::fs::File;
use std::io::Write;
use std::path::Path;

pub fn get_sanitized_fasta_writer(
    name: &str,
    path: &std::path::Path,
) -> Result<noodles::fasta::Writer<File>, VircovError> {
    let sanitized_name = name.replace(' ', "_");
    let file_path = path.join(sanitized_name).with_extension("fasta");
    let file_handle = std::fs::File::create(file_path.as_path())?;

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
pub fn get_grouped_segments(
    tags: Vec<Coverage>,
    segment_field: Option<String>,
    group_sep: String,
) -> Result<BTreeMap<String, Vec<Coverage>>, VircovError> {
    match segment_field {
        Some(seg_field) => {
            let mut grouped_segments: BTreeMap<String, Vec<Coverage>> = BTreeMap::new();
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
            }
            Ok(grouped_segments)
        }
        None => Err(VircovError::GroupSelectSplitError),
    }
}

// Select the best segments from grouped references of the same segment identifier by reads or coverage
pub fn get_segment_selections(
    grouped_segments: BTreeMap<String, Vec<Coverage>>,
    select_by: String,
) -> Result<BTreeMap<String, Coverage>, VircovError> {
    let mut selected_segments: BTreeMap<String, Coverage> = BTreeMap::new();

    for (segment, cov_fields) in grouped_segments {
        let selected = match select_by.as_str() {
            "reads" => cov_fields.iter().max_by_key(|x| x.reads),
            "coverage" => cov_fields.iter().max_by_key(|x| OrderedFloat(x.coverage)),
            _ => return Err(VircovError::GroupSelectByError),
        };

        match selected {
            Some(selected_cov_field) => {
                let sanitized_name = selected_cov_field.name.replace(' ', "_");
                let sanitized_name = sanitized_name.trim_matches(';');
                let segment_sanitized = segment.trim_matches(' ').to_string(); // trim whitespaces from segment field
                let seg_name = format!("{}_{}", sanitized_name, segment_sanitized);
                selected_segments
                    .entry(seg_name)
                    .or_insert_with(|| selected_cov_field.clone());
            }
            None => return Err(VircovError::GroupSelectReference),
        }
    }

    Ok(selected_segments)
}

pub fn init_logger() {
    Builder::new()
        .format(|buf, record| {
            let timestamp = buf.timestamp();

            let mut red_style = buf.style();
            red_style.set_color(Color::Red).set_bold(true);
            let mut green_style = buf.style();
            green_style.set_color(Color::Green).set_bold(true);
            let mut white_style = buf.style();
            white_style.set_color(Color::White).set_bold(false);
            let mut orange_style = buf.style();
            orange_style
                .set_color(Color::Rgb(255, 102, 0))
                .set_bold(true);
            let mut apricot_style = buf.style();
            apricot_style
                .set_color(Color::Rgb(255, 195, 0))
                .set_bold(true);

            let msg = match record.level() {
                Level::Warn => (
                    orange_style.value(record.level()),
                    orange_style.value(record.args()),
                ),
                Level::Info => (
                    green_style.value(record.level()),
                    white_style.value(record.args()),
                ),
                Level::Debug => (
                    apricot_style.value(record.level()),
                    apricot_style.value(record.args()),
                ),
                Level::Error => (
                    red_style.value(record.level()),
                    red_style.value(record.args()),
                ),
                _ => (
                    white_style.value(record.level()),
                    white_style.value(record.args()),
                ),
            };

            writeln!(
                buf,
                "{} [{}] - {}",
                white_style.value(timestamp),
                msg.0,
                msg.1
            )
        })
        .filter(None, LevelFilter::Info)
        .init();
}

pub trait CompressionExt {
    fn from_path<S: AsRef<OsStr> + ?Sized>(p: &S) -> Self;
}

/// Attempts to infer the compression type from the file extension.
/// If the extension is not known, then Uncompressed is returned.
impl CompressionExt for niffler::compression::Format {
    fn from_path<S: AsRef<OsStr> + ?Sized>(p: &S) -> Self {
        let path = Path::new(p);
        match path.extension().map(|s| s.to_str()) {
            Some(Some("gz")) => Self::Gzip,
            Some(Some("bz") | Some("bz2")) => Self::Bzip,
            Some(Some("lzma") | Some(".xz")) => Self::Lzma,
            _ => Self::No,
        }
    }
}
