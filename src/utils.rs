use crate::error::VircovError;
use crate::vircov::Vircov;
use anyhow::Result;
use csv::{Reader, ReaderBuilder, Writer, WriterBuilder};
use env_logger::{fmt::Color, Builder};
use log::{Level, LevelFilter};
use needletail::{parse_fastx_file, FastxReader};
use niffler::{get_reader, get_writer};
use serde::{Deserialize, Serialize};
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};



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
            Some(Some("lzma") | Some("xz")) => Self::Lzma,
            _ => Self::No,
        }
    }
}

/// Enum to specify the type of file component to retrieve
pub enum FileComponent {
    /// The full file name including the extension
    FileName,
    /// The file name without the extension
    FileStem,
}

/// Extracts the specified file component from a `PathBuf` and returns it as a `String`.
///
/// # Arguments
///
/// * `path` - A `PathBuf` representing the file path.
/// * `component` - A `FileComponent` specifying whether to get the file name or the file stem.
///
/// # Returns
///
/// * `Result<String, DatabaseError>` - A `Result` containing the specified file component as a `String`
///   if successful, or a `DatabaseError` if an error occurs.
///
/// # Examples
///
/// ```
/// use std::path::PathBuf;
/// use cipher::{get_file_component, FileComponent};
///
/// let path = PathBuf::from("/some/path/to/file.txt");
/// match get_file_component(path, FileComponent::FileName) {
///     Ok(file_name) => println!("File name: {}", file_name),
///     Err(e) => eprintln!("Error: {}", e),
/// }
/// ```
///
/// ```
/// use std::path::PathBuf;
/// use cipher::{get_file_component, FileComponent};
///
/// let path = PathBuf::from("/some/path/to/file.txt");
/// match get_file_component(path, FileComponent::FileStem) {
///     Ok(file_stem) => println!("File stem: {}", file_stem),
///     Err(e) => eprintln!("Error: {}", e),
/// }
/// ```
pub fn get_file_component(path: &PathBuf, component: FileComponent) -> Result<String, VircovError> {
    match component {
        FileComponent::FileName => {
            path.file_name()
                .ok_or(VircovError::FileNameConversionError)
                .and_then(|os_str| os_str.to_str().map(String::from).ok_or(VircovError::FileNameConversionError))
        }
        FileComponent::FileStem => {
            path.file_stem()
                .ok_or(VircovError::FileNameConversionError)
                .and_then(|os_str| os_str.to_str().map(String::from).ok_or(VircovError::FileNameConversionError))
        }
    }
}

pub fn get_niffler_fastx_writer(output: &PathBuf) -> Result<Box<dyn std::io::Write>, VircovError> {

    let file: File = File::create(output)?;
    let file_handle = BufWriter::new(file);
    let fmt = niffler::Format::from_path(&output);

    let writer = niffler::get_writer(
        Box::new(file_handle), 
        fmt,
        niffler::compression::Level::Six
    )?;

    Ok(writer)
}

// Utility function to get a Needletail reader and Niffler compressed/uncompressed writer
pub fn get_niffler_fastx_reader_writer(
    input: &PathBuf,
    output: &PathBuf,
    compression_level: niffler::compression::Level,
    output_format: Option<niffler::compression::Format>,
) -> Result<(Box<dyn FastxReader>, Box<dyn std::io::Write>), VircovError> {

    // Input output of read files includes compression detection
    let reader = parse_fastx_file(input)?;

    let file: File = File::create(output)?;
    let file_handle = BufWriter::new(file);
    let fmt = match output_format {
        None => niffler::Format::from_path(output),
        Some(format) => format,
    };

    let writer = niffler::get_writer(Box::new(file_handle), fmt, compression_level)?;

    Ok((reader, writer))
        
}

pub fn is_file_empty<P: AsRef<Path>>(path: P) -> Result<bool, VircovError> {
    let file = File::open(&path)?;
    
    // Use niffler to get a reader for the (possibly compressed) file
    let (mut reader, _format) = match get_reader(Box::new(file)) {
        Ok(reader_format) => reader_format,
        Err(niffler::Error::FileTooShort) => return Ok(true),
        Err(e) => return Err(VircovError::NifflerError(e)),
    };
    // Try to read the first byte
    let mut buffer = [0; 1];
    match reader.read(&mut buffer) {
        Ok(0) => Ok(true),
        Ok(_) => Ok(false), // Successfully read a byte, file is not empty
        Err(e) => Err(VircovError::IOError(e))
    }
}

pub fn parse_fastx_file_with_check<P: AsRef<Path>>(path: P) -> Result<Option<Box<dyn FastxReader>>, VircovError> {
    if is_file_empty(&path)? {
        Ok(None)
    } else {
        Ok(Some(parse_fastx_file(&path)?))
    }
}

pub fn get_tsv_reader(file: &Path, flexible: bool, header: bool) -> Result<Reader<Box<dyn Read>>, VircovError> {

    let buf_reader = BufReader::new(File::open(&file)?);
    let (reader, _format) = get_reader(Box::new(buf_reader))?;

    let csv_reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(header)
        .flexible(flexible) // Allows records with a different number of fields
        .from_reader(reader);

    Ok(csv_reader)
}

pub fn get_tsv_writer(file: &Path,) -> Result<Writer<Box<dyn Write>>, VircovError> {
    
    let buf_writer = BufWriter::new(File::create(&file)?);
    let writer = get_writer(Box::new(buf_writer), niffler::Format::from_path(file), niffler::compression::Level::Six)?;

    let csv_writer = WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(writer);

    Ok(csv_writer)
}

pub fn write_tsv<T: Serialize>(data: &Vec<T>, file: &Path) -> Result<(), VircovError> {

    let mut writer = get_tsv_writer(file)?;

    for value in data {
        // Serialize each value in the vector into the writer
        writer.serialize(&value)?;
    }

    // Flush and complete writing
    writer.flush()?;
    Ok(())
}

pub fn read_tsv<T: for<'de>Deserialize<'de>>(file: &Path, flexible: bool, header: bool) -> Result<Vec<T>, VircovError> {

    let mut reader = get_tsv_reader(file, flexible, header)?;

    let mut records = Vec::new();
    for record in reader.deserialize() {
        records.push(record?)
    }

    Ok(records)
}


pub fn get_seq_id(id: &[u8]) -> Result<String, VircovError> {
    let header = std::str::from_utf8(id)?;
    let header_components = header.split_whitespace().collect::<Vec<&str>>();
    let id = header_components[0].to_string();

    Ok(id)
}
