use crate::error::VircovError;
use crate::vircov::Vircov;
use anyhow::Result;
use env_logger::{fmt::Color, Builder};
use log::{Level, LevelFilter};
use needletail::{parse_fastx_file, FastxReader};
use niffler::get_reader;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufWriter, Write};
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
            Some(Some("lzma") | Some(".xz")) => Self::Lzma,
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

