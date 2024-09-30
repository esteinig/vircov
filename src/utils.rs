use crate::error::VircovError;
use anyhow::Result;
use env_logger::{fmt::Color, Builder};
use log::{Level, LevelFilter};
use std::ffi::OsStr;
use std::io::Write;
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

