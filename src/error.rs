use thiserror::Error;

#[derive(Error, Debug)]
pub enum VircovError {
    #[error("Subtype database error")]
    SubtypeDatabase(#[from] crate::subtype::SubtypeDatabaseError),
    #[error("Read alignment error")]
    ReadAlignment(#[from] crate::align::ReadAlignmentError),
    #[error(transparent)]
    IOError(#[from] std::io::Error),
}