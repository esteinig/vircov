use thiserror::Error;

use crate::alignment::Aligner;

#[derive(Error, Debug)]
pub enum VircovError {
    #[error("Subtype database error")]
    SubtypeDatabase(#[from] crate::subtype::SubtypeDatabaseError),
    #[error(transparent)]
    IOError(#[from] std::io::Error),
    /// Represents all other cases of `niffler::Error`.
    #[error(transparent)]
    NifflerError(#[from] niffler::Error),
    
    /// Represents an error when a command execution fails.
    #[error("Failed to execute command '{0}': {1}")]
    CommandExecutionFailed(String, String),
    /// Represents an error when a command exits with a non-zero status code.
    #[error("Command '{0}' exited with status code: {1}")]
    CommandFailed(String, i32),

    /// Represents an error when no aligner is configured.
    #[error("No aligner configured.")]
    MissingAligner,
    /// Represents an error when the specified aligner cannot be executed, possibly due to it not being installed.
    #[error("Aligner `{0}` cannot be executed - is it installed?")]
    AlignerDependencyMissing(Aligner),
    /// Represents an error when alignment index is not set while aligner is configured.
    #[error("Alignment index must be set when aligner is configured.")]
    MissingAlignmentIndex,
    /// Represents an error when alignment output is not set while alignment is configured.
    #[error("Alignment output must be set when alignment is configured.")]
    MissingAlignment,
    /// Represents an error when the alignment format is not explicitly set and not recognized from extension
    #[error("Unable to recognize alignment input format from extension.")]
    AlignmentInputFormatNotRecognized,
    /// Represents an error when the alignment format is explicitly set and not recognized
    #[error("Unable to recognize alignment input format - is this version compiled with 'htslib'?")]
    AlignmentInputFormatInvalid,
    /// Represents an error when no preset is configured.
    #[error("Minimap2 was set as aligner but no preset was configured.")]
    MissingMinimap2Preset,
    /// Represents an error when no preset is configured.
    #[error("Minigraph was set as aligner but no preset was configured.")]
    MissingMinigraphPreset,

    /// Indicates failure to read file from command line option
    #[error("failed to read file from option")]
    FileInputError,
    /// Indicates a failure to get sequence length require for coverage plots
    #[error("failed to get sequence length for coverage plot")]
    CovPlotSeqLengthError(),
    /// Indicates failure with the coverage plot module
    #[error("failed to generate data for coverage plot")]
    CovPlot(#[from] crate::covplot::CovPlotError),
    /// Indicates failure to parse a BAM file
    #[error("failed to parse records from BAM")]
    HTSLIBError(#[from] rust_htslib::errors::Error),
    /// Indicates failure to parse a record name from BAM file
    #[error("failed to parse record name from BAM")]
    UTF8Error(#[from] std::str::Utf8Error),
    /// Indicates failure to parse a target name from BAM file
    #[error("failed to parse a valid record target name from BAM")]
    TIDError(#[from] std::num::TryFromIntError),
    /// Indicates failure to parse an u64 from PAF
    #[error("failed to parse a valid integer from PAF")]
    PafRecordIntError(#[from] std::num::ParseIntError),
    #[error("failed to parse a valid integer from record")]
    FloatError(#[from] std::num::ParseFloatError),
    /// Indicates failure to conduct grouping because no reference sequences were parsed
    #[error("failed to group outputs due to missing reference sequences")]
    GroupSequenceError,
    /// Indicates failure to plot coverage when data is grouped
    #[error("coverage plots are not enabled when grouping output")]
    GroupCovPlotError,
    /// Indicates failure to infer or identify file format from explicit option
    #[error("failed to parse a valid input format")]
    InputFormatError,
    /// Indicates failure when no grouping options are provided and selection is specified
    #[error("failed to use the group selection options becuase no grouping is specified")]
    GroupSelectSplitError,
    /// Indicates failure when no group select by is provided (should not occurr)
    #[error("failed to group select by reads or coverage")]
    GroupSelectByError,
    /// Indicates failure when no group select by is provided (should not occurr)
    #[error("failed to provide a negative segment field")]
    SegmentFieldNaNError,
    /// Indicates failure when no best value from the grouped coverage fields could be selected
    #[error("failed to select a best reference sequence")]
    GroupSelectReference,
    /// Indicates failure when no reference name could be selected from the grouped identifier
    #[error("failed to extract a reference name from the grouped identifier")]
    GroupSelectReferenceName,
    /// Indicates failure when no coverage value could be extracted from a group of coverage fields
    #[error("failed to extract the highest coverage value from the grouped fields")]
    GroupSelectCoverage,
    /// Indicates failure when zero-value reference sequences should be included, but no reference sequences were provided
    #[error("no reference sequences found, zero-value records cannot be included")]
    ZeroReferenceSequences,
    
}
