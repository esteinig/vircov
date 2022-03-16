use anyhow::Result;
use rust_lapper::{Interval, Lapper};
use thiserror::Error;

use std::io::stdout;
use crossterm::execute;
use crossterm::style::{Print, SetForegroundColor, SetAttribute, ResetColor, Color, Attribute};


/*
========================
Custom error definitions
========================
*/

#[derive(Error, Debug)]
pub enum CovPlotError {
    /// Indicates failure to specify maximum plot width
    #[error("sequence width cannot be zero")]
    WidthError(),
    /// Indicates failure to execute a crossterm macro
    #[error("failed to write a crossterm command")]
    CrosstermError(#[from] std::io::Error)
}


/*
=================
Coverage plotter
=================
*/

#[derive(Debug, Clone)]
pub struct CovPlot {
    pub segments: Vec<Interval<usize, u64>>
}

impl CovPlot {
    /// Create a new coverage plot instance and
    /// computes the conversions from sequence 
    /// to output string
    pub fn new(targets: Lapper<usize, String>, seq_length: u64, width: u64) -> Result<Self, CovPlotError> {

        let segment_length = match width {
            0 => return Err(CovPlotError::WidthError()),
            _ => (seq_length / width) as usize
        };

        let segment_intervals = (0..seq_length).step_by(segment_length).map(
            |x| Interval{start: x as usize, stop: x as usize + segment_length, val: 0}
        ).collect::<Vec<Interval<usize, u64>>>();

        // For each segment interval, find the merged overlap intervals from the coverage alignment
        // that overlap start .. stop of the segment
        let mut merged_targets = targets.clone();
        merged_targets.merge_overlaps();

        // For each segment interval, use the merged target coverage regions to determine whether they are
        // contained within the segment interval: 
        //      0 indicates it is not contained
        //      1 indicates it is contained
        //      2 indicates a coundary between two merged target regions
        let tagged_segments = segment_intervals.iter().map(|x| {
            let _count = merged_targets.find(x.start, x.stop).count();
            let mut _x = x.clone();
            if _count > 0 {
                _x.val = _count as u64;
            }
            _x
        }).collect::<Vec<Interval<usize, u64>>>();

        Ok(Self { segments: tagged_segments })
    }

    // Prints the five prime end to console in default style
    pub fn to_console(&self, seq_name: String, seq_length: u64, coverage_color: Color) -> Result<(), CovPlotError>{
        println!("{}  - {} bp", seq_name, seq_length);
        execute!(
            stdout(),
            Print("5'")
        )?;
        for segment in &self.segments {
            let (segment_color, segment_attribute) = match segment.val {
                0 => (Color::Reset, Attribute::Reset),
                1 => (coverage_color, Attribute::Reset),
                2 => (coverage_color, Attribute::Bold),
                _ => (Color::Yellow, Attribute::Bold)
            };
            execute!(
                stdout(),
                SetAttribute(segment_attribute),
                SetForegroundColor(segment_color),
                Print("-")
            )?;
        }
        execute!(
            stdout(),
            ResetColor,
            Print("3'")
        )?;

        println!();
        println!();
        Ok(())
    }
    
}