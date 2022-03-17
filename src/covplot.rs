use anyhow::Result;
use rust_lapper::{Interval, Lapper};
use thiserror::Error;

use crossterm::execute;
use crossterm::style::{Attribute, Color, Print, ResetColor, SetAttribute, SetForegroundColor};
use std::io::stdout;

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
    CrosstermError(#[from] std::io::Error),
}

/*
=================
Coverage plotter
=================
*/

#[derive(Debug, Clone)]
pub struct CovPlot {
    pub segments: Vec<Interval<usize, u64>>,
}

impl CovPlot {
    /// Create a new coverage plot instance and
    /// computes the conversions from sequence
    /// to output string
    pub fn new(
        targets: Lapper<usize, String>,
        seq_length: u64,
        width: u64,
    ) -> Result<Self, CovPlotError> {
        let segment_length = match width {
            0 => return Err(CovPlotError::WidthError()),
            _ => (seq_length / (width - 1)) as usize, // make the segments  slighlty larger to account for integer division
        };

        let segment_intervals = (0..seq_length)
            .step_by(segment_length)
            .map(|x| Interval {
                start: x as usize,
                stop: x as usize + segment_length,
                val: 0,
            })
            .collect::<Vec<Interval<usize, u64>>>();

        // For each segment interval, find the merged overlap intervals from the coverage alignment
        // that overlap start .. stop of the segment
        let mut merged_targets = targets;
        merged_targets.merge_overlaps();

        // For each segment interval, use the merged target coverage regions to determine whether they are
        // contained within the segment interval:
        //      0 indicates it is not contained
        //      1 indicates it is contained
        //      2 indicates a coundary between two merged target regions
        let tagged_segments = segment_intervals
            .iter()
            .map(|x| {
                let _count = merged_targets.find(x.start, x.stop).count();
                let mut _x = x.clone();
                if _count > 0 {
                    _x.val = _count as u64;
                }
                _x
            })
            .collect::<Vec<Interval<usize, u64>>>();

        Ok(Self {
            segments: tagged_segments,
        })
    }

    #[cfg(not(tarpaulin_include))]
    // Prints the segmented coverage plot
    pub fn to_console(
        &self,
        seq_name: String,
        seq_length: u64,
        coverage_color: Color,
    ) -> Result<(), CovPlotError> {
        println!("{}  - {} bp", seq_name, seq_length);
        execute!(stdout(), Print("5'"))?;
        for segment in &self.segments {
            let (segment_color, segment_attribute) = match segment.val {
                0 => (Color::Reset, Attribute::Reset),
                1 => (coverage_color, Attribute::Reset),
                2 => (coverage_color, Attribute::Bold),
                _ => (Color::Yellow, Attribute::Bold),
            };
            execute!(
                stdout(),
                SetAttribute(segment_attribute),
                SetForegroundColor(segment_color),
                Print("-")
            )?;
        }
        execute!(stdout(), ResetColor, Print("3'"))?;

        println!();
        println!();
        Ok(())
    }
}

#[cfg(test)]
#[cfg(not(tarpaulin_include))]
mod tests {

    use super::*;

    struct TestCases {
        // PAF test alignment intervals for L segment sequence
        paf_test_intervals_l_segment: Vec<Interval<usize, String>>,
    }

    impl TestCases {
        fn new() -> Self {
            Self {
                paf_test_intervals_l_segment: vec![
                    Interval {
                        start: 1786,
                        stop: 1834,
                        val: "FS10001392:17:BPL20314-1135:1:1113:8980:1660".to_string(),
                    },
                    Interval {
                        start: 4538,
                        stop: 4665,
                        val: "FS10001392:17:BPL20314-1135:1:1101:5600:2170".to_string(),
                    },
                    Interval {
                        start: 4758,
                        stop: 4904,
                        val: "FS10001392:17:BPL20314-1135:1:1101:5600:2170".to_string(),
                    },
                ],
            }
        }
    }

    #[test]
    fn covplot_create_new_ok() {
        let test_cases = TestCases::new();
        let lapper = Lapper::new(test_cases.paf_test_intervals_l_segment);
        let covplot = CovPlot::new(lapper, 7194, 100).unwrap();

        let _covered = covplot.segments.iter().filter(|x| x.val > 0).count();

        assert_eq!(covplot.segments.len(), 100);
        assert_eq!(
            covplot.segments[0],
            Interval {
                start: 0,
                stop: 72,
                val: 0
            }
        );
        assert_eq!(_covered, 7)
    }

    #[test]
    #[should_panic]
    fn covplot_create_new_width_fail() {
        let test_cases = TestCases::new();
        let lapper = Lapper::new(test_cases.paf_test_intervals_l_segment);
        let _ = CovPlot::new(lapper, 7194, 0).unwrap();
    }
}
