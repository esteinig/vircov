use anyhow::Result;
use rust_lapper::{Interval, Lapper};
use thiserror::Error;

use std::io::{stdout, Write};
use crossterm::{execute, Result as CrosstermResult};
use crossterm::style::{Print, SetForegroundColor, SetBackgroundColor, ResetColor, Color, Attribute};


/*
========================
Custom error definitions
========================
*/

#[derive(Error, Debug)]
pub enum CovPlotError {
}


/*
=================
Coverage plotter
=================
*/

#[derive(Debug, Clone)]
pub struct CovPlot {
    
}

impl CovPlot {
    /// Create a new coverage plot instance and
    /// computes the conversions from sequence 
    /// to output string
    pub fn new(seq_length: u64, max_width: u64) -> Result<()> {




        let _ = execute!(
            stdout(),
            Print("5'"),
            Print("------------------------------------------------"),
            // Blue foreground
            SetForegroundColor(Color::DarkRed),
            // Print text
            Print("--".to_string()),
            // Reset to default colors
            ResetColor,
            Print("------------------------------------------------"),
            Print("3'"),
        );
        println!();
        Ok(())
    }
    // Computes the segment length for terminal output
    fn get_segment_length(seq_length: u64, max_width: u64){
        
    }

}
