//! `PreQC-Pack` is a suite of programs for pre-quality control.
#[macro_use]
extern crate lazy_static;

pub mod qc;
pub mod util;

use regex::Regex;

pub fn is_fastq_file(filepath: &str) -> bool {
  // Import at the crate root - preqc-pack.rs
  lazy_static! {
    static ref RE: Regex = Regex::new(".*(.fq|.fastq)$").unwrap();
  }

  RE.is_match(filepath)
}

pub fn is_fastq_gz_file(filepath: &str) -> bool {
  // Import at the crate root - preqc-pack.rs
  lazy_static! {
    static ref RE: Regex = Regex::new(".*(.fq.gz|.fastq.gz)$").unwrap();
  }

  RE.is_match(filepath)
}
