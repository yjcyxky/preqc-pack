//! `PreQC-Pack` is a suite of programs for pre-quality control.
#[macro_use]
extern crate lazy_static;

pub mod fastqc;
pub mod hasher;
pub mod util;
pub mod qc;

use regex::Regex;
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
pub struct QCPack {
  fastqc: fastqc::FastQC,
  filemeta: hasher::Meta,
}

impl QCPack {
  pub fn new() -> QCPack {
    return QCPack {
      fastqc: fastqc::init_fastqc(0),
      filemeta: hasher::init_meta(),
    };
  }

  pub fn set_fastqc(&mut self, fastqc: fastqc::FastQC) {
    self.fastqc = fastqc;
  }

  pub fn set_filemeta(&mut self, filemeta: hasher::Meta) {
    self.filemeta = filemeta;
  }
}

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