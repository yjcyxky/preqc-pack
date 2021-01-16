extern crate exitcode;
extern crate preqc_pack;
extern crate regex;
extern crate serde_json;

use log::*;
use preqc_pack::{fastqc, hasher};
use regex::Regex;
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::Path;
use structopt::StructOpt;

use blake2::Blake2b;
use md5::Md5;

/// A collection of metadata, such as file size, md5sum
#[derive(StructOpt, PartialEq, Debug)]
#[structopt(setting=structopt::clap::AppSettings::ColoredHelp, name="PreQC Tool Suite - Hasher", author="Jingcheng Yang <yjcyxky@163.com>")]
pub struct Arguments {
  /// Bam file to process
  #[structopt(name = "FILE")]
  input: String,

  /// A hash algorithms for output file.
  #[structopt(name="algorithm", short="m", long="algorithm", possible_values=&["md5sum", "blake2b"], default_value="md5sum")]
  algorithm: String,
}

#[derive(Serialize, Deserialize)]
struct QCPack {
  fastqc: fastqc::FastQC,
  filemeta: hasher::Meta,
}

fn is_fastq_file(filepath: &str) -> bool {
  // Import at the crate root - preqc-pack.rs
  lazy_static! {
    static ref RE: Regex = Regex::new(".*(.fq|.fastq)$").unwrap();
  }

  RE.is_match(filepath)
}

fn is_fastq_gz_file(filepath: &str) -> bool {
  // Import at the crate root - preqc-pack.rs
  lazy_static! {
    static ref RE: Regex = Regex::new(".*(.fq.gz|.fastq.gz)$").unwrap();
  }

  RE.is_match(filepath)
}

pub fn run(args: &Arguments) {
  if Path::new(&args.input).exists() {
    // TODO: Multi threads?
    // Get filemeta
    let mut file = fs::File::open(&args.input).unwrap();
    let meta = match &args.algorithm[..] {
      "blake2b" => hasher::process::<Blake2b, _>(&mut file),
      _ => hasher::process::<Md5, _>(&mut file),
    };

    let mut fastqc_metrics = fastqc::init_fastqc(0);
    if is_fastq_file(&args.input) {
      // Generate fastqc metrics
      fastqc_metrics = fastqc::compute_data_size(&args.input);
    } else if is_fastq_gz_file(&args.input) {
      // fastqc_metrics = fastqc::compute_gz_data_size(&args.input);
      fastqc_metrics = fastqc::compute_data_size_par(&args.input);
    } else {
      error!("Not a valid fastq/fastq.gz file")
    }

    let qc_pack = QCPack {
      fastqc: fastqc_metrics,
      filemeta: meta,
    };

    println!("{}", serde_json::to_string(&qc_pack).unwrap());
  } else {
    error!("{} - Not Found: {:?}", module_path!(), args.input);
  }
}
