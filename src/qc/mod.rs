pub mod fastqc;
pub mod util;

use fastq::parse_path;
use serde_json;
use std::io::Error;

pub fn init_fastqc() -> fastqc::FastQC {
  return fastqc::FastQC::new();
}

pub fn run_fastqc(fastq_path: &str) -> fastqc::FastQC {
  match parse_path(Some(fastq_path), |parser| {
    let result: Result<Vec<_>, Error> = parser.parallel_each(5, move |record_sets| {
      let mut qc = fastqc::FastQC::new();
      for record_set in record_sets {
        for record in record_set.iter() {
          qc.process_sequence(&record.to_owned_record());
        }
      }

      return qc;
    });

    return match result {
      Ok(o) => {
        if o.len() > 0 {
          let mut first = o[0].clone();
          first.merge(&o[1..]);
          first
          // println!("{}", serde_json::to_string_pretty(&first).unwrap());
        } else {
          init_fastqc()
        }
      }
      Err(msg) => {
        panic!("Parse fastq file error: {}", msg);
      }
    };
  }) {
    Err(msg) => {
      panic!("Cannot parse fastq file: {}", msg);
    }
    Ok(o) => o,
  }
}
