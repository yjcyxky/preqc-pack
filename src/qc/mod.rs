pub mod fastqc;
pub mod util;

use fastq::parse_path;
use serde_json;
use std::io::Error;

pub fn run_fastqc(fastq_path: &str) {
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

    match result {
      Ok(o) => {
        if o.len() > 0 {
          let mut first = o[0].clone();
          first.merge(&o[1..]);
          println!("{}", serde_json::to_string_pretty(&first).unwrap());
        } else {
          println!("No sequences")
        }
      }
      Err(msg) => {
        eprintln!("Parse fastq file error: {}", msg);
      }
    };
    return true;
  }) {
    Err(msg) => {
      eprintln!("Cannot parse fastq file: {}", msg);
    }
    _ => {}
  };
}
