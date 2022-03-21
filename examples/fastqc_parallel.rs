extern crate preqc_pack;
use fastq::parse_path;
use preqc_pack::qc::fastqc::FastQC;
use std::io::Error;
use std::sync::{Arc, Mutex};

fn test_process_sequence() {
  let fastq_path = "/Users/choppy/Downloads/R20059384_merge_R1.fastq.gz";
  parse_path(Some(fastq_path), |parser| {
    let result: Result<Vec<_>, Error> = parser.parallel_each(5, move |record_sets| {
      let mut qc = FastQC::new();
      for record_set in record_sets {
        for record in record_set.iter() {
          qc.process_sequence(&record.to_owned_record());
        }
      }

      return qc;
    });

    match result {
      Ok(o) => {
        println!("{:?}", o);
      }
      Err(msg) => {}
    };
    return true;
  })
  .unwrap();
}

fn main() {
  test_process_sequence();
}
