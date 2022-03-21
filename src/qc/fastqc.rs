use fastq::{OwnedRecord, Record};
use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize)]
pub struct FastQC {
  total_reads: usize,
  total_bases: usize,
  t_count: usize,
  c_count: usize,
  g_count: usize,
  a_count: usize,
  n_count: usize,
}

impl FastQC {
  pub fn new() -> FastQC {
    return FastQC {
      total_reads: 0,
      total_bases: 0,
      t_count: 0,
      c_count: 0,
      g_count: 0,
      a_count: 0,
      n_count: 0,
    };
  }

  pub fn process_sequence(&mut self, record: &OwnedRecord) {
    self.total_reads += 1;
    self.total_bases += record.seq().len();

    for base in record.seq() {
      match char::from(base.clone()).to_uppercase().to_string().as_str() {
        "T" => self.t_count += 1,
        "C" => self.c_count += 1,
        "G" => self.g_count += 1,
        "A" => self.a_count += 1,
        "N" => self.n_count += 1,
        _ => {}
      }
    }
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use fastq::parse_path;

  #[test]
  fn test_process_sequence() {
    let fastq_path = "/Users/choppy/Downloads/test.fq.gz";
    let mut qc = FastQC::new();
    parse_path(Some(fastq_path), |parser| {
      parser
        .each(|record| {
          qc.process_sequence(&record.to_owned_record());
          return true;
        })
        .unwrap()
    })
    .unwrap();

    println!("{:?}", qc);
  }
}
