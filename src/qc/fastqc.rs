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
      // NOTE: to_string is expensive than from_utf8
      match std::str::from_utf8(&[base.clone()]).unwrap() {
      // match char::from(base.clone()).to_uppercase().to_string().as_str() {
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
  use fastq::OwnedRecord;

  #[test]
  fn test_process_sequence() {
    let mut qc = FastQC::new();
    let read1 = OwnedRecord {
      head: b"some_name".to_vec(),
      seq: b"GTCGCACTGATCTGGGTTAGGCGCGGAGCCGAGGGTTGCACCATTTTTCATTATTGAATGCCAAGATA".to_vec(),
      qual: b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII".to_vec(),
      sep: None,
    };
    qc.process_sequence(&read1);

    assert_eq!(qc.total_bases, 68);
    assert_eq!(qc.total_reads, 1);
    assert_eq!(qc.g_count, 20);
    assert_eq!(qc.a_count, 15);
    assert_eq!(qc.c_count, 14);
    assert_eq!(qc.t_count, 19);
    assert_eq!(qc.n_count, 0);
  }
}
