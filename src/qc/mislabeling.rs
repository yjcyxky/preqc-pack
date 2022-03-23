use bson::Document;
use fastq::{OwnedRecord, Record};
use serde::{Deserialize, Serialize};
use std::fs::File;

const PATTERN_LENGTH: usize = 21;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct VAFMatrix {
  #[serde(skip_serializing)]
  length: usize,
  #[serde(skip_serializing)]
  patterns: Document,
  reference: Vec<i32>,
  alternative: Vec<i32>,
}

impl VAFMatrix {
  /// Read SNP patterns from bson file.
  ///
  /// NOTE: BSON Document -> HashMap<String, [usize; 2]>
  ///
  /// # Arguments
  ///
  /// * `pattern_file`: where is the pattern file.
  ///
  /// # Example
  ///
  /// ```
  /// extern crate preqc_pack;
  /// use preqc_pack::qc::mislabeling::VAFMatrix;
  /// use bson::Bson;
  ///
  /// let (patterns, count) = VAFMatrix::read_patterns("data/patterns.bson");
  /// assert_eq!(&Bson::Array(vec![Bson::Int32(0), Bson::Int32(0)]), patterns.get("TCCTTGTCATATGTTTTTCTG").unwrap());
  /// ```
  ///
  pub fn read_patterns(pattern_file: &str) -> (Document, usize) {
    let mut f = match File::open(pattern_file) {
      Ok(f) => f,
      Err(msg) => panic!("Cannot open {} - {}", pattern_file, msg),
    };

    return match Document::from_reader(&mut f) {
      Ok(fcontent) => {
        let count = fcontent.get("count").unwrap().as_i32().unwrap();
        let data = fcontent.get("data").unwrap();
        let patterns = data.as_document().unwrap();
        (patterns.to_owned(), count as usize)
      }
      Err(msg) => {
        panic!("Cannot read pattern file {} - {}", pattern_file, msg)
      }
    };
  }

  pub fn new(patterns: Document, count: usize) -> VAFMatrix {
    return VAFMatrix {
      length: count,
      patterns: patterns,
      reference: vec![-1; count],
      alternative: vec![-1; count],
    };
  }

  pub fn from_pattern_file(pattern_file: &str) -> VAFMatrix {
    let (patterns, count) = VAFMatrix::read_patterns(pattern_file);
    return VAFMatrix {
      length: count,
      patterns: patterns,
      reference: vec![-1; count],
      alternative: vec![-1; count],
    };
  }

  pub fn process_sequence(&mut self, fastq_record: &OwnedRecord) {
    let seq = fastq_record.seq();
    let length = seq.len();

    for i in 0..length {
      if PATTERN_LENGTH + i <= length {
        let substr = &seq[i..PATTERN_LENGTH + i];
        match std::str::from_utf8(substr) {
          Err(msg) => {
            panic!("Cannot parse the sequence: {}", msg)
          }
          Ok(s) => match self.patterns.get(s) {
            Some(matched) => {
              let matched_array = matched.as_array().unwrap();
              // matched_array: 0 - index; 1 - ref_or_alt;
              let index = matched_array[0].as_i32().unwrap() as usize;
              match matched_array[1].as_i32().unwrap() {
                // 0 = ref, 1 = alt
                1 => {
                  if self.alternative[index] == -1 {
                    self.alternative[index] = 1;
                  } else {
                    self.alternative[index] += 1;
                  }
                }
                0 => {
                  if self.reference[index] == -1 {
                    self.reference[index] = 1;
                  } else {
                    self.reference[index] += 1;
                  }
                }
                _ => {}
              }
            }

            None => {}
          },
        }
      }
    }
  }
}
