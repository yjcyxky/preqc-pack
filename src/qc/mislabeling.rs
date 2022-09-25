use bson::Document;
use digest::{Digest, Output};
use fastq::{OwnedRecord, Record};
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::Read;

const PATTERN_LENGTH: usize = 21;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct VAFMatrix {
    #[serde(skip_serializing)]
    length: usize,
    indexes: Vec<usize>,
    reference: Vec<Option<usize>>,
    alternative: Vec<Option<usize>>,
    #[serde(skip_serializing)]
    seq_ref_hited: Vec<usize>,
    #[serde(skip_serializing)]
    seq_alt_hited: Vec<usize>,
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
    /// let (patterns, indexes, count) = VAFMatrix::read_patterns("data/patterns.bson");
    /// assert_eq!(&Bson::Array(vec![Bson::Int32(0), Bson::Int32(0)]), patterns.get("TCCTTGTCATATGTTTTTCTG").unwrap());
    /// ```
    ///
    pub fn read_patterns(pattern_file: &str) -> (Document, Vec<usize>, usize) {
        let f = match File::open(pattern_file) {
            Ok(f) => f,
            Err(msg) => panic!("Cannot open {} - {}", pattern_file, msg),
        };

        return VAFMatrix::read_patterns_with_reader(f);
    }

    pub fn read_patterns_with_reader<R: Read>(mut reader: R) -> (Document, Vec<usize>, usize) {
        return match Document::from_reader(&mut reader) {
            Ok(fcontent) => {
                let count = fcontent.get("count").unwrap().as_i32().unwrap();
                let data = fcontent.get("data").unwrap();
                let indexes = fcontent.get("indexes").unwrap().as_array().unwrap();
                let patterns = data.as_document().unwrap();
                let indexes = indexes
                    .into_iter()
                    .map(|i| i.as_i32().unwrap() as usize)
                    .collect::<Vec<usize>>();

                (patterns.to_owned(), indexes.to_owned(), count as usize)
            }
            Err(msg) => {
                panic!("Cannot read pattern file {}", msg)
            }
        };
    }

    pub fn new(indexes: &Vec<usize>, count: usize) -> VAFMatrix {
        let mut init_values: Vec<Option<usize>> = vec![None; count];
        for i in indexes {
            init_values[*i] = Some(0);
        }

        return VAFMatrix {
            length: count,
            indexes: (0..count).into_iter().collect(),
            reference: init_values.clone(),
            alternative: init_values.clone(),
            seq_alt_hited: vec![],
            seq_ref_hited: vec![],
        };
    }

    pub fn merge(&mut self, vaf_matrixes: &[VAFMatrix]) {
        for i in vaf_matrixes {
            for j in 0..self.alternative.len() {
                if self.alternative[j] == None && i.alternative[j] == None {
                    self.alternative[j] = None;
                } else if self.alternative[j] == None {
                    self.alternative[j] = i.alternative[j];
                } else {
                    self.alternative[j] =
                        Some(self.alternative[j].unwrap() + i.alternative[j].unwrap());
                }
            }

            for k in 0..self.reference.len() {
                if self.reference[k] == None && i.reference[k] == None {
                    self.reference[k] = None;
                } else if self.reference[k] == None {
                    self.reference[k] = i.reference[k];
                } else {
                    self.reference[k] = Some(self.reference[k].unwrap() + i.reference[k].unwrap());
                }
            }
        }
    }

    pub fn from_pattern_file(pattern_file: &str) -> (VAFMatrix, Document) {
        let (patterns, indexes, count) = VAFMatrix::read_patterns(pattern_file);
        let mut init_values = vec![None; count];
        for i in indexes {
            init_values[i] = Some(0);
        }

        return (
            VAFMatrix {
                length: count,
                indexes: (0..count).into_iter().collect(),
                reference: init_values.clone(),
                alternative: init_values.clone(),
                seq_alt_hited: vec![],
                seq_ref_hited: vec![],
            },
            patterns,
        );
    }

    fn reset_seq_hited(&mut self) {
        self.seq_alt_hited = vec![];
        self.seq_ref_hited = vec![];
    }

    fn exists_in_seq_ref_hited(&self, index: usize) -> bool {
        match self.seq_ref_hited.iter().position(|&r| r == index) {
            None => false,
            _ => true,
        }
    }

    fn exists_in_seq_alt_hited(&self, index: usize) -> bool {
        match self.seq_alt_hited.iter().position(|&r| r == index) {
            None => false,
            _ => true,
        }
    }

    fn hash<D: Digest + Default>(&self, seq: &[u8]) -> String
    where
        Output<D>: core::fmt::LowerHex,
    {
        let mut sh = D::new();
        sh.update(&seq);
        format!("{:x}", sh.finalize())
    }

    pub fn process_sequence_unsafe(&mut self, patterns: &Document, fastq_record: &OwnedRecord) {
        let seq = fastq_record.seq();
        let length = seq.len();

        self.reset_seq_hited();
        for i in 0..length {
            if PATTERN_LENGTH + i <= length {
                let substr = &seq[i..PATTERN_LENGTH + i];

                unsafe {
                    // It assume the sequence which is from the fastq file has valid ASCII and UTF8 characters.
                    let s = std::str::from_utf8_unchecked(substr);
                    match patterns.get(s) {
                        Some(matched) => {
                            let matched_array = matched.as_array().unwrap();
                            // matched_array: 0 - index; 1 - ref_or_alt;
                            let index = matched_array[0].as_i32().unwrap() as usize;
                            // For Debug
                            // println!(
                            //   "Seq: {}, Index: {}, Matched Hash: {}, RefAlt: {}",
                            //   self.hash::<Md5>(seq),
                            //   index,
                            //   s,
                            //   matched_array[1].as_i32().unwrap()
                            // );
                            match matched_array[1].as_i32().unwrap() {
                                // 0 = ref, 1 = alt
                                1 => {
                                    if !self.exists_in_seq_alt_hited(index) {
                                        self.alternative[index] = match self.alternative[index] {
                                            None => Some(1),
                                            Some(num) => Some(num + 1),
                                        };
                                        self.seq_alt_hited.push(index);
                                    }
                                }
                                0 => {
                                    if !self.exists_in_seq_ref_hited(index) {
                                        self.reference[index] = match self.reference[index] {
                                            None => Some(1),
                                            Some(num) => Some(num + 1),
                                        };
                                        self.seq_ref_hited.push(index);
                                    }
                                }
                                _ => {}
                            }
                        }

                        None => {
                            // println!("No Matched Hash: {}", s);
                        }
                    }
                }
            }
        }

        // Only need to record once even if it (the same SNP) is matched more than once in the same sequence.
        self.reset_seq_hited();
    }

    pub fn process_sequence(&mut self, patterns: &Document, fastq_record: &OwnedRecord) {
        let seq = fastq_record.seq();
        let length = seq.len();

        self.reset_seq_hited();
        for i in 0..length {
            if PATTERN_LENGTH + i <= length {
                let substr = &seq[i..PATTERN_LENGTH + i];
                match std::str::from_utf8(substr) {
                    Err(msg) => {
                        panic!("Cannot parse the sequence: {}", msg)
                    }
                    Ok(s) => match patterns.get(s) {
                        Some(matched) => {
                            let matched_array = matched.as_array().unwrap();
                            // matched_array: 0 - index; 1 - ref_or_alt;
                            let index = matched_array[0].as_i32().unwrap() as usize;
                            // For Debug
                            // println!(
                            //   "Seq: {}, Index: {}, Matched Hash: {}, RefAlt: {}",
                            //   self.hash::<Md5>(seq),
                            //   index,
                            //   s,
                            //   matched_array[1].as_i32().unwrap()
                            // );
                            match matched_array[1].as_i32().unwrap() {
                                // 0 = ref, 1 = alt
                                1 => {
                                    if !self.exists_in_seq_alt_hited(index) {
                                        self.alternative[index] = match self.alternative[index] {
                                            None => Some(1),
                                            Some(num) => Some(num + 1),
                                        };
                                        self.seq_alt_hited.push(index);
                                    }
                                }
                                0 => {
                                    if !self.exists_in_seq_ref_hited(index) {
                                        self.reference[index] = match self.reference[index] {
                                            None => Some(1),
                                            Some(num) => Some(num + 1),
                                        };
                                        self.seq_ref_hited.push(index);
                                    }
                                }
                                _ => {}
                            }
                        }

                        None => {
                            // println!("No Matched Hash: {}", s);
                        }
                    },
                }
            }
        }

        // Only need to record once even if it (the same SNP) is matched more than once in the same sequence.
        self.reset_seq_hited();
    }
}
