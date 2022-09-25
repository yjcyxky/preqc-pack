pub mod fastqc;
pub mod hasher;
pub mod mislabeling;
pub mod util;

use serde::{Deserialize, Serialize};

use fastq::parse_path;
// use serde_json;
use std::io::Error;
use std::path::Path;

const PATTERN_FILE: &[u8] = include_bytes!("../../data/patterns.bson");

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct QCResults {
    filemeta: Option<hasher::Meta>,
    fastqc: fastqc::FastQC,
    vaf_matrix: mislabeling::VAFMatrix,
}

impl QCResults {
    pub fn fastqc(&self) -> &fastqc::FastQC {
        return &self.fastqc;
    }

    pub fn vaf_matrix(&self) -> &mislabeling::VAFMatrix {
        return &self.vaf_matrix;
    }

    pub fn set_filemeta(&mut self, filemeta: Option<hasher::Meta>) {
        self.filemeta = filemeta;
    }

    pub fn run_fastqc(fastq_path: &str, pattern_path: &str, n_threads: usize) -> QCResults {
        let (patterns, indexes, count) = if pattern_path.len() > 0 {
            mislabeling::VAFMatrix::read_patterns(pattern_path)
        } else {
            mislabeling::VAFMatrix::read_patterns_with_reader(PATTERN_FILE)
        };

        match parse_path(Some(fastq_path), |parser| {
            let result: Result<Vec<_>, Error> =
                parser.parallel_each(n_threads, move |record_sets| {
                    let mut qc = fastqc::FastQC::new();
                    let mut vaf_matrix =
                        mislabeling::VAFMatrix::new(patterns.clone(), indexes.clone(), count);
                    for record_set in record_sets {
                        for record in record_set.iter() {
                            qc.process_sequence(&record.to_owned_record());
                            vaf_matrix.process_sequence(&record.to_owned_record());
                        }
                    }

                    QCResults {
                        filemeta: None,
                        fastqc: qc,
                        vaf_matrix: vaf_matrix,
                    }
                });

            return match result {
                Ok(qc_results) => {
                    let mut merged_qc = qc_results[0].fastqc().to_owned();
                    let mut merged_vaf_matrix = qc_results[0].vaf_matrix().to_owned();
                    for i in 1..qc_results.len() {
                        merged_qc.merge(&[qc_results[i].fastqc().to_owned()]);
                    }

                    for i in 1..qc_results.len() {
                        merged_vaf_matrix.merge(&[qc_results[i].vaf_matrix().to_owned()]);
                    }

                    let filename = Path::new(fastq_path).file_name().unwrap().to_str().unwrap();
                    QCResults {
                        filemeta: None,
                        fastqc: merged_qc.update_name(filename),
                        vaf_matrix: merged_vaf_matrix,
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
}
