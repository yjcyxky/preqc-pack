pub mod fastqc;
pub mod hasher;
pub mod mislabeling;
pub mod util;

use serde::{Deserialize, Serialize};

use fastq::parse_path;
// use std::collections::HashMap;
use hashbrown::HashMap;
use std::io::Error;
use std::path::Path;
use std::sync::Arc;

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

    pub fn run_qc_par(
        fastq_path: &str,
        adapters: Arc<String>,
        contaminants: Arc<String>,
        patterns: Arc<HashMap<String, [usize; 2]>>,
        count_vec: Arc<Vec<Option<usize>>>,
        count: usize,
        n_threads: usize,
        which: Arc<String>,
    ) -> QCResults {
        match parse_path(Some(fastq_path), |parser| {
            let result: Result<Vec<_>, Error> =
                parser.parallel_each(n_threads, move |record_sets| {
                    let mut qc = fastqc::FastQC::new(&contaminants, &adapters);
                    let mut vaf_matrix = mislabeling::VAFMatrix::new(count, &count_vec);
                    for record_set in record_sets {
                        for record in record_set.iter() {
                            let which_step = &which[..];
                            if which_step == "fastqc" || which_step == "all" {
                                qc.process_sequence(&record);
                            }
                            
                            if which_step == "checkmate" || which_step == "all" {
                                vaf_matrix.process_sequence_unsafe(&patterns, &record);
                            }
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
                        merged_vaf_matrix.merge(&[qc_results[i].vaf_matrix().to_owned()]);
                    }

                    merged_qc.finish();

                    let filename = Path::new(fastq_path).file_name().unwrap().to_str().unwrap();
                    let fastqc = merged_qc.update_name(filename);

                    merged_vaf_matrix.finish();

                    QCResults {
                        filemeta: None,
                        fastqc: fastqc,
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

    pub fn run_qc(
        fastq_path: &str,
        adapters: Arc<String>,
        contaminants: Arc<String>,
        patterns: &HashMap<String, [usize; 2]>,
        count_vec: &Vec<Option<usize>>,
        count: usize,
        which: &str,
    ) -> QCResults {
        match parse_path(Some(fastq_path), |parser| {
            let mut qc = fastqc::FastQC::new(&contaminants, &adapters);
            let mut vaf_matrix = mislabeling::VAFMatrix::new(count, &count_vec);
            parser
                .each(|record| {
                    if which == "fastqc" || which == "all" {
                        qc.process_sequence(&record);
                    } else if which == "checkmate" || which == "all" {
                        vaf_matrix.process_sequence_unsafe(&patterns, &record);
                    }
                    
                    return true;
                })
                .expect("Invalid fastq file");

            let filename = Path::new(fastq_path).file_name().unwrap().to_str().unwrap();

            qc.finish();
            vaf_matrix.finish();

            QCResults {
                filemeta: None,
                fastqc: qc.update_name(filename),
                vaf_matrix: vaf_matrix,
            }
        }) {
            Err(msg) => {
                panic!("Cannot parse fastq file: {}", msg);
            }
            Ok(o) => o,
        }
    }
}
