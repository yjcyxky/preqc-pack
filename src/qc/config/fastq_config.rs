use crate::qc::fastq::{fastqc, hasher, mislabeling};
use serde::{Deserialize, Serialize};

use fastq::parse_path;
use std::collections::HashMap;
// use hashbrown::HashMap;
use std::io::Error;
use std::path::Path;
use std::sync::Arc;

#[allow(dead_code)]
pub struct FastQCConstants {}
impl FastQCConstants {
    pub const DEFAULT_NUMBER_OF_WINDOWS: i32 = 400;
    pub const DEFAULT_CHUNK_SIZE: i32 = 1000;

    pub const SANGER_ENCODING_OFFSET: usize = 33;
    pub const ILLUMINA_1_3_ENCODING_OFFSET: usize = 64;
    pub const SANGER_ILLUMINA_1_9: &str = "Sanger / Illumina 1.9";
    pub const ILLUMINA_1_3: &str = "Illumina 1.3";
    pub const ILLUMINA_1_5: &str = "Illumina 1.5";

    pub const FORWARD_TYPE: usize = 1;
    pub const REVERSE_TYPE: usize = 2;

    pub const FASTQC_CONFIG_DUP_LENGTH: usize = 0;
    pub const FASTQC_CONFIG_KMER_SIZE: usize = 0;

    pub const INDICATOR_CONFIG_TILE_IGNORE: usize = 0;
    pub const INDICATOR_CONFIG_OVERREPESENTED_WARN: f64 = 0.1;

    pub const DUPLICATION_WARNING: f64 = 70.0;
    pub const DUPLICATION_ERROR: f64 = 50.0;

    pub const KMER_CONTENT_WARNING: f64 = 2.0;
    pub const KMER_CONTENT_ERROR: f64 = 5.0;

    pub const N_CONTENT_WARNING: f64 = 5.0;
    pub const N_CONTENT_ERROR: f64 = 20.0;

    pub const OVEREPRESENTED_WARNING: f64 = 0.1;
    pub const OVEREPRESENTED_ERROR: f64 = 1.0;

    pub const QUALITY_BASE_LOWER_WARNING: f64 = 10.0;
    pub const QUALITY_BASE_LOWER_ERROR: f64 = 5.0;
    pub const QUALITY_BASE_MEDIAN_WARNING: f64 = 25.0;
    pub const QUALITY_BASE_MEDIAN_ERROR: f64 = 20.0;

    pub const SEQUENCE_WARNING: f64 = 10.0;
    pub const SEQUENCE_ERROR: f64 = 20.0;

    pub const GC_SEQUENCE_WARNING: f64 = 15.0;
    pub const GC_SEQUENCE_ERROR: f64 = 30.0;

    pub const QUALITY_SEQUENCE_WARNING: usize = 20;
    pub const QUALITY_SEQUENCE_ERROR: usize = 27;

    pub const TILE_WARNING: f64 = 5.0;
    pub const TILE_ERROR: f64 = 10.0;

    pub const SEQUENCE_LENGTH_WARNING: usize = 1;
    pub const SEQUENCE_LENGTH_ERROR: usize = 1;

    pub const ADAPTER_CONTENT_WARNING: f64 = 5.0;
    pub const ADAPTER_CONTENT_ERROR: f64 = 10.0;

    pub const DUPLICATION_IGNORE: usize = 0;
    pub const KMER_IGNORE: usize = 0;
    pub const N_CONTENT_IGNORE: usize = 0;
    pub const OVEREPRESENTED_IGNORE: usize = 0;
    pub const QUALITY_BASE_IGNORE: usize = 0;
    pub const SEQUENCE_IGNORE: usize = 0;
    pub const GC_QUALITY_IGNORE: usize = 0;
    pub const QUALITY_SEQUENCE_IGNORE: usize = 0;
    pub const TILE_IGNORE: usize = 0;
    pub const SEQUENCE_LENGTH_IGNORE: usize = 0;
    pub const ADAPTER_IGNORE: usize = 0;
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct FastQCConfig {
    overrepresented_max_unique_seq_count: Option<usize>,
    kmer_ignore_smapling_interval: Option<usize>,
    tile_continuous_sampling_boundary: Option<usize>,
    tile_ignore_smapling_interval: Option<usize>,
    adapters: String,
    contaminants: String,
}

impl FastQCConfig {
    pub fn new(
        adapters: String,
        contaminants: String,
        overrepresented_max_unique_seq_count: Option<usize>,
        kmer_ignore_smapling_interval: Option<usize>,
        tile_continuous_sampling_boundary: Option<usize>,
        tile_ignore_smapling_interval: Option<usize>,
    ) -> FastQCConfig {
        return FastQCConfig {
            overrepresented_max_unique_seq_count,
            kmer_ignore_smapling_interval,
            tile_continuous_sampling_boundary,
            tile_ignore_smapling_interval,
            adapters,
            contaminants,
        };
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct MislabelingConfig {
    patterns: HashMap<String, [usize; 2]>,
    count_vec: Vec<Option<usize>>,
    count: usize,
}

impl MislabelingConfig {
    pub fn new(
        patterns: HashMap<String, [usize; 2]>,
        count_vec: Vec<Option<usize>>,
        count: usize,
    ) -> MislabelingConfig {
        return MislabelingConfig {
            patterns,
            count_vec,
            count,
        };
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct QCResults {
    filemeta: Option<hasher::Meta>,
    fastqc: Option<fastqc::FastQC>,
    vaf_matrix: Option<mislabeling::VAFMatrix>,
}

impl QCResults {
    pub fn fastqc(&self) -> &Option<fastqc::FastQC> {
        return &self.fastqc;
    }

    pub fn vaf_matrix(&self) -> &Option<mislabeling::VAFMatrix> {
        return &self.vaf_matrix;
    }

    pub fn set_filemeta(&mut self, filemeta: Option<hasher::Meta>) {
        self.filemeta = filemeta;
    }

    pub fn run_qc_par(
        fastq_path: &str,
        n_threads: usize,
        which: Arc<String>,
        fastqc_config: Arc<FastQCConfig>,
        mislabeling_config: Arc<MislabelingConfig>,
    ) -> QCResults {
        match parse_path(Some(fastq_path), |parser| {
            let which_arc = Arc::clone(&which);
            let result: Result<Vec<_>, Error> =
                parser.parallel_each(n_threads, move |record_sets| {
                    let which_step = &which_arc[..];
                    let mut qc = fastqc::FastQC::new(
                        &fastqc_config.contaminants,
                        &fastqc_config.adapters,
                        fastqc_config.overrepresented_max_unique_seq_count,
                        fastqc_config.kmer_ignore_smapling_interval,
                        fastqc_config.tile_continuous_sampling_boundary,
                        fastqc_config.tile_ignore_smapling_interval,
                    );

                    let mut vaf_matrix = mislabeling::VAFMatrix::new(
                        mislabeling_config.count,
                        &mislabeling_config.count_vec,
                    );

                    for record_set in record_sets {
                        for record in record_set.iter() {
                            if which_step == "fastqc" || which_step == "all" {
                                qc.process_sequence(&record);
                            }

                            if which_step == "checkmate" || which_step == "all" {
                                vaf_matrix
                                    .process_sequence_unsafe(&mislabeling_config.patterns, &record);
                            }
                        }
                    }

                    QCResults {
                        filemeta: None,
                        fastqc: Some(qc),
                        vaf_matrix: Some(vaf_matrix),
                    }
                });

            return match result {
                Ok(qc_results) => {
                    let which_step = &which[..];
                    let mut merged_qc = qc_results[0].fastqc().to_owned().unwrap();
                    let mut merged_vaf_matrix = qc_results[0].vaf_matrix().to_owned().unwrap();
                    for i in 1..qc_results.len() {
                        if which_step == "fastqc" || which_step == "all" {
                            if let Some(fastqc) = qc_results[i].fastqc().to_owned() {
                                merged_qc.merge(&[fastqc]);
                            }
                        }

                        if which_step == "checkmate" || which_step == "all" {
                            if let Some(vaf_matrix) = qc_results[i].vaf_matrix().to_owned() {
                                merged_vaf_matrix.merge(&[vaf_matrix]);
                            }
                        }
                    }

                    let fastqc = if which_step == "fastqc" || which_step == "all" {
                        merged_qc.finish();
                        let filename = Path::new(fastq_path).file_name().unwrap().to_str().unwrap();
                        Some(merged_qc.update_name(filename))
                    } else {
                        None
                    };

                    let vaf_matrix = if which_step == "checkmate" || which_step == "all" {
                        merged_vaf_matrix.finish();
                        Some(merged_vaf_matrix)
                    } else {
                        None
                    };

                    QCResults {
                        filemeta: None,
                        fastqc: fastqc,
                        vaf_matrix: vaf_matrix,
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
        which: &str,
        fastqc_config: &FastQCConfig,
        mislabeling_config: &MislabelingConfig,
    ) -> QCResults {
        match parse_path(Some(fastq_path), |parser| {
            let mut qc = fastqc::FastQC::new(
                &fastqc_config.contaminants,
                &fastqc_config.adapters,
                fastqc_config.overrepresented_max_unique_seq_count,
                fastqc_config.kmer_ignore_smapling_interval,
                fastqc_config.tile_continuous_sampling_boundary,
                fastqc_config.tile_ignore_smapling_interval,
            );
            let mut vaf_matrix = mislabeling::VAFMatrix::new(
                mislabeling_config.count,
                &mislabeling_config.count_vec,
            );
            parser
                .each(|record| {
                    if which == "fastqc" || which == "all" {
                        qc.process_sequence(&record);
                    }

                    if which == "checkmate" || which == "all" {
                        vaf_matrix.process_sequence_unsafe(&mislabeling_config.patterns, &record);
                    }

                    return true;
                })
                .expect("Invalid fastq file");

            let fastqc = if which == "fastqc" || which == "all" {
                let filename = Path::new(fastq_path).file_name().unwrap().to_str().unwrap();
                qc.finish();
                Some(qc.update_name(filename))
            } else {
                None
            };

            let vaf_matrix = if which == "checkmate" || which == "all" {
                vaf_matrix.finish();
                Some(vaf_matrix)
            } else {
                None
            };

            QCResults {
                filemeta: None,
                fastqc: fastqc,
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
