use log::*;
use preqc_pack::qc::fastq::{fastqc,hasher,mislabeling};
use preqc_pack::qc::config::fastq_config::{ FastQCConfig, MislabelingConfig, QCResults};
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::sync::Arc;
use serde::{Deserialize, Serialize};
use structopt::StructOpt;
use std::thread;

const PATTERN_FILE: &[u8] = include_bytes!("../../../../data/patterns.json");
const ADAPTER_LIST: &[u8] = include_bytes!("../../../../data/adapter_list.txt");
const CONTAMINANT_LIST: &[u8] = include_bytes!("../../../../data/contaminant_list.txt");

/// A collection of fastq qc metrics, such as file size, md5sum, fastqc, checkmate.Run fastq -h for more usage.
#[derive(StructOpt, PartialEq, Debug)]
#[structopt(
    setting=structopt::clap::AppSettings::ColoredHelp, 
    name="PreQC Tool Suite", 
    author="Jingcheng Yang <yjcyxky@163.com>; Haonan Chen <haonanchen0815@163.com>"
)]

pub struct Arguments {
    /// Fastq file(s) to process
    #[structopt(name = "FILE")]
    input: Vec<String>,

    /// Which module will be called.
    #[structopt(name="which", short="w", long="which", possible_values=&["checksum", "fastqc", "checkmate", "all"], default_value="all")]
    which: String,

    /// The number of green threads for each file.
    #[structopt(name = "nthreads", short = "n", long = "nthreads", default_value = "1")]
    nthreads: usize,

    /// Output directory.
    #[structopt(name = "output", short = "o", long = "output", default_value = "")]
    output: String,

    /// [checksum] A hash algorithms.
    #[structopt(name="algorithm", short="m", long="algorithm", possible_values=&["md5sum", "blake2b"], default_value="md5sum")]
    algorithm: String,

    /// [checkmate] SNP pattern file (format: JSON).
    #[structopt(
        name = "pattern-file",
        short = "p",
        long = "pattern-file",
        default_value = ""
    )]
    pattern_file: String,

    /// [fastqc] Adapter file (format: txt).
    #[structopt(
        name = "adapter-file",
        short = "A",
        long = "adapter-file",
        default_value = ""
    )]
    adapter_file: String,

    /// [fastqc] Contaminant file (format: txt).
    #[structopt(
        name = "contaminant-file",
        short = "C",
        long = "contaminant-file",
        default_value = ""
    )]
    contaminant_file: String,

    /// [fastqc] The number of different sequences we want to track in overrepresented module. it will be unlimit when you specify 0.
    #[structopt(
        name = "overrepresented-musc",
        long = "overrepresented-musc",
        default_value = "100000"
    )]
    overrepresented_musc: usize,

    /// [fastqc] Ignore one read out of  every a  specified number of reads in 'kmer' module. it will be unlimit when you specify 0.
    #[structopt(name = "kmer-isi", long = "kmer-isi", default_value = "50")]
    kmer_isi: usize,

    /// [fastqc] The max number of sequences we want to process continuously in 'tile quality' module. it will be unlimit when you specify 0.
    #[structopt(name = "tile-csb", long = "tile-csb", default_value = "10000")]
    tile_csb: usize,

    /// [fastqc] Ignore one read out of  every a  specified number of reads when crossing boundaries in 'tile quality' module. it will be unlimit when you specify 0.
    #[structopt(name = "tile-isi", long = "tile-isi", default_value = "10")]
    tile_isi: usize,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct MetricsConfig {
    which: String,
    algorithm: String,
    nthreads: usize,
    fastqc_config: FastQCConfig,
    mislabeling_config: MislabelingConfig
}

impl MetricsConfig {
    pub fn new(args: &Arguments ) -> MetricsConfig {
        let fastqc_config = MetricsConfig::construct_fastqc_config(&args.contaminant_file, &args.adapter_file, args.overrepresented_musc, args.kmer_isi, args.tile_csb, args.tile_isi);
        let mislabeling_config = MetricsConfig::construct_mislabeling_config(&args.pattern_file);

        MetricsConfig { nthreads:args.nthreads, which: args.which.clone(), algorithm: args.algorithm.clone(), fastqc_config, mislabeling_config }
    }

    pub fn construct_fastqc_config( contaminant_file: &str, adapter_file: &str, overrepresented_musc: usize, kmer_isi: usize, tile_csb: usize, tile_isi: usize) ->FastQCConfig {
        info!("Started reading contaminants file");
        let contaminants = if contaminant_file.len() > 0 {
            fastqc::OverRepresentedSeqs::read_contaminants_file(contaminant_file)
        } else {
            fastqc::OverRepresentedSeqs::read_contaminants_list(CONTAMINANT_LIST)
        };
        info!("Finished reading contaminants file");

        info!("Started reading adapter file");
        let adapters = if adapter_file.len() > 0 {
            fastqc::AdapterContent::read_adapter_file(adapter_file)
        } else {
            fastqc::AdapterContent::read_adapter_list(ADAPTER_LIST)
        };
        info!("Finished reading adapter file");

        let overrepresented_max_unique_seq_count = if overrepresented_musc == 0 {
            Some(2 ^ 64 - 1)
        } else {
            Some(overrepresented_musc)
        };

        let kmer_ignore_sampling_interval = if kmer_isi == 0 {
            Some(2 ^ 64 - 1)
        } else {
            Some(kmer_isi)
        };

        let tile_continuous_sampling_boundary = if tile_csb == 0 {
            Some(2 ^ 64 - 1)
        } else {
            Some(tile_csb)
        };

        let tile_ignore_sampling_interval = if tile_isi == 0 {
            Some(2 ^ 64 - 1)
        } else {
            Some(tile_isi)
        };

        FastQCConfig::new(
            adapters,
            contaminants,
            overrepresented_max_unique_seq_count,
            kmer_ignore_sampling_interval,
            tile_continuous_sampling_boundary,
            tile_ignore_sampling_interval,
        )
    }

    pub fn construct_mislabeling_config( pattern_file: &str) ->MislabelingConfig {
        info!("Started reading patternfile");
        let (patterns, indexes, count) = if pattern_file.len() > 0 {
            mislabeling::VAFMatrix::read_patterns(pattern_file)
        } else {
            mislabeling::VAFMatrix::read_patterns_with_reader(PATTERN_FILE)
        };

        info!("Finished reading patternfile");

        let mut count_vec: Vec<Option<usize>> = vec![None; count];
        for i in indexes {
            count_vec[i] = Some(0);
        }

        info!("Finished building count array");

        MislabelingConfig::new(patterns, count_vec, count)
    }
}

pub fn run(args: &Arguments) {
    if args.input.len() == 0 {
        error!("The input file should not be null.");
        return;
    }

    info!("Run with {:?} threads", args.nthreads);
    if Path::new(&args.output).is_dir() || &args.output == "" {
        let config = MetricsConfig::new(args);

        if args.input.len() > 1 {
            let inputs = args.input.to_owned();
            let mut handles = vec![];
            let output_arc = Arc::new(args.output.to_owned());
            let config_arc = Arc::new(config.clone());

            for input in inputs {
                let output_arc_ = output_arc.clone();
                let config_arc_ = config_arc.clone();
                handles.push(
                    thread::spawn(move|| {
                        run_with_args(&input, &output_arc_, &config_arc_);
                    })
                )
            }

            for handle in handles {
                handle.join().unwrap();
            }
        } else  {
            run_with_args(&args.input[0], &args.output, &config)
        } 
    } else {
        error!("The output ({:?}) need to be a directory.", &args.output);
    }
}

pub fn run_with_args(input: &str, output: &str, config: &MetricsConfig) {


    let results = if Path::new(input).exists() {
        // TODO: Multi threads?
        if config.which == "checksum" {
            info!("Run checksum on {:?}...", input);
            let md5sum = hasher::checksum(input, &config.algorithm);
            format!("{}", serde_json::to_string(&md5sum).unwrap())
        } else {
            if config.which == "fastqc" {
                info!("Run fastqc on {:?}...", input);
            } else if config.which == "checkmate" {
                info!("Run checkmate on {:?}...", input)
            } else {
                info!("Run checksum, fastqc and checkmate on {:?}...", input);
            }

            let mut qc = if config.nthreads == 1 {
                QCResults::run_qc(
                    input,
                    &config.which,
                    &config.fastqc_config,
                    &config.mislabeling_config,
                )
            } else {
                let which = Arc::new(config.which.clone());

                let fastqc_config_arc = Arc::new(config.fastqc_config.clone());
                let mislabeling_config_arc = Arc::new(config.mislabeling_config.clone());

                QCResults::run_qc_par(
                    input,
                    config.nthreads,
                    which,
                    fastqc_config_arc,
                    mislabeling_config_arc,
                )
            };

            if config.which != "all" {
                format!("{}", serde_json::to_string(&qc).unwrap())
            } else {
                qc.set_filemeta(Some(hasher::checksum(input, &config.algorithm)));
                format!("{}", serde_json::to_string(&qc).unwrap())
            }
        }
    } else {
        error!("{} - Not Found: {:?}", module_path!(), input);
        std::process::exit(1);
    };

    // xxx.fq.gz/xxx.fastq.gz -> xxx
    // xxx.fq/xxx.fastq -> xxx
    let basename = Path::new(input).file_stem().unwrap();
    let basename = Path::new(basename).file_stem().unwrap();

    let filepath = if output.len() > 0 {
        Path::new(output).join(format!("{}.json", basename.to_str().unwrap()))
    } else {
        Path::new(".").join(format!("{}.json", basename.to_str().unwrap()))
    };

    let mut f = File::create(filepath).unwrap();
    f.write(results.as_bytes()).unwrap();
}
