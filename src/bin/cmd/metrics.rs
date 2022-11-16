use log::*;
use preqc_pack::qc;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::sync::Arc;
use structopt::StructOpt;

const PATTERN_FILE: &[u8] = include_bytes!("../../../data/patterns.json");
const ADAPTER_LIST: &[u8] = include_bytes!("../../../data/adapter_list.txt");
const CONTAMINANT_LIST: &[u8] = include_bytes!("../../../data/contaminant_list.txt");

/// A collection of qc metrics, such as file size, md5sum, fastqc, checkmate
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

    /// The number of green threads.
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

pub fn run(args: &Arguments) {
    if Path::new(&args.output).is_dir() || &args.output == "" {
        for input in &args.input {
            run_with_args(input, args);
        }
    } else {
        error!("The output ({:?}) need to be a directory.", &args.output);
    }
}

pub fn run_with_args(input: &str, args: &Arguments) {
    let output = if Path::new(input).exists() {
        // TODO: Multi threads?
        if args.which == "checksum" {
            info!("Run checksum...");
            let md5sum = qc::hasher::checksum(input, &args.algorithm);
            format!("{}", serde_json::to_string(&md5sum).unwrap())
        } else {
            info!("Started reading patternfile");
            let (patterns, indexes, count) = if args.pattern_file.len() > 0 {
                qc::mislabeling::VAFMatrix::read_patterns(&args.pattern_file[..])
            } else {
                qc::mislabeling::VAFMatrix::read_patterns_with_reader(PATTERN_FILE)
            };

            info!("Finished reading patternfile");

            let mut count_vec: Vec<Option<usize>> = vec![None; count];
            for i in indexes {
                count_vec[i] = Some(0);
            }

            info!("Finished building count array");

            info!("Started reading contaminants file");
            let contaminants = if args.contaminant_file.len() > 0 {
                qc::fastqc::OverRepresentedSeqs::read_contaminants_file(&args.contaminant_file[..])
            } else {
                qc::fastqc::OverRepresentedSeqs::read_contaminants_list(CONTAMINANT_LIST)
            };
            info!("Finished reading contaminants file");

            info!("Started reading adapter file");
            let adapters = if args.adapter_file.len() > 0 {
                qc::fastqc::AdapterContent::read_adapter_file(&args.adapter_file[..])
            } else {
                qc::fastqc::AdapterContent::read_adapter_list(ADAPTER_LIST)
            };
            info!("Finished reading adapter file");

            if args.which == "fastqc" {
                info!("Run fastqc...");
            } else if args.which == "checkmate" {
                info!("Run checkmate...")
            } else {
                info!("Run checksum, fastqc and checkmate...");
            }

            let overrepresented_max_unique_seq_count = if args.overrepresented_musc == 0 {
                Some(2 ^ 64 - 1)
            } else {
                Some(args.overrepresented_musc)
            };

            let kmer_ignore_sampling_interval = if args.kmer_isi == 0 {
                Some(2 ^ 64 - 1)
            } else {
                Some(args.kmer_isi)
            };

            let tile_continuous_sampling_boundary = if args.tile_csb == 0 {
                Some(2 ^ 64 - 1)
            } else {
                Some(args.tile_csb)
            };

            let tile_ignore_sampling_interval = if args.tile_isi == 0 {
                Some(2 ^ 64 - 1)
            } else {
                Some(args.tile_isi)
            };

            let fastqc_config = qc::FastQCConfig::new(
                adapters,
                contaminants,
                overrepresented_max_unique_seq_count,
                kmer_ignore_sampling_interval,
                tile_continuous_sampling_boundary,
                tile_ignore_sampling_interval,
            );

            let mislabeling_config = qc::MislabelingConfig::new(patterns, count_vec, count);

            let mut qc = if args.nthreads == 1 {
                qc::QCResults::run_qc(
                    input,
                    &args.which,
                    &fastqc_config,
                    &mislabeling_config,
                )
            } else {
                info!("Run with {:?} threads", args.nthreads);
                let which = Arc::new(args.which.to_string());

                let fastqc_config_arc = Arc::new(fastqc_config);
                let mislabeling_config_arc = Arc::new(mislabeling_config);

                qc::QCResults::run_qc_par(
                    input,
                    args.nthreads,
                    which,
                    fastqc_config_arc,
                    mislabeling_config_arc,
                )
            };

            if args.which != "all" {
                format!("{}", serde_json::to_string(&qc).unwrap())
            } else {
                qc.set_filemeta(Some(qc::hasher::checksum(input, &args.algorithm)));
                format!("{}", serde_json::to_string(&qc).unwrap())
            }
        }
    } else {
        error!("{} - Not Found: {:?}", module_path!(), args.input);
        std::process::exit(1);
    };

    // xxx.fq.gz/xxx.fastq.gz -> xxx
    // xxx.fq/xxx.fastq -> xxx
    let basename = Path::new(input).file_stem().unwrap();
    let basename = Path::new(basename).file_stem().unwrap();

    let filepath = if args.output.len() > 0 {
        Path::new(&args.output).join(format!("{}.json", basename.to_str().unwrap()))
    } else {
        Path::new(".").join(format!("{}.json", basename.to_str().unwrap()))
    };

    let mut f = File::create(filepath).unwrap();
    f.write(output.as_bytes()).unwrap();
}
