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

/// A collection of metadata, such as file size, md5sum, fastqc, ngscheckmate
#[derive(StructOpt, PartialEq, Debug)]
#[structopt(setting=structopt::clap::AppSettings::ColoredHelp, name="PreQC Tool Suite - Hasher", author="Jingcheng Yang <yjcyxky@163.com>")]
pub struct Arguments {
    /// Fastq file to process
    #[structopt(name = "FILE")]
    input: String,

    /// SNP pattern file (format: JSON).
    #[structopt(
        name = "pattern-file",
        short = "p",
        long = "pattern-file",
        default_value = ""
    )]
    pattern_file: String,

    /// Adapter file (format: txt).
    #[structopt(
        name = "adapter-file",
        short = "A",
        long = "adapter-file",
        default_value = ""
    )]
    adapter_file: String,

    /// Contaminant file (format: txt).
    #[structopt(
        name = "contaminant-file",
        short = "C",
        long = "contaminant-file",
        default_value = ""
    )]
    contaminant_file: String,

    /// A hash algorithms for output file.
    #[structopt(name="algorithm", short="m", long="algorithm", possible_values=&["md5sum", "blake2b"], default_value="md5sum")]
    algorithm: String,

    /// Which module will be called.
    #[structopt(name="which", short="w", long="which", possible_values=&["checksum", "fastqc", "checkmate", "all"], default_value="all")]
    which: String,

    /// The number of green threads.
    #[structopt(name = "nthreads", short = "n", long = "nthreads", default_value = "1")]
    nthreads: usize,

    /// Output file.
    #[structopt(name = "output", short = "o", long = "output", default_value = "")]
    output: String,
}

pub fn run(args: &Arguments) {
    let output = if Path::new(&args.input).exists() {
        // TODO: Multi threads?
        if args.which == "checksum" {
            eprintln!("Run checksum...");
            let md5sum = qc::hasher::checksum(&args.input, &args.algorithm);
            format!("{}", serde_json::to_string(&md5sum).unwrap())
        } else {
            eprintln!("Started reading patternfile");
            let (patterns, indexes, count) = if args.pattern_file.len() > 0 {
                qc::mislabeling::VAFMatrix::read_patterns(&args.pattern_file[..])
            } else {
                qc::mislabeling::VAFMatrix::read_patterns_with_reader(PATTERN_FILE)
            };
            
            eprintln!("Finished reading patternfile");

            let mut count_vec: Vec<Option<usize>> = vec![None; count];
            for i in indexes {
                count_vec[i] = Some(0);
            }
        
            eprintln!("Finished building count array");

            eprintln!("Started reading contaminants file");
            let contaminants = if args.contaminant_file.len() > 0 {
                qc::fastqc::OverRepresentedSeqs::read_contaminants_file(&args.contaminant_file[..])
            } else {
                qc::fastqc::OverRepresentedSeqs::read_contaminants_list(CONTAMINANT_LIST)
            };
            eprintln!("Finished reading contaminants file");

            eprintln!("Started reading adapter file");
            let adapters = if args.adapter_file.len() > 0 {
                qc::fastqc::AdapterContent::read_adapter_file(&args.adapter_file[..])
            } else {
                qc::fastqc::AdapterContent::read_adapter_list(ADAPTER_LIST)
            };
            eprintln!("Finished reading adapter file");

            if args.which == "fastqc" {
                eprintln!("Run fastqc...");
            } else if args.which == "checkmate" {
                eprintln!("Run checkmate...")
            } else {
                eprintln!("Run checksum, fastqc and checkmate...");
            }

            let adapters = Arc::new(adapters);
            let contaminants = Arc::new(contaminants);

            let mut qc = if args.nthreads == 1 {
                qc::QCResults::run_qc(
                    &args.input,
                    adapters,
                    contaminants,
                    &patterns,
                    &count_vec,
                    count,
                    &args.which,
                )
            } else {
                eprintln!("Run with {:?} threads", args.nthreads);
                let patterns = Arc::new(patterns);
                let count_vec = Arc::new(count_vec);
                let which = Arc::new(args.which.to_string());
                qc::QCResults::run_qc_par(
                    &args.input,
                    adapters,
                    contaminants,
                    patterns,
                    count_vec,
                    count,
                    args.nthreads,
                    which,
                )
            };

            if args.which != "all" {
                format!("{}", serde_json::to_string(&qc).unwrap())
            } else {
                qc.set_filemeta(Some(qc::hasher::checksum(&args.input, &args.algorithm)));
                format!("{}", serde_json::to_string(&qc).unwrap())
            }
        }
    } else {
        error!("{} - Not Found: {:?}", module_path!(), args.input);
        std::process::exit(1);
    };

    if args.output.len() > 0 {
        let mut f = File::create(&args.output).unwrap();
        f.write(output.as_bytes()).unwrap();
    } else {
        print!("{}", output);
    }
}
