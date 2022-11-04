use log::*;
use preqc_pack::qc;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::sync::Arc;
use structopt::StructOpt;

const PATTERN_FILE: &[u8] = include_bytes!("../../../data/patterns.bson");

/// A collection of metadata, such as file size, md5sum
#[derive(StructOpt, PartialEq, Debug)]
#[structopt(setting=structopt::clap::AppSettings::ColoredHelp, name="PreQC Tool Suite - Hasher", author="Jingcheng Yang <yjcyxky@163.com>")]
pub struct Arguments {
    /// Fastq file to process
    #[structopt(name = "FILE")]
    input: String,

    /// SNP pattern file (format: BSON).
    #[structopt(
        name = "pattern-file",
        short = "p",
        long = "pattern-file",
        default_value = ""
    )]
    pattern_file: String,

    /// A hash algorithms for output file.
    #[structopt(name="algorithm", short="m", long="algorithm", possible_values=&["md5sum", "blake2b"], default_value="md5sum")]
    algorithm: String,

    /// Which module will be called.
    #[structopt(name="which", short="w", long="which", possible_values=&["checksum", "fastqc", "all"], default_value="all")]
    which: String,

    /// The number of green threads.
    #[structopt(name = "nthreads", short = "n", long = "nthreads", default_value = "1")]
    nthreads: usize,

    /// Output file.
    #[structopt(name = "output", short = "o", long = "output", default_value = "")]
    output: String,
}

pub fn run(args: &Arguments) {
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

    let output = if Path::new(&args.input).exists() {
        // TODO: Multi threads?
        if args.which == "checksum" {
            eprintln!("Run checksum...");
            let md5sum = qc::hasher::checksum(&args.input, &args.algorithm);
            format!("{}", serde_json::to_string(&md5sum).unwrap())
        } else {
            if args.which == "fastqc" {
                eprintln!("Run fastqc...");
            } else {
                eprintln!("Run checksum and fastqc...");
            }

            let mut qc = if args.nthreads == 1 {
                qc::QCResults::run_fastqc(&args.input, &patterns, &count_vec, count)
            } else {
                eprintln!("Run with {:?} threads", args.nthreads);
                let patterns = Arc::new(patterns);
                let count_vec = Arc::new(count_vec);
                qc::QCResults::run_fastqc_par(
                    &args.input,
                    patterns,
                    count_vec,
                    count,
                    args.nthreads,
                )
            };
            if args.which == "fastqc" {
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
