use log::*;
use preqc_pack::{qc, util};
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::thread;
use structopt::StructOpt;

/// A collection of metadata, such as file size, md5sum
#[derive(StructOpt, PartialEq, Debug)]
#[structopt(setting=structopt::clap::AppSettings::ColoredHelp, name="PreQC Tool Suite - Hasher", author="Jingcheng Yang <yjcyxky@163.com>")]
pub struct Arguments {
    /// FastQ R1 file to process (fastq.gz/fq.gz/fastq/fq)
    #[structopt(name = "FASTQ R1")]
    fastq_r1: String,

    /// FastQ R2 file to process (fastq.gz/fq.gz/fastq/fq)
    #[structopt(name = "FASTQ R2")]
    fastq_r2: String,

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
    #[structopt(name = "nthreads", short = "n", long = "nthreads", default_value = "5")]
    nthreads: u64,

    /// Output file.
    #[structopt(name = "output", short = "o", long = "output", default_value = "")]
    output: String,
}

pub fn run(
    input: String,
    pattern_file: String,
    algorithm: String,
    nthreads: u64,
    which: String,
) -> String {
    let output = if Path::new(&input[..]).exists() {
        // TODO: Multi threads?
        if which == "checksum" {
            eprintln!("Run checksum on {}...", &input[..]);
            let md5sum = qc::hasher::checksum(&input[..], &algorithm[..]);
            format!("{}", serde_json::to_string(&md5sum).unwrap())
        } else if which == "fastqc" {
            eprintln!("Run fastqc & NGSCheckMate on {}...", &input[..]);
            let qc = qc::QCResults::run_fastqc(&input[..], &pattern_file[..], nthreads as usize);
            format!("{}", serde_json::to_string(&qc).unwrap())
        } else {
            eprintln!("Run checksum, fastqc & NGSCheckMate on {}...", &input[..]);
            let mut qc_results =
                qc::QCResults::run_fastqc(&input[..], &pattern_file[..], nthreads as usize);
            qc_results.set_filemeta(Some(qc::hasher::checksum(&input[..], &algorithm[..])));
            format!("{}", serde_json::to_string(&qc_results).unwrap())
        }
    } else {
        error!("{} - Not Found: {:?}", module_path!(), input);
        std::process::exit(1);
    };

    format!("{}", output)
}

pub fn batch_run(args: &Arguments) {
    let mut outputs = vec![];
    let fastq_r1 = args.fastq_r1.clone();
    let fastq_r2 = args.fastq_r2.clone();
    let pattern_file = args.pattern_file.clone();
    let algorithm = args.algorithm.clone();
    let nthreads = args.nthreads.clone();
    let which = args.which.clone();
    let handler = thread::spawn(move || run(fastq_r1, pattern_file, algorithm, nthreads, which));

    outputs.push(run(
        fastq_r2,
        args.pattern_file.clone(),
        args.algorithm.clone(),
        args.nthreads.clone(),
        args.which.clone(),
    ));
    outputs.push(handler.join().unwrap());

    let outputs = format!("[{}, {}]", outputs[0], outputs[1]);
    if args.output.len() > 0 {
        let mut f = File::create(&args.output).unwrap();
        f.write(outputs.as_bytes()).unwrap();
    } else {
        print!("{}", outputs);
    }
}
