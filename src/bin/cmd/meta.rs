use log::*;
use preqc_pack::{qc, util};
use std::fs::File;
use std::io::Write;
use std::path::Path;
use structopt::StructOpt;

/// A collection of metadata, such as file size, md5sum
#[derive(StructOpt, PartialEq, Debug)]
#[structopt(setting=structopt::clap::AppSettings::ColoredHelp, name="PreQC Tool Suite - Hasher", author="Jingcheng Yang <yjcyxky@163.com>")]
pub struct Arguments {
  /// Fastq file to process
  #[structopt(name = "FILE")]
  input: String,

  /// SNP pattern file (format: BSON).
  #[structopt(name = "pattern-file", short="p", long="pattern-file")]
  pattern_file: String,

  /// A hash algorithms for output file.
  #[structopt(name="algorithm", short="m", long="algorithm", possible_values=&["md5sum", "blake2b"], default_value="md5sum")]
  algorithm: String,

  /// Which region?
  #[structopt(name="region", long="region", possible_values=&["cn-shanghai"], default_value="cn-shanghai")]
  region: String,

  /// Use internal network to get remote file.
  #[structopt(short = "i", long = "internal")]
  internal: bool,

  /// Which module will be called.
  #[structopt(name="which", short="w", long="which", possible_values=&["checksum", "fastqc", "all"], default_value="all")]
  which: String,

  /// The number of green threads.
  #[structopt(
    name = "nthreads",
    short = "n",
    long = "nthreads",
    default_value = "12"
  )]
  nthreads: u64,

  /// Chunk size (Bytes).
  #[structopt(
    name = "chunksize",
    short = "c",
    long = "chunksize",
    default_value = "8388608"
  )]
  chunk_size: u64,

  /// Output file.
  #[structopt(name = "output", short = "o", long = "output", default_value = "")]
  output: String,
}

pub async fn run(args: &Arguments) {
  let md5sum: qc::hasher::Meta;
  let mut output = String::from("");

  if util::is_remote_file(&args.input) {
    md5sum = qc::hasher::checksum_remote(
      &args.input,
      &args.algorithm,
      &args.region,
      args.internal,
      args.nthreads,
      args.chunk_size,
    )
    .await;
    output = format!("{}", serde_json::to_string(&md5sum).unwrap());
  } else {
    if Path::new(&args.input).exists() {
      // TODO: Multi threads?
      if args.which == "checksum" {
        let md5sum = qc::hasher::checksum(&args.input, &args.algorithm);
        output = format!("{}", serde_json::to_string(&md5sum).unwrap());
      } else if args.which == "fastqc" {
        let qc = qc::QCResults::run_fastqc(&args.input, &args.pattern_file);
        output = format!("{}", serde_json::to_string(&qc).unwrap());
      } else {
        let mut qc_results = qc::QCResults::run_fastqc(&args.input, &args.pattern_file);
        qc_results.set_filemeta(Some(qc::hasher::checksum(&args.input, &args.algorithm)));
        output = format!("{}", serde_json::to_string(&qc_results).unwrap());
      }
    } else {
      error!("{} - Not Found: {:?}", module_path!(), args.input);
    }
  }

  if args.output.len() > 0 {
    let mut f = File::create(&args.output).unwrap();
    f.write(output.as_bytes()).unwrap();
  } else {
    print!("{}", output);
  }
}
