use log::*;
use preqc_pack::{fastqc, hasher, util, QCPack};
use s3::creds::Credentials;
use s3::region::Region;
use std::env;
use std::fs::{self, File};
use std::io::Write;
use std::path::Path;
use structopt::StructOpt;

use blake2::Blake2b;
use md5::Md5;

/// A collection of metadata, such as file size, md5sum
#[derive(StructOpt, PartialEq, Debug)]
#[structopt(setting=structopt::clap::AppSettings::ColoredHelp, name="PreQC Tool Suite - Hasher", author="Jingcheng Yang <yjcyxky@163.com>")]
pub struct Arguments {
  /// Bam file to process
  #[structopt(name = "FILE")]
  input: String,

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
  #[structopt(name = "nthreads", short = "n", long = "nthreads", default_value = "12")]
  nthreads: u64,

  /// Chunk size (Bytes).
  #[structopt(name = "chunksize", short = "c", long = "chunksize", default_value = "8388608")]
  chunk_size: u64,

  /// Output file.
  #[structopt(name = "output", short = "o", long = "output", default_value = "")]
  output: String,
}

async fn checksum_remote(
  input: &str,
  algorithm: &str,
  region: &str,
  internal: bool,
  n_threads: u64,
  chunk_size: u64
) -> hasher::Meta {
  // Get filemeta
  let endpoint = if internal {
    format!("https://oss-{}-internal.aliyuncs.com", region).to_owned()
  } else {
    format!("https://oss-{}.aliyuncs.com", region).to_owned()
  };

  let region = Region::Custom {
    region: region.to_owned(),
    endpoint: endpoint,
  };

  let credentials = Credentials::new(
    Some(&env::var("OSS_ACCESS_KEY").unwrap()),
    Some(&env::var("OSS_ACCESS_KEY_SECRET").unwrap()),
    None,
    None,
    None,
  )
  .unwrap();

  let meta = match algorithm {
    "blake2b" => hasher::process_remote::<Blake2b>(input, region, credentials, n_threads, chunk_size).await,
    _ => hasher::process_remote::<Md5>(input, region, credentials, n_threads, chunk_size).await,
  };
  
  meta
}

fn checksum(input: &str, algorithm: &str) -> hasher::Meta {
  let mut file = fs::File::open(input).unwrap();
  let meta = match algorithm {
    "blake2b" => hasher::process::<Blake2b, _>(&mut file),
    _ => hasher::process::<Md5, _>(&mut file),
  };
  meta
}

fn checksum_mmap(input: &str, algorithm: &str) -> hasher::Meta {
  let file = fs::File::open(input).unwrap();
  let meta = match algorithm {
    "blake2b" => hasher::process_mmap::<Blake2b>(&file),
    _ => hasher::process_mmap::<Md5>(&file),
  };
  meta
}

fn fastqc(input: &str) -> fastqc::FastQC {
  // Generate fastqc metrics
  let fastqc_metrics = fastqc::compute_data_size(input);

  fastqc_metrics
}

pub async fn run(args: &Arguments) {
  let md5sum: hasher::Meta;
  let mut output = String::from("");

  if util::is_remote_file(&args.input) {
    md5sum = checksum_remote(&args.input, &args.algorithm, &args.region, args.internal, args.nthreads, args.chunk_size).await;
    output = format!("{}", serde_json::to_string(&md5sum).unwrap());
  } else {
    if Path::new(&args.input).exists() {
      // TODO: Multi threads?
      if args.which == "checksum" {
        let md5sum = checksum(&args.input, &args.algorithm);
        output = format!("{}", serde_json::to_string(&md5sum).unwrap());
      } else if args.which == "fastqc" {
        output = format!("{}", serde_json::to_string(&fastqc(&args.input)).unwrap());
      } else {
        let mut qc_pack = QCPack::new();
        qc_pack.set_filemeta(checksum(&args.input, &args.algorithm));
        qc_pack.set_fastqc(fastqc(&args.input));
        output = format!("{}", serde_json::to_string(&qc_pack).unwrap());
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
