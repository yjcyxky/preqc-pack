use log::*;
use preqc_pack::util;
use std::path::Path;
use structopt::StructOpt;

/// Merge several fastq files to one.
#[derive(StructOpt, PartialEq, Debug)]
#[structopt(setting=structopt::clap::AppSettings::ColoredHelp, name="PreQC Tool Suite - Merge", author="Jingcheng Yang <yjcyxky@163.com>")]
pub struct Arguments {
    /// Fastq files
    #[structopt(name = "FILE", multiple = true, takes_value = true)]
    inputs: Vec<String>,

    /// Output file.
    #[structopt(name = "output", short = "o", long = "output")]
    output: String,
}

fn exists(files: &[String]) -> bool {
    let items = files.iter().filter(|file| !Path::new(&file).exists());
    if items.count() > 0 {
        return false;
    }

    true
}

pub fn run(args: &Arguments) {
    let outputs: Vec<String> = vec![args.output.clone()];
    if !exists(&outputs) {
        if exists(&args.inputs) {
            // TODO: Multi threads?
            println!("Merge all intputs {:?} to {}", args.inputs, args.output);
            for input in args.inputs.clone() {
                util::zcat(&input, &args.output)
            }
        } else {
            error!("{} - Not Found: {:?}", module_path!(), args.inputs);
        }
    } else {
        error!("{} exists", args.output);
    }
}
