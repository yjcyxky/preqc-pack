use log::*;
use preqc_pack::qc;
use preqc_pack::qc::config::bam_config::{QualimapConfig,QCResults};
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::sync::Arc;
use serde::{Deserialize, Serialize};
use structopt::StructOpt;
use std::thread;


/// A collection of bam metrics, such as fastq screen, qualimap. Run bam -h for more usage.
#[derive(StructOpt, PartialEq, Debug)]
#[structopt(
    setting=structopt::clap::AppSettings::ColoredHelp, 
    name="PreQC Tool Suite", 
    author="Jingcheng Yang <yjcyxky@163.com>; Haonan Chen <haonanchen0815@163.com>"
)]

pub struct Arguments {
    /// Bam file(s) to process
    #[structopt(name = "FILE")]
    input: Vec<String>,

    /// Which module will be called.
    #[structopt(name="which", short="w", long="which", possible_values=&["fastqscreen", "qualimap",  "all"], default_value="all")]
    which: String,

    /// The number of green threads for each file.
    #[structopt(name = "nthreads", short = "n", long = "nthreads", default_value = "1")]
    nthreads: usize,

    /// Output directory.
    #[structopt(name = "output", short = "o", long = "output", default_value = "")]
    output: String,

}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct MetricsConfig {
    which: String,
    nthreads: usize,
    qualimap_config:QualimapConfig
}

impl MetricsConfig {
    pub fn new(_which: &str,_nthreads: usize) -> Self {
        Self {
            which:_which.to_string(),
            nthreads:_nthreads,
            qualimap_config:QualimapConfig::new(0, 0, 0)
        }
    }
}

pub fn run(args: &Arguments) {
    info!("Run with {:?} threads", args.nthreads);
    if Path::new(&args.output).is_dir() || &args.output == "" {
        let config = MetricsConfig::new(
            &args.which, 
            args.nthreads, 
        );

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
        } else {
            run_with_args(&args.input[0], &args.output, &config)
        }
    } else {
        error!("The output ({:?}) need to be a directory.", &args.output);
    }
}

pub fn run_with_args(input: &str, output: &str, config: &MetricsConfig) {
    let results = if Path::new(input).exists() {
        // TODO: Multi threads?
        if config.which == "qualimap" {
            info!("Run qualimap on {:?}...", input);
        } else if config.which == "fastqscreen" {
            info!("Run checkmate on {:?}...", input)
        } else {
            info!("Run qualimap, fastqscreen on {:?}...", input);
        }

        let mut qc = if config.nthreads == 1 {
            QCResults::run_qc(
                input,
                &config.which,
                &config.qualimap_config,
            )
        } else {
            // multi-threads
            QCResults::new()
        };

        format!("{}", serde_json::to_string(&qc).unwrap())
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
