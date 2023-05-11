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


/// A collection of bam qc metrics, such as fastq screen, qualimap. Run bam -h for more usage.
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
    #[structopt( short = "n", long = "nthreads", default_value = "1")]
    nthreads: usize,

    /// Output directory.
    #[structopt( short = "o", long = "output", default_value = "")]
    output: String,

    /// [qualimap] Number of sample windows across the genome.
    #[structopt(name = "nwindows", short = "s",long = "nw",  default_value = "400")]
    nwindows:usize,

    /// [qualimap] Number of reads analyzed in a chunk.
    #[structopt(name = "nchunks", short = "c",long = "nc", default_value = "1000")]
    nchunks:usize,

    /// [qualimap]  Feature file with regions of interest in GFF/GTF or BED format.
    #[structopt( short = "f",long = "gff",default_value ="")]
    gff:String,

    /// [qualimap]  Activate this option to collect statistics of overlapping paired-end reads.
    #[structopt( short = "a",long = "collect-overlap-pairs")]
    collect_overlap_flag:bool,

    /// [qualimap]  Activate this option to skip duplicate alignments from the analysis.
    #[structopt( short = "d",long = "skip-duplicated")]
    skip_dup_flag: bool,

    /// [qualimap]  Activate this option to report information for the regions outside.
    #[structopt( short = "s",long = "outside-stats")]
    outside_stats_flag:bool,

    /// [qualimap]  Specific type of duplicated alignments to skip (if this option is activated).
    #[structopt( short = "u",long = "skip-dup-mode",default_value = "flagged",possible_values=&["flagged", "estimated", "both"])]
    skip_dup_mode:String,

    /// [qualimap]  Species to compare with genome GC distribution.
    #[structopt( short = "g",long = "genome-gc-distr",default_value = "",possible_values=&["","HUMAN(hg19)", "MOUSE(mm9)", "MOUSE(mm10)"])]
    gc_genome:String,

    /// [qualimap] Sequencing library protocol.
    #[structopt( short = "p", long = "sequencing-protocol",default_value = "non-strand-specific",possible_values=&["non-strand-specific", "strand-specific-forward",  "strand-specific-reverse"] )]
    protocol:String,

    /// [qualimap] Minimum size for a homopolymer to be considered in indel analysis.
    #[structopt(name = "min homopolymer size",short = "m",long = "hm",  default_value = "3")]
    min_homopolymer_size:usize,


}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct MetricsConfig {
    which: String,
    nthreads: usize,
    qualimap_config:QualimapConfig
}

impl MetricsConfig {
    pub fn new(args: &Arguments ) -> Self {
        let qualimap_config = MetricsConfig::construct_qualimap_config(args);

        Self {
            which:args.which.clone(),
            nthreads:args.nthreads,
            qualimap_config:qualimap_config
        }
    }

    pub fn construct_qualimap_config(args: &Arguments) ->QualimapConfig{
        let thread_num = args.nthreads;
        let window_num = args.nwindows;
        let bun_size = args.nchunks;
        let min_homopolymer_size = args.min_homopolymer_size;
        let feature_file = args.gff.clone();
        let outside_analyze_flag = args.outside_stats_flag;
        let lib_protocal = args.protocol.clone();
        let overlap_analyze_flag = args.collect_overlap_flag;
        let skip_dup_flag = args.skip_dup_flag;
        let skip_dup_mode = args.skip_dup_mode.clone();
        let gc_genome = args.gc_genome.clone();
        
        let mut config =  QualimapConfig::new();
        // set advanced options
        config.set_thread_num(thread_num);
        config.set_window_num(window_num);
        config.set_bunch_size(bun_size);
        config.set_min_homopolymer_size(min_homopolymer_size);

        // set region analyze
        if !feature_file.is_empty() {
            config.set_feature_file(feature_file);
            config.set_outside_region_analyze_flag(outside_analyze_flag)
        }

        // selt overlap detection flag
        config.set_collect_overlap_flag(overlap_analyze_flag);

        // set skip duplicates option
        config.set_skip_duplicate_flag(skip_dup_flag);
        config.set_skip_duplicate_mode(skip_dup_mode);

        // set reference gc content
        if !gc_genome.is_empty() {
            config.set_gc_genome(gc_genome);
        }
        config
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
