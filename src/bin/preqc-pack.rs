extern crate log;
extern crate stderrlog;
extern crate structopt;

mod cmd;

use cmd::merge;
use cmd::metrics::{bam_metrics, fastq_metrics};
use structopt::StructOpt;

/// A suite of qc programs for interacting with fastq/bam/vcf/exp file
#[derive(StructOpt, Debug)]
#[structopt(setting=structopt::clap::AppSettings::ColoredHelp, name = "PreQC Tool Suite (preqc-pack)", author="Jingcheng Yang <yjcyxky@163.com>; Haonan Chen <haonanchen0815@163.com>")]
struct Opt {
    /// A flag which control whether show more messages, true if used in the command line
    #[structopt(short = "q", long = "quiet")]
    quiet: bool,

    /// The number of occurrences of the `v/verbose` flag
    /// Verbose mode (-v/Debug, -vv/Trace, etc.)
    #[structopt(short = "v", long = "verbose", parse(from_occurrences))]
    verbose: usize,

    /// Timestamp (sec, ms, ns, none)
    #[structopt(short = "t", long = "timestamp")] // Option makes it optional
    ts: Option<stderrlog::Timestamp>,

    #[structopt(subcommand)] // Note that we mark a field as a subcommand
    cmd: SubCommands,
}

#[derive(Debug, PartialEq, StructOpt)]
enum SubCommands {
    #[structopt(name = "metrics")]
    Meta(Metrics),
    #[structopt(name = "merge")]
    Merge(merge::Arguments),
}

/// Run qc at different stages, such as fastq, bam, vcf
#[derive(Debug, PartialEq, StructOpt)]
#[structopt(setting=structopt::clap::AppSettings::ColoredHelp)]
struct Metrics {
    #[structopt(subcommand)]
    metrics_type: MetricsType,
}

#[derive(Debug, PartialEq, StructOpt)]
enum MetricsType {
    #[structopt(name = "fastq")]
    Fastq(fastq_metrics::Arguments),
    #[structopt(name = "bam")]
    Bam(bam_metrics::Arguments),
}

fn main() {
    let opt = Opt::from_args();

    stderrlog::new()
        .module(module_path!())
        .quiet(opt.quiet)
        .show_module_names(true)
        .verbosity(opt.verbose + 2)
        .timestamp(opt.ts.unwrap_or(stderrlog::Timestamp::Second))
        .init()
        .unwrap();

    match opt.cmd {
        SubCommands::Meta(metrics) => {
            // metrics::run(&arguments);
            match metrics.metrics_type {
                MetricsType::Fastq(arguments) => {
                    fastq_metrics::run(&arguments);
                }
                MetricsType::Bam(arguments) => {
                    bam_metrics::run(&arguments);
                }
            }
        }
        SubCommands::Merge(arguments) => {
            merge::run(&arguments);
        }
    }
}
