extern crate log;
extern crate stderrlog;
extern crate structopt;

mod cmd;

use cmd::merge;
use cmd::metrics;
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
    #[structopt(short = "t", long = "timestamp")]
    ts: Option<stderrlog::Timestamp>,

    #[structopt(subcommand)]
    cmd: SubCommands,
}

#[derive(Debug, PartialEq, StructOpt)]
enum SubCommands {
    #[structopt(name = "metrics")]
    Meta(metrics::Arguments),
    #[structopt(name = "merge")]
    Merge(merge::Arguments),
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
        SubCommands::Meta(arguments) => {
            metrics::run(&arguments);
        }
        SubCommands::Merge(arguments) => {
            merge::run(&arguments);
        }
    }
}
