extern crate log;
extern crate structopt;
extern crate stderrlog;

mod cmd;

use structopt::StructOpt;
use cmd::meta;
use cmd::merge;

/// A suite of programs for interacting with fastq file
#[derive(StructOpt, Debug)]
#[structopt(setting=structopt::clap::AppSettings::ColoredHelp, name = "PreQC Tool Suite (preqc-pack)", author="Jingcheng Yang <yjcyxky@163.com>")]
struct Opt {
  /// A flag which control whether show more messages, true if used in the command line
  #[structopt(short="q", long="quiet")]
  quiet: bool,

  /// The number of occurrences of the `v/verbose` flag
  /// Verbose mode (-v, -vv, -vvv, etc.)
  #[structopt(short="v", long="verbose", parse(from_occurrences))]
  verbose: usize,

  #[structopt(subcommand)]
  cmd: SubCommands
}

#[derive(Debug, PartialEq, StructOpt)]
enum SubCommands {
  #[structopt(name="meta")]
  Meta(meta::Arguments),
  #[structopt(name="merge")]
  Merge(merge::Arguments)
}

fn main() {
  let opt = Opt::from_args();

  stderrlog::new()
    .module(module_path!())
    .modules(vec!["preqc_pack"])
    .quiet(opt.quiet)
    .verbosity(opt.verbose)
    .init()
    .unwrap();

  match opt.cmd {
    SubCommands::Meta(arguments) => {
      meta::run(&arguments);
    },
    SubCommands::Merge(arguments) => {
      merge::run(&arguments);
    }
  }
}
