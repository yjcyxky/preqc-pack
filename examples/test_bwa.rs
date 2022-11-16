extern crate bwa;

use bwa::BwaAligner;
use fastq::{parse_path, Record};
use std::io::Error;

fn main() {
    parse_path(Some("/data/preqc-pack/data/anc_R1.fastq.gz"), |parser| {
        let _result: Result<Vec<_>, Error> = parser.parallel_each(3, |record_sets| {
            let bwa =
                BwaAligner::from_path(&"/data/preqc-pack/references/fastq_screen/human/genome.fa")
                    .unwrap();
            for record_set in record_sets {
                for r in record_set.iter() {
                    let name = std::str::from_utf8(r.head()).unwrap();
                    let (r1_alns, _r2_alns) =
                        bwa.align_read_pair(r.head(), r.seq(), r.qual(), r.seq(), r.qual());
                    for i in r1_alns {
                        println!("{:?}\t{:?}", name, i.flags());
                    }
                }
            }
            true
        });
    })
    .expect("Not found the fastq file.");
}
