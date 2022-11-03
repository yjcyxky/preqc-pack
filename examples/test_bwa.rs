extern crate bwa;
use bwa::BwaAligner;

fn main() {
    let bwa = BwaAligner::from_path(&"/data/preqc-pack/references/fastq_screen/human/genome.fa").unwrap();

    let r1 = b"GATGGCTGCGCAAGGGTTCTTACTGATCGCCACGTTTTTACTGGTGTTAATGGTGCTGGCGCGTCCTTTAGGCAGCGGG";
    let q1 = b"2222222222222222222222222222222222222222222222222222222222222222222222222222222";
    let r2 = b"TGCTGCGTAGCAGATCGACCCAGGCATTCCCTAGCGTGCTCATGCTCTGGCTGGTAAACGCACGGATGAGGGCAAAAAT";
    let q2 = b"2222222222222222222222222222222222222222222222222222222222222222222222222222222";

    let (r1_alns, _r2_alns) = bwa.align_read_pair(b"read_name", r1, q1, r2, q2);
    println!("r1 mapping -- tid: {}, pos: {}", r1_alns[0].tid(), r1_alns[0].pos());
}