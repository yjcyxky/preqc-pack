extern crate preqc_pack;

fn test_process_sequence() {
  // let fastq_path = "/Users/choppy/Downloads/test.fq.gz";
  let fastq_path = "/Users/choppy/Downloads/R20059384_merge_R1.fastq.gz";
  preqc_pack::qc::run_fastqc(fastq_path);
}

fn main() {
  test_process_sequence();
}
