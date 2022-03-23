extern crate preqc_pack;

fn test_process_sequence() {
  let fastq_path = "examples/test.fastq.gz";
  // let fastq_path = "/Users/choppy/Downloads/R20059384_merge_R1.fastq.gz";
  let qc_results = preqc_pack::qc::QCResults::run_fastqc(fastq_path, "data/patterns.bson");
  println!("{}", serde_json::to_string(&qc_results.vaf_matrix()).unwrap());
}

fn main() {
  test_process_sequence();
}
