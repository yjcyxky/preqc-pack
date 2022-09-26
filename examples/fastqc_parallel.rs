extern crate preqc_pack;
use std::sync::Arc;

fn test_process_sequence() {
    let fastq_path = "./examples/test.fastq.gz";
    let pattern_file = "../data/patterns.bson";
    let (patterns, indexes, count) =
        preqc_pack::qc::mislabeling::VAFMatrix::read_patterns(pattern_file);
    let qc_results = preqc_pack::qc::QCResults::run_fastqc(
        fastq_path,
        Arc::new(patterns),
        Arc::new(indexes),
        count,
        1,
    );
    println!("{}", serde_json::to_string(&qc_results).unwrap());
}

fn main() {
    test_process_sequence();
}
