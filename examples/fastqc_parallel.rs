extern crate preqc_pack;
use std::sync::Arc;

fn test_process_sequence() {
    let fastq_path = "./examples/test.fastq.gz";
    let pattern_file = "../data/patterns.bson";
    let (patterns, indexes, count) =
        preqc_pack::qc::mislabeling::VAFMatrix::read_patterns(pattern_file);
    let mut count_vec: Vec<Option<usize>> = vec![None; count];
    for i in indexes {
        count_vec[i] = Some(0);
    }
    let qc_results = preqc_pack::qc::QCResults::run_qc_par(
        fastq_path,
        Arc::new(patterns),
        Arc::new(count_vec),
        count,
        1,
        Arc::new("all".to_string())
    );
    println!("{}", serde_json::to_string(&qc_results).unwrap());
}

fn main() {
    test_process_sequence();
}
