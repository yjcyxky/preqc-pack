extern crate preqc_pack;
use preqc_pack::qc::fastqc::BasicStats;
use std::fs::File;
use std::io::prelude::*;
use std::sync::Arc;

use fastq::OwnedRecord;

fn test_process_sequence() {
    let fastq_path = "examples/test.fastq.gz";
    let pattern_file = "data/patterns.json";
    let (patterns, indexes, count) =
        preqc_pack::qc::mislabeling::VAFMatrix::read_patterns(pattern_file);
    let mut count_vec: Vec<Option<usize>> = vec![None; count];
    for i in indexes {
        count_vec[i] = Some(0);
    }

    let adapter_file = "data/adapter_list.txt";
    let contaminant_file = "data/contaminant_list.txt";
    let contaminants =
        preqc_pack::qc::fastqc::OverRepresentedSeqs::read_contaminants_file(contaminant_file);
    let adapters = preqc_pack::qc::fastqc::AdapterContent::read_adapter_file(adapter_file);

    let fastqc_config =
        preqc_pack::qc::FastQCConfig::new(adapters, contaminants, None, None, None, None);

    let mislabeling_config = preqc_pack::qc::MislabelingConfig::new(patterns, count_vec, count);

    let qc_results = preqc_pack::qc::QCResults::run_qc_par(
        fastq_path,
        1,
        Arc::new("all".to_string()),
        Arc::new(fastqc_config),
        Arc::new(mislabeling_config),
    );
    let re = serde_json::to_string(&qc_results).unwrap();
    println!("{}", serde_json::to_string(&qc_results).unwrap());
    let mut f = File::create("./examples/result.json").unwrap();
    f.write_all(re.as_bytes()).unwrap();
}

fn test_process_sequence1() {
    let fastq_path = "examples/test.fastq.gz";
    let pattern_file = "data/patterns.json";
    let (patterns, indexes, count) =
        preqc_pack::qc::mislabeling::VAFMatrix::read_patterns(pattern_file);
    let mut count_vec: Vec<Option<usize>> = vec![None; count];
    for i in indexes {
        count_vec[i] = Some(0);
    }

    let adapter_file = "data/adapter_list.txt";
    let contaminant_file = "data/contaminant_list.txt";
    let contaminants =
        preqc_pack::qc::fastqc::OverRepresentedSeqs::read_contaminants_file(contaminant_file);
    let adapters = preqc_pack::qc::fastqc::AdapterContent::read_adapter_file(adapter_file);

    let fastqc_config =
        preqc_pack::qc::FastQCConfig::new(adapters, contaminants, None, None, None, None);

    let mislabeling_config = preqc_pack::qc::MislabelingConfig::new(patterns, count_vec, count);

    let qc_results = preqc_pack::qc::QCResults::run_qc_par(
        fastq_path,
        10,
        Arc::new("all".to_string()),
        Arc::new(fastqc_config),
        Arc::new(mislabeling_config),
    );
    let re = serde_json::to_string(&qc_results).unwrap();
    println!("{}", serde_json::to_string(&qc_results).unwrap());
    let mut f = File::create("./examples/result1.json").unwrap();
    f.write_all(re.as_bytes()).unwrap();
}

fn test_basic_statistics() {
    let read1 = OwnedRecord {
        head: b"some_name".to_vec(),
        seq: b"GTCGCACTGATCTGGGTTAGGCGCGGAGCCGAGGGTTGCACCATTTTTCATTATTGAATGCCAAGATA".to_vec(),
        qual: b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII".to_vec(),
        sep: None,
    };
    let mut a = BasicStats::new();

    let mut b = BasicStats::new();

    a.process_sequence(&read1);
    b.process_sequence(&read1);

    // a = b + a;
    // println!("{}", a.name);
}
fn main() {
    test_process_sequence();
    test_process_sequence1();
    test_basic_statistics();
}
