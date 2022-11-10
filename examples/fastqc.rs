extern crate preqc_pack;
use fastq::parse_path;
use preqc_pack::qc::fastqc::FastQC;

fn test_process_sequence() {
    let fastq_path = "examples/test.fastq.gz";

    let adapter_file = "data/adapter_list.txt";
    let contaminant_file = "data/contaminant_list.txt";
    let contaminants =
        preqc_pack::qc::fastqc::OverRepresentedSeqs::read_contaminants_file(contaminant_file);
    let adapters = preqc_pack::qc::fastqc::AdapterContent::read_adapter_file(adapter_file);

    let mut qc = FastQC::new(&contaminants, &adapters);

    parse_path(Some(fastq_path), |parser| {
        parser
            .each(|record| {
                qc.process_sequence(&record.to_owned_record());
                return true;
            })
            .unwrap()
    })
    .unwrap();

    println!("{:?}", qc);
}

fn main() {
    test_process_sequence();
}
