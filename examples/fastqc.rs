extern crate preqc_pack;
use fastq::parse_path;
use preqc_pack::qc::fastq::fastqc::{AdapterContent, FastQC, OverRepresentedSeqs};

fn run_fastqc() {
    let fastq_path = "./examples/test.fastq.gz";
    let adapter_file = "data/adapter_list.txt";
    let contaminant_file = "data/contaminant_list.txt";
    let contaminants = OverRepresentedSeqs::read_contaminants_file(contaminant_file);
    let adapters = AdapterContent::read_adapter_file(adapter_file);

    let mut qc = FastQC::new(&contaminants, &adapters, None, None, None, None);

    parse_path(Some(fastq_path), |parser| {
        parser
            .each(|record| {
                qc.process_sequence(&record.to_owned_record());
                return true;
            })
            .unwrap()
    })
    .unwrap();
    qc.finish();

    // export basics results
    qc.export_to_txt("./examples/");
}

fn main() {
    
}
