extern crate preqc_pack;
use fastq::parse_path;
use preqc_pack::qc::fastqc::FastQC;

fn test_process_sequence() {
  let fastq_path = "/Users/choppy/Downloads/test.fq.gz";
  let mut qc = FastQC::new();
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
