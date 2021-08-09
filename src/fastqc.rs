use fastq::parse_path;
use serde::{Deserialize, Serialize};
use std::fs::{File, OpenOptions};
use std::io::{BufReader, BufWriter, Read, Write};
use std::str;

#[derive(Serialize, Deserialize)]
pub struct FastQC {
  datasize: usize,
}

pub fn init_fastqc(datasize: usize) -> FastQC {
  return FastQC { datasize: datasize };
}

#[allow(unused_variables)]
pub fn compute_data_size(filepath: &str) -> FastQC {
  let mut datasize: usize = 0;
  parse_path(Some(filepath), |parser| {
    let results: Vec<usize> = parser
      .parallel_each(1, |record_sets| {
        let mut thread_total = 0;
        for record_set in record_sets {
          thread_total += record_set.len();
        }

        thread_total
      })
      .expect("Not a valid fastq/fastq.gz file");
    datasize = results.iter().sum::<usize>();
  })
  .expect("Not a valid fastq/fastq.gz file");

  return FastQC { datasize: datasize };
}

pub fn zcat(infile: &str, output: &str) {
  let input = File::open(infile).unwrap();
  let mut reader = BufReader::new(input);
  let f = OpenOptions::new().append(true).create(true).open(output);
  let mut output = f.map(BufWriter::new).unwrap();

  let mut in_buf = [0; 1024 * 64];

  while let Ok(n) = reader.read(&mut in_buf) {
    if n == 0 {
      break;
    }

    output.write_all(&in_buf[..n]).unwrap();
  }
}
