extern crate flate2;

use flate2::bufread::GzDecoder;
use seq_io::fastq::{Reader, Record};
use seq_io::parallel::parallel_fastq;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::BufReader;
use std::str;

const BUFFER_SIZE: usize = 512000;

#[derive(Serialize, Deserialize)]
pub struct FastQC {
  datasize: usize,
}

pub fn init_fastqc(datasize: usize) -> FastQC {
  return FastQC { datasize: datasize };
}

pub fn compute_data_size_par(filepath: &str) -> FastQC {
  let gz_reader = GzDecoder::new(BufReader::with_capacity(BUFFER_SIZE, File::open(filepath).unwrap()))
    .expect("Not valid gzip format.");
  let reader = Reader::with_capacity(gz_reader, BUFFER_SIZE);

  let mut datasize: usize = 0;
  parallel_fastq(
    reader,
    4,
    2,
    |record, found| {
      // runs in worker
      *found = record.seq().len() > 0;
    },
    |_, found| {
      // runs in main thread
      if *found {
        datasize += 1;
      }

      // Some(value) will stop the reader, and the value will be returned.
      // In the case of never stopping, we need to give the compiler a hint about the
      // type parameter, thus the special 'turbofish' notation is needed,
      // hoping on progress here: https://github.com/rust-lang/rust/issues/27336
      None::<()>
    },
  )
  .unwrap();

  return FastQC {
    datasize: datasize * 150,
  };
}

#[allow(unused_variables)]
pub fn compute_gz_data_size(filepath: &str) -> FastQC {
  let gz_reader = GzDecoder::new(File::open(filepath).map(BufReader::new).unwrap())
    .expect("Not valid gzip format.");
  let mut reader = Reader::new(gz_reader);

  let mut datasize: usize = 0;
  while let Some(record) = reader.next() {
    // let record = record.expect("Error reading record");
    // println!("{}", str::from_utf8(record.seq()).unwrap());
    datasize += 1;
  }

  return FastQC {
    datasize: datasize * 150,
  };
}

#[allow(unused_variables)]
pub fn compute_data_size(filepath: &str) -> FastQC {
  let mut reader = Reader::from_path(filepath).unwrap();

  let mut datasize: usize = 0;
  while let Some(record) = reader.next() {
    let record = record.expect("Error reading record");
    // TODO: More actions?
    // str::from_utf8(record.seq()).unwrap();
    datasize += 1;
  }

  return FastQC {
    datasize: datasize * 150,
  };
}
