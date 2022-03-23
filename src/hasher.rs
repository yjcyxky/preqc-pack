use crate::util;
use digest::{Digest, Output};
use log::*;
use memmap2::{MmapAsRawDesc, MmapOptions};
use s3::bucket::Bucket;
use s3::creds::Credentials;
use s3::region::Region;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::Read;
use std::time::Instant;

const BUFFER_SIZE: usize = 51200;
#[derive(Serialize, Deserialize)]
pub struct Meta {
  md5sum: String,
  filesize: usize,
}

impl Meta {
  pub fn new(md5sum: &str, filesize: usize) -> Meta {
    Meta {
      md5sum: String::from(md5sum),
      filesize: filesize,
    }
  }
}

pub fn init_meta() -> Meta {
  return Meta {
    md5sum: String::from(""),
    filesize: 0,
  };
}

/// Compute digest value for given `Reader` and print it
/// On any error simply return without doing anything
pub fn process_mmap<D: Digest + Default>(reader: &File) -> Meta
where
  Output<D>: core::fmt::LowerHex,
{
  let mut sh = D::new();
  let mmap = unsafe { MmapOptions::new().map(reader).unwrap() };
  let filesize: usize = mmap.len();
  let n = filesize / BUFFER_SIZE;

  for i in 1..n + 1 {
    sh.update(&mmap[(i - 1) * BUFFER_SIZE..i * BUFFER_SIZE]);

    if filesize - i * BUFFER_SIZE < BUFFER_SIZE {
      sh.update(&mmap[i * BUFFER_SIZE..filesize])
    }
  }

  return Meta {
    filesize: filesize,
    md5sum: format!("{:x}", sh.finalize()),
  };
}

/// Compute digest value for given `Reader` and print it
/// On any error simply return without doing anything
pub fn process<D: Digest + Default, R: Read>(reader: &mut R) -> Meta
where
  Output<D>: core::fmt::LowerHex,
{
  let mut sh = D::new();
  let mut filesize: usize = 0;
  let mut buffer = [0u8; BUFFER_SIZE];

  loop {
    let n = match reader.read(&mut buffer) {
      Ok(n) => n,
      Err(_) => {
        return Meta {
          filesize: 0,
          md5sum: String::new(),
        }
      }
    };

    sh.update(&buffer[..n]);
    filesize += n;

    if n == 0 || n < BUFFER_SIZE {
      break;
    }
  }

  return Meta {
    filesize: filesize,
    md5sum: format!("{:x}", sh.finalize()),
  };
}

/// Computes the hash value for a remote file (such as locating on oss or s3 service).
///
/// # Example:
///
/// ```no_run
/// use preqc_pack::hasher;
/// use std::env;
/// use s3::creds::Credentials;
/// use s3::region::Region;
/// use md5::Md5;
///
/// let remote_path = "oss://choppy-app-example-data/RNAseq/test_fq/test_R2.fq.gz";
/// let region = Region::Custom {
///   region: "cn-shang".to_owned(),
///   endpoint: "https://oss-cn-shanghai.aliyuncs.com".to_owned(),
/// };
/// let credentials = Credentials::new(
///   Some(&env::var("OSS_ACCESS_KEY").unwrap()),
///   Some(&env::var("OSS_ACCESS_KEY_SECRET").unwrap()),
///   None,
///   None,
///   None,
///  ).unwrap();
///  hasher::process_remote::<Md5>(remote_path, region, credentials, 12, 8_388_608);
///
/// ```
pub async fn process_remote<D: Digest + Default>(
  remote_path: &str,
  region: Region,
  credentials: Credentials,
  n_threads: u64,
  chunk_size: u64,
) -> Meta
where
  Output<D>: core::fmt::LowerHex,
{
  let (_protocol, bucket_name, filepath) = util::parse_remote_path(remote_path);

  let bucket = Bucket::new(&bucket_name[..], region, credentials).unwrap();

  let (head, _code) = bucket.head_object(&filepath[..]).await.unwrap();
  let file_size: u64 = head.content_length.unwrap() as u64;
  let mut start: u64 = 0;
  let mut end: u64 = chunk_size;

  if chunk_size > file_size {
    end = file_size;
  }

  let mut sh: D = Digest::new();

  let nums = file_size / (chunk_size * n_threads) + 1;

  for n in 1..nums + 1 {
    let mut tasks = Vec::new();

    let join_all_time = Instant::now();
    for idx in 1..n_threads + 1 {
      print!(
        "start: {:?}, end: {:?}, n: {:?}, nums: {:?}, idx: {:?}, n_threads: {:?}, ",
        start, end, n, nums, idx, n_threads
      );

      let task = bucket.get_object_range(&filepath[..], start, Some(end));
      tasks.push(task);
      print!("The length of Vector: {:?}\n", tasks.len());

      if end == file_size {
        break;
      }

      start = idx * n * chunk_size + 1;
      if (idx + 1) * n * chunk_size > file_size {
        end = file_size;
      } else {
        end = (idx + 1) * n * chunk_size;
      }
    }

    for handle in tasks {
      let (temp, _code) = handle.await.unwrap();
      sh.update(&temp);
    }

    println!(
      "Joining all handles: {}ms",
      join_all_time.elapsed().as_millis()
    );
  }

  return Meta {
    filesize: file_size as usize,
    md5sum: format!("{:x}", sh.finalize()),
  };
}
