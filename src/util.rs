use regex::Regex;
use std::path::{Path};
use std::fs::{self, File, OpenOptions};
use std::io::{BufReader, BufWriter, Read, Write};
use std::str;

// These are values of experience.
const GZIP_COMPRESSION_RATIO: f64 = 5.5;
const ONE_READS_SIZE: u64 = 400;

pub fn is_remote_file(remote_path: &str) -> bool {
    let re = Regex::new(r"(oss|s3)://.*$").unwrap();
    re.is_match(remote_path)
}

pub fn parse_remote_path(remote_path: &str) -> (String, String, String) {
    // TODO: how to deal with exception when the filepath is not similar with oss://<bucket-name>/<filepath>
    let re = Regex::new(r"(oss|s3)://([^/]+)(/.*)$").unwrap();
    let cap = re.captures(remote_path).unwrap();
    return (
        String::from(&cap[1]),
        String::from(&cap[2]),
        String::from(&cap[3]),
    );
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

pub fn guess_nreads(fpath: &str) -> u64 {
    let path = Path::new(fpath);
    if path.exists() {
        let metadata = fs::metadata(fpath).unwrap();
        let fsize = metadata.len();
        let fext = path.extension().unwrap();
        if fext == "gz" {
            let total_fsize = fsize * GZIP_COMPRESSION_RATIO as u64;
            let nreads = total_fsize / ONE_READS_SIZE;
            return nreads;
        }
        return fsize;
    } else {
        return 0;
    }
}
