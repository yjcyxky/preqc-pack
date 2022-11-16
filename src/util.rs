use regex::Regex;
use std::fs::{self, File, OpenOptions};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::str;

pub fn is_remote_file(remote_path: &str) -> bool {
    let re = Regex::new(r"(oss|s3)://.*$").unwrap();
    re.is_match(remote_path)
}

pub fn parse_remote_path(remote_path: &str) -> (String, String, String) {
    // TODO: how to deal with exception when the filepath is not similar with oss://<bucket-name>/<filepath>
    let re = Regex::new(r"(oss|s3)://([^/]+)(/.*)$").unwrap();
    let cap = re.captures(remote_path).unwrap();
    (
        String::from(&cap[1]),
        String::from(&cap[2]),
        String::from(&cap[3]),
    )
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
        // Bytes
        let fsize = metadata.len();
        let fext = path.extension().unwrap();
        if fext == "gz" {
            let nreads = 1968523.0 + 0.01164386 * fsize as f64;
            return nreads as u64;
        }

        fsize
    } else {
        0
    }
}
