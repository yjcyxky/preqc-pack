extern crate serde;

use digest::{Digest, Output};
use serde::{Deserialize, Serialize};
use std::io::Read;

const BUFFER_SIZE: usize = 51200;

#[derive(Serialize, Deserialize)]
pub struct Meta {
    md5sum: String,
    filesize: usize,
}

pub fn init_meta() -> Meta {
    return Meta {
        md5sum: String::from(""),
        filesize: 0,
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
