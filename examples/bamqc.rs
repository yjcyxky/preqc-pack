extern crate rust_htslib;

use preqc_pack::qc::bam::qualimap::Qualimap;
use preqc_pack::qc::bam::fastbam::FastBam;
use preqc_pack::qc::config::bam_config::QualimapConfig;
use time::Instant;
use std::error::Error;

fn run_qualimap(bam_path:&str) -> Result<(), Box<dyn Error>> {
    let start_time = Instant::now();
    let config = QualimapConfig::new();
    let mut bam_qc = Qualimap::new();
    bam_qc.run(bam_path, &config);
    let _ =bam_qc.export("./examples");
    
    let end_time = Instant::now();
    println!("Time to run for qualimap: {:?}", end_time- start_time);

    Ok(())
}

fn run_fastbam(bam_path:&str) {
    let start_time = Instant::now();
    let mut bam_qc = FastBam::new();
    let config = QualimapConfig::new();
    bam_qc.run(bam_path, &config);

    let end_time = Instant::now();
    println!("Time to run for fastbam: {:?}", end_time- start_time);
}
fn main() {
    let bam_path = "examples/out_sorted_unique.bam";
    let _ =run_qualimap(bam_path);

    // run_fastbam(bam_path);
}
