extern crate rust_htslib;

use preqc_pack::qc::bam::qualimap::Qualimap;
use preqc_pack::qc::config::bam_config::QualimapConfig;

fn main() {
    let bam_path = "examples/out_sorted_unique.bam";
    let config = QualimapConfig::new();
    let mut bam_qc = Qualimap::new();
    bam_qc.run(bam_path, &config);
}
