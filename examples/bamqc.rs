use preqc_pack::qc::bam::qualimap::Qualimap;

fn main() {
    let bam_path = "examples/out_sorted_unique.bam";
    let mut bam_qc = Qualimap::new();
    bam_qc.run(bam_path.to_string());
}
