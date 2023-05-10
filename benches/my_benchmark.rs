use criterion::{black_box, criterion_group, criterion_main, Criterion};
use preqc_pack::qc::bamqc::BamStatsAnalysis;
fn fibonacci(n: u64) -> u64 {
    let mut a = 0;
    let mut b = 1;

    match n {
        0 => b,
        _ => {
            for _ in 0..n {
                let c = a + b;
                a = b;
                b = c;
            }
            b
        }
    }
}

fn test() {
    let bam_path = "examples/out_sorted_unique.bam";
    let bam_path01: &str = "examples/gambusino_reads.sorted.bam";
    let bam_path02: &str = "examples/pl_1_and_2.bam";
    let bam_path03: &str =
        "examples/R21045930-pool-1-dVF1-D5-1_combined_test_qdp_wes_rly_20221109_LCL5.sorted.deduped.bam";
    // Create bam qc instance
    let mut bam_stats_analysis = BamStatsAnalysis::new(bam_path03.to_string());
    // Run the analysis of bam qc
    bam_stats_analysis.run();
}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("PreQC-Pack Bam QC (1.1GB)", |b| b.iter(|| test()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
