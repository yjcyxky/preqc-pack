extern crate preqc_pack;
use preqc_pack::qc::mislabeling::VAFMatrix;

fn test_read_patterns() {
  let (patterns, indexes, count) = VAFMatrix::read_patterns("data/patterns.bson");
  println!("{:?}, {:?}", count, patterns.get("TCCTTGTCATATGTTTTTCTG"));
}

fn main() {
  test_read_patterns();
}
