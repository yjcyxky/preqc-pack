extern crate preqc_pack;
use preqc_pack::qc::util::read_patterns;

fn test_read_patterns() {
  let content = read_patterns("data/patterns.bson");
  println!("{:?}", content.get("TCCTTGTCATATGTTTTTCTG"));
}

fn main() {
  test_read_patterns();
}
