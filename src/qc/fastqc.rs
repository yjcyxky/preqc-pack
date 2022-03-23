use super::util::{BaseGroup, QualityCount};
use fastq::{OwnedRecord, Record};
use serde::{Deserialize, Serialize};

const SANGER_ENCODING_OFFSET: usize = 32;
const ILLUMINA_1_3_ENCODING_OFFSET: usize = 64;
const SANGER_ILLUMINA_1_9: &str = "Sanger / Illumina 1.9";
const ILLUMINA_1_3: &str = "Illumina 1.3";
const ILLUMINA_1_5: &str = "Illumina 1.5";

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct PhredEncoding {
  name: String,
  offset: usize,
}

impl PhredEncoding {
  pub fn new(name: &str, offset: usize) -> PhredEncoding {
    return PhredEncoding {
      name: name.to_string(),
      offset: offset,
    };
  }

  pub fn get_fastq_encoding_offset(acscii_num: usize) -> PhredEncoding {
    let lowest_char = char::from_u32(acscii_num as u32).unwrap();
    if acscii_num < 33 {
      panic!(
        "No known encodings with chars < 33 (Yours was {} with value {})",
        lowest_char, acscii_num
      );
    } else if acscii_num < 64 {
      return PhredEncoding::new(SANGER_ILLUMINA_1_9, SANGER_ENCODING_OFFSET);
    } else if acscii_num == ILLUMINA_1_3_ENCODING_OFFSET + 1 {
      return PhredEncoding::new(ILLUMINA_1_3, ILLUMINA_1_3_ENCODING_OFFSET);
    } else if acscii_num <= 126 {
      return PhredEncoding::new(ILLUMINA_1_5, ILLUMINA_1_3_ENCODING_OFFSET);
    }

    panic!(
      "No Known encodings with chars > 126 (Yours was {} with value {})",
      lowest_char, acscii_num
    );
  }

  pub fn convert_sanger_phred_to_probability(phred: usize) -> f32 {
    let base_10 = 10.0_f32;
    return base_10.powf(phred as f32 / -10.0);
  }

  pub fn convert_old_illumina_phred_to_probability(phred: usize) -> f32 {
    let base_10 = 10.0_f32;
    return base_10.powf((phred as f32 / phred as f32 + 1.0) / -10.0);
  }

  pub fn convert_probability_to_sanger_phred(p: f32) -> usize {
    return (-10.0_f32 * f32::log10(p)) as usize;
  }

  pub fn convert_probability_to_old_illumina_phred(p: f32) -> usize {
    return (-10.0_f32 * f32::log10(p / (1.0 - p))) as usize;
  }

  pub fn name(&self) -> String {
    return self.name.clone();
  }

  pub fn offset(&self) -> usize {
    return self.offset;
  }
}

#[cfg(test)]
mod phred_encoding_tests {
  use super::*;

  #[test]
  fn test_phred_encoding() {
    let phred = PhredEncoding::new("Illumina 1.3", 33);
    assert_eq!(phred.name, "Illumina 1.3".to_string());
    assert_eq!(phred.offset, 33);
  }

  #[test]
  fn test_convert_probability_to_old_illumina_phred() {
    let phred_score = PhredEncoding::convert_probability_to_old_illumina_phred(0.01);
    assert_eq!(phred_score, 19);
  }

  #[test]
  fn test_get_fastq_encoding_offset() {
    let phred = PhredEncoding::get_fastq_encoding_offset('A' as usize);
    assert_eq!(phred.offset, 64);
  }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct PerBaseSeqQuality {
  #[serde(skip_serializing)]
  quality_counts: Vec<QualityCount>,
  #[serde(skip_serializing)]
  base_pos: Vec<usize>,
  mean: Vec<f32>,
  median: Vec<f32>,
  lower_quartile: Vec<f32>,
  upper_quartile: Vec<f32>,
  lowest: Vec<f32>,
  highest: Vec<f32>,
  xlabels: Vec<String>,
}

impl PerBaseSeqQuality {
  pub fn new() -> PerBaseSeqQuality {
    return PerBaseSeqQuality {
      quality_counts: vec![],
      base_pos: vec![],
      mean: vec![],
      median: vec![],
      lower_quartile: vec![],
      upper_quartile: vec![],
      lowest: vec![],
      highest: vec![],
      xlabels: vec![],
    };
  }

  pub fn add_quality_counts(&mut self, quality_counts: &Vec<QualityCount>) {
    for i in 0..self.quality_counts.len() {
      self.quality_counts[i].add_quality_count(&quality_counts[i]);
    }
  }

  pub fn get_percentages(&mut self, offset: usize) {
    let groups: Vec<BaseGroup> = BaseGroup::make_base_groups(self.quality_counts.len());
    let length = groups.len();

    self.base_pos = (1..length + 1).collect();

    self.mean = vec![0.0; length];
    self.median = vec![0.0; length];

    self.lowest = vec![0.0; length];
    self.highest = vec![0.0; length];

    self.lower_quartile = vec![0.0; length];
    self.upper_quartile = vec![0.0; length];

    self.xlabels = vec!["".to_string(); length];

    for i in 0..length {
      let group = &groups[i];
      self.xlabels[i] = group.name();
      let min_base = group.lower_count();
      let max_base = group.upper_count();
      self.lowest[i] = self.get_percentile(min_base, max_base, offset, 10);
      self.highest[i] = self.get_percentile(min_base, max_base, offset, 90);
      self.mean[i] = self.get_mean(min_base, max_base, offset);
      self.median[i] = self.get_percentile(min_base, max_base, offset, 50);
      self.lower_quartile[i] = self.get_percentile(min_base, max_base, offset, 25);
      self.upper_quartile[i] = self.get_percentile(min_base, max_base, offset, 75);
    }
  }

  pub fn process_qual(&mut self, qual: &Vec<u8>) {
    let quality_counts_len = self.quality_counts.len();
    let qual_len = qual.len();
    if quality_counts_len < qual_len {
      for _ in quality_counts_len..qual_len {
        self.quality_counts.push(QualityCount::new());
      }
    }

    for i in 0..qual_len {
      self.quality_counts[i].add_value(qual[i] as usize);
    }
  }

  fn get_percentile(&self, minbp: usize, maxbp: usize, offset: usize, percentile: usize) -> f32 {
    let mut count: usize = 0;
    let mut total: usize = 0;

    for i in (minbp - 1)..maxbp {
      if self.quality_counts[i].total_counts() > 100 {
        count += 1;
        total += self.quality_counts[i].get_percentile(offset, percentile);
      }
    }

    if count > 0 {
      return total as f32 / count as f32;
    } else {
      // TODO: What value should select?
      return 0.0;
    }
  }

  fn get_mean(&self, minbp: usize, maxbp: usize, offset: usize) -> f32 {
    let mut count: usize = 0;
    let mut total: f32 = 0.0;

    for i in (minbp - 1)..maxbp {
      if self.quality_counts[i].total_counts() > 0 {
        count += 1;
        total += self.quality_counts[i].get_mean(offset);
      }
    }

    if count > 0 {
      return total / count as f32;
    }

    return 0.0;
  }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct BasicStats {
  name: String,
  total_reads: usize,
  total_bases: usize,
  t_count: usize,
  c_count: usize,
  g_count: usize,
  a_count: usize,
  n_count: usize,
  gc_percentage: f32,
  lowest_char: usize,
  highest_char: usize,
  file_type: String,
  min_length: usize,
  max_length: usize,
  phred: PhredEncoding,
}

impl BasicStats {
  fn new() -> BasicStats {
    return BasicStats {
      name: "".to_string(),
      total_reads: 0,
      total_bases: 0,
      t_count: 0,
      c_count: 0,
      g_count: 0,
      a_count: 0,
      n_count: 0,
      gc_percentage: 0.0,
      lowest_char: 126,
      highest_char: 0,
      file_type: "".to_string(),
      // We guess that the length of a sequence is impossible to be greater than 1000
      min_length: 1000,
      max_length: 0,
      phred: PhredEncoding::new("", 0),
    };
  }

  pub fn total_bases(&self) -> usize {
    return self.total_bases;
  }

  pub fn total_reads(&self) -> usize {
    return self.total_reads;
  }

  fn add_total_reads(&mut self, total_reads: usize) {
    self.total_reads += total_reads;
  }

  fn add_total_bases(&mut self, total_bases: usize) {
    self.total_bases += total_bases;
  }

  fn add_to_a_count(&mut self, a_count: usize) {
    self.a_count += a_count;
  }

  fn add_to_t_count(&mut self, t_count: usize) {
    self.t_count += t_count;
  }

  fn add_to_c_count(&mut self, c_count: usize) {
    self.c_count += c_count;
  }

  fn add_to_g_count(&mut self, g_count: usize) {
    self.g_count += g_count;
  }

  fn add_to_n_count(&mut self, n_count: usize) {
    self.n_count += n_count;
  }

  fn add_to_count(
    &mut self,
    a_count: usize,
    t_count: usize,
    c_count: usize,
    g_count: usize,
    n_count: usize,
  ) {
    self.a_count += a_count;
    self.t_count += t_count;
    self.c_count += c_count;
    self.g_count += g_count;
    self.n_count += n_count;
  }

  fn set_lowest_char(&mut self, c: usize) {
    self.lowest_char = c;
  }

  fn set_highest_char(&mut self, c: usize) {
    self.highest_char = c;
  }

  fn set_min_len(&mut self, seq_len: usize) {
    if seq_len < self.min_length {
      self.min_length = seq_len;
    }
  }

  fn set_max_len(&mut self, seq_len: usize) {
    if seq_len > self.max_length {
      self.max_length = seq_len;
    }
  }

  /// Guess the phred encoding based on the lowest char.
  ///
  /// NOTE: You must set the lowest char before running the set_phred method.
  ///
  fn set_phred(&mut self) {
    self.phred = PhredEncoding::get_fastq_encoding_offset(self.lowest_char);
  }

  /// Compute the gc percentage based on total_bases, g_count and c_count.
  ///
  /// NOTE: You must set the atcg base count before running the set_gc_percentage method.
  ///
  fn set_gc_percentage(&mut self) {
    self.gc_percentage = (self.g_count + self.c_count) as f32 / self.total_bases as f32;
  }

  fn finish(&mut self) {
    self.set_phred();
    self.set_gc_percentage();
  }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct FastQC {
  pub basic_stats: BasicStats,
  pub per_base_seq_quality: PerBaseSeqQuality,
}

impl FastQC {
  pub fn new() -> FastQC {
    return FastQC {
      basic_stats: BasicStats::new(),
      per_base_seq_quality: PerBaseSeqQuality::new(),
    };
  }

  /// Finish method is crucial, don't forget it.
  pub fn finish(&mut self) {
    self.basic_stats.finish();
    self
      .per_base_seq_quality
      .get_percentages(self.basic_stats.phred.offset);
  }

  pub fn set_highest_lowest_char(&mut self, qual: &Vec<u8>) {
    for c in qual {
      let num = c.clone() as usize;
      if self.basic_stats.lowest_char > num {
        self.basic_stats.set_lowest_char(num);
      }

      if self.basic_stats.highest_char < num {
        self.basic_stats.set_highest_char(num);
      }
    }
  }

  /// Process sequence one by one, and update the statistics data.
  ///
  /// A `record` contains head, seq, sep, qual fields.
  ///
  /// # Examples
  ///
  /// Basic usage:
  ///
  /// ```
  /// extern crate preqc_pack;
  /// use preqc_pack::qc::fastqc::FastQC;
  /// use fastq::OwnedRecord;
  ///
  /// let read1 = OwnedRecord {
  ///   head: b"some_name".to_vec(),
  ///   seq: b"GTCGCACTGATCTGGGTTAGGCGCGGAGCCGAGGGTTGCACCATTTTTCATTATTGAATGCCAAGATA".to_vec(),
  ///   qual: b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII".to_vec(),
  ///   sep: None,
  /// };
  ///
  /// let mut qc = FastQC::new();
  /// qc.process_sequence(&read1);
  ///
  /// assert_eq!(qc.basic_stats.total_bases(), 68);
  /// assert_eq!(qc.basic_stats.total_reads(), 1);
  /// // assert_eq!(qc.basic_stats.g_count, 20);
  /// // assert_eq!(qc.basic_stats.a_count, 15);
  /// // assert_eq!(qc.basic_stats.c_count, 14);
  /// // assert_eq!(qc.basic_stats.t_count, 19);
  /// // assert_eq!(qc.basic_stats.n_count, 0);
  /// ```
  ///
  pub fn process_sequence(&mut self, record: &OwnedRecord) {
    for base in record.seq() {
      // NOTE: to_string is expensive than from_utf8
      match std::str::from_utf8(&[base.clone()]).unwrap() {
        // match char::from(base.clone()).to_uppercase().to_string().as_str() {
        "T" => self.basic_stats.t_count += 1,
        "C" => self.basic_stats.c_count += 1,
        "G" => self.basic_stats.g_count += 1,
        "A" => self.basic_stats.a_count += 1,
        "N" => self.basic_stats.n_count += 1,
        _ => {}
      }
    }

    let seq_len = record.seq().len();
    self.basic_stats.total_bases += seq_len;
    self.basic_stats.set_min_len(seq_len);
    self.basic_stats.set_max_len(seq_len);

    self.set_highest_lowest_char(&record.qual);
    self.basic_stats.total_reads += 1;

    self.per_base_seq_quality.process_qual(&record.qual);
  }

  /// Merge several FastQC instances.
  ///
  /// You may get several FastQC instances When you handle fastq data with several threads.
  ///
  /// # Examples
  ///
  /// Basic usage:
  ///
  /// ```
  /// extern crate preqc_pack;
  /// use preqc_pack::qc::fastqc::FastQC;
  /// use fastq::OwnedRecord;
  ///
  /// let read1 = OwnedRecord {
  ///   head: b"some_name".to_vec(),
  ///   seq: b"GTCGCACTGATCTGGGTTAGGCGCGGAGCCGAGGGTTGCACCATTTTTCATTATTGAATGCCAAGATA".to_vec(),
  ///   qual: b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII".to_vec(),
  ///   sep: None,
  /// };
  ///
  /// let mut qc = FastQC::new();
  /// qc.process_sequence(&read1);
  ///
  /// let mut qc2 = FastQC::new();
  /// qc2.process_sequence(&read1);
  ///
  /// qc.merge(&[qc2]);
  /// assert_eq!(qc.basic_stats.total_bases(), 136);
  /// ```
  ///
  pub fn merge(&mut self, fastqc_vec: &[FastQC]) {
    for i in fastqc_vec {
      self.basic_stats.add_to_count(
        i.basic_stats.a_count,
        i.basic_stats.t_count,
        i.basic_stats.c_count,
        i.basic_stats.g_count,
        i.basic_stats.n_count,
      );

      self.basic_stats.add_total_bases(i.basic_stats.total_bases);
      self.basic_stats.add_total_reads(i.basic_stats.total_reads);
      self.basic_stats.set_min_len(i.basic_stats.min_length);
      self.basic_stats.set_max_len(i.basic_stats.max_length);
      self.basic_stats.set_lowest_char(i.basic_stats.lowest_char);
      self
        .per_base_seq_quality
        .add_quality_counts(&i.per_base_seq_quality.quality_counts);
    }

    // Finish method is crucial, don't forget it.
    self.finish();
  }
}

pub type FilteredFastQC = FastQC;
