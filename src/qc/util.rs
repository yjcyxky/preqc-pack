use bson::Document;
use core::panic;
use std::fs::File;

use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct QualityCount {
  actual_counts: Vec<usize>,
  total_counts: usize,
}

impl QualityCount {
  pub fn new() -> QualityCount {
    return QualityCount {
      actual_counts: vec![0; 150],
      total_counts: 0,
    };
  }

  pub fn add_value(&mut self, c_ascii: usize) {
    self.total_counts += 1;
    self.actual_counts[c_ascii] += 1;
  }

  pub fn add_quality_count(&mut self, quality_count: &QualityCount) {
    self.total_counts += quality_count.total_counts();
    self.actual_counts = (0..self.actual_counts.len())
      .map(|i| self.actual_counts[i] + quality_count.actual_counts[i])
      .collect();
  }

  pub fn total_counts(&self) -> usize {
    return self.total_counts;
  }

  pub fn get_min_char(&self) -> char {
    for i in 0..self.actual_counts.len() {
      if self.actual_counts[i] > 0 {
        return char::from_u32(i as u32).unwrap();
      }
    }

    return char::from_u32(1000).unwrap();
  }

  pub fn get_max_char(&self) -> char {
    let length = self.actual_counts.len();
    for i in 0..length {
      let idx = length - 1 - i;
      if self.actual_counts[idx] > 0 {
        return char::from_u32(i as u32).unwrap();
      }
    }

    return char::from_u32(1000).unwrap();
  }

  pub fn get_mean(&self, offset: usize) -> f32 {
    let mut total: usize = 0;
    let mut count: usize = 0;
    let mut i = offset;

    while i < self.actual_counts.len() {
      total += self.actual_counts[i] * (i - offset);
      count += self.actual_counts[i];

      i += 1;
    }

    return (total / count) as f32;
  }

  pub fn get_percentile(&self, offset: usize, percentile: usize) -> usize {
    let mut total = self.total_counts;
    total *= percentile;
    total /= 100;

    let mut count: usize = 0;
    let mut i = offset;
    while i < self.actual_counts.len() {
      count += self.actual_counts[i];
      if count >= total {
        return i - offset;
      }

      i += 1;
    }

    return 0;
  }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct BaseGroup {
  name: String,
  lower_count: usize,
  upper_count: usize,
}

impl BaseGroup {
  pub fn new(lower_count: usize, upper_count: usize) -> BaseGroup {
    let name = if lower_count == upper_count {
      format!("{}", lower_count)
    } else {
      format!("{}-{}", lower_count, upper_count)
    };

    return BaseGroup {
      name: name,
      lower_count: lower_count,
      upper_count: upper_count,
    };
  }

  pub fn name(&self) -> String {
    return self.name.clone();
  }

  pub fn lower_count(&self) -> usize {
    return self.lower_count.clone();
  }

  pub fn upper_count(&self) -> usize {
    return self.upper_count.clone();
  }

  pub fn make_ungrouped_groups(max_length: usize) -> Vec<BaseGroup> {
    let mut starting_base: usize = 1;
    let interval: usize = 1;

    let mut groups: Vec<BaseGroup> = vec![];

    while starting_base <= max_length {
      let mut end_base = starting_base + (interval - 1);
      if end_base > max_length {
        end_base = max_length;
      }

      let bg: BaseGroup = BaseGroup::new(starting_base, end_base);
      groups.push(bg);

      starting_base += interval;
    }

    return groups;
  }

  pub fn make_base_groups(max_length: usize) -> Vec<BaseGroup> {
    return BaseGroup::make_linear_base_groups(max_length);
  }

  pub fn make_exponential_base_groups(max_length: usize) -> Vec<BaseGroup> {
    let mut starting_base: usize = 1;
    let mut interval: usize = 1;

    let mut groups: Vec<BaseGroup> = vec![];

    while starting_base <= max_length {
      let mut end_base = starting_base + (interval - 1);
      if end_base > max_length {
        end_base = max_length;
      }

      let bg = BaseGroup::new(starting_base, end_base);
      groups.push(bg);

      starting_base += interval;

      if starting_base == 10 && max_length > 75 {
        interval = 5;
      }

      if starting_base == 50 && max_length > 200 {
        interval = 10;
      }

      if starting_base == 100 && max_length > 300 {
        interval = 50;
      }

      if starting_base == 500 && max_length > 1000 {
        interval = 100;
      }

      if starting_base == 1000 && max_length > 2000 {
        interval = 500;
      }
    }

    return groups;
  }

  pub fn get_linear_interval(length: usize) -> usize {
    let base_values: Vec<usize> = vec![2, 5, 10];
    let mut multiplier: usize = 1;

    loop {
      for i in 0..base_values.len() {
        let interval = base_values[i] * multiplier;
        let mut group_count = 9 + (length - 9) / interval;
        if (length - 9) % interval != 0 {
          group_count += 1;
        }

        if group_count < 75 {
          return interval;
        }
      }

      multiplier *= 10;

      if multiplier == 10000000 {
        panic!(
          "Couldn't find a sensible interval grouping for length {}",
          length
        );
      }
    }
  }

  pub fn make_linear_base_groups(max_length: usize) -> Vec<BaseGroup> {
    if max_length <= 75 {
      return BaseGroup::make_ungrouped_groups(max_length);
    }

    let interval = BaseGroup::get_linear_interval(max_length);

    let mut starting_base = 1;
    let mut groups: Vec<BaseGroup> = vec![];

    while starting_base <= max_length {
      let mut end_base = starting_base + (interval - 1);

      if starting_base < 10 {
        end_base = starting_base;
      }

      if starting_base == 10 && interval > 10 {
        end_base = interval - 1;
      }

      if end_base > max_length {
        end_base = max_length;
      }

      let bg = BaseGroup::new(starting_base, end_base);
      groups.push(bg);

      if starting_base < 10 {
        starting_base += 1;
      } else if starting_base == 10 && interval > 10 {
        starting_base = interval;
      } else {
        starting_base += interval;
      }
    }

    return groups;
  }
}

/// Read SNP patterns from bson file.
/// 
/// NOTE: BSON Document -> HashMap<String, [usize; 2]>
/// 
/// # Arguments
/// 
/// * `pattern_file`: where is the pattern file.
/// 
/// # Example
/// 
/// ```
/// extern crate preqc_pack;
/// use preqc_pack::qc::util::read_patterns;
/// use bson::Bson;
/// 
/// let content = read_patterns("data/patterns.bson");
/// assert_eq!(&Bson::Array(vec![Bson::Int32(0), Bson::Int32(0)]), content.get("TCCTTGTCATATGTTTTTCTG").unwrap());
/// ```
/// 
pub fn read_patterns(pattern_file: &str) -> Document {
  let mut f = match File::open(pattern_file) {
    Ok(f) => f,
    Err(msg) => panic!("Cannot open {} - {}", pattern_file, msg),
  };

  return match Document::from_reader(&mut f) {
    Ok(content) => content,
    Err(msg) => {
      panic!("Cannot read pattern file {} - {}", pattern_file, msg)
    }
  };
}
