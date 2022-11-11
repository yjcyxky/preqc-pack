use fastq::Record;
use probability::prelude::*;
use regex::Regex;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::io::Read;
use std::{
    cmp,
    collections::HashMap,
    f64::consts::{E, PI},
    str::from_utf8,
    vec,
};

const SANGER_ENCODING_OFFSET: usize = 33;
const ILLUMINA_1_3_ENCODING_OFFSET: usize = 64;
const SANGER_ILLUMINA_1_9: &str = "Sanger / Illumina 1.9";
const ILLUMINA_1_3: &str = "Illumina 1.3";
const ILLUMINA_1_5: &str = "Illumina 1.5";

const FORWARD_TYPE: usize = 1;
const REVERSE_TYPE: usize = 2;

const FASTQC_CONFIG_DUP_LENGTH: usize = 0;
const FASTQC_CONFIG_KMER_SIZE: usize = 0;
const FASTQC_CONFIG_CONTAMINANT_FILE: &str = "";
const FASTQC_CONFIG_ADAPTER_FILE: &str = "";

const INDICATOR_CONFIG_TILE_IGNORE: usize = 0;
const INDICATOR_CONFIG_OVERREPESENTED_WARN: f64 = 0.1;

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

    pub fn merge(&mut self, other: &QualityCount) {
        self.total_counts += other.total_counts;

        for i in 0..150 {
            self.actual_counts[i] += other.actual_counts[i];
        }
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
                return char::from_u32(idx as u32).unwrap();
            }
        }

        return char::from_u32(1000).unwrap();
    }

    pub fn get_mean(&self, offset: usize) -> f64 {
        let mut total: usize = 0;
        let mut count: usize = 0;
        let mut i = offset;

        while i < self.actual_counts.len() {
            total += self.actual_counts[i] * (i - offset);
            count += self.actual_counts[i];

            i += 1;
        }

        return total as f64 / count as f64;
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

    pub fn convert_sanger_phred_to_probability(phred: usize) -> f64 {
        let base_10 = 10.0_f64;
        return base_10.powf(phred as f64 / -10.0);
    }

    pub fn convert_old_illumina_phred_to_probability(phred: usize) -> f64 {
        let base_10 = 10.0_f64;
        return base_10.powf((phred as f64 / phred as f64 + 1.0) / -10.0);
    }

    pub fn convert_probability_to_sanger_phred(p: f64) -> usize {
        return (-10.0_f64 * f64::log10(p)) as usize;
    }

    pub fn convert_probability_to_old_illumina_phred(p: f64) -> usize {
        return (-10.0_f64 * f64::log10(p / (1.0 - p))) as usize;
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
    mean: Vec<f64>,
    median: Vec<f64>,
    lower_quartile: Vec<f64>,
    upper_quartile: Vec<f64>,
    lowest: Vec<f64>,
    highest: Vec<f64>,
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

    pub fn process_sequence(&mut self, record: &impl Record) {
        let qual = record.qual();
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

    fn get_percentile(&self, minbp: usize, maxbp: usize, offset: usize, percentile: usize) -> f64 {
        let mut count: usize = 0;
        let mut total: usize = 0;

        for i in (minbp - 1)..maxbp {
            if self.quality_counts[i].total_counts() > 100 {
                count += 1;
                total += self.quality_counts[i].get_percentile(offset, percentile);
            }
        }

        if count > 0 {
            return total as f64 / count as f64;
        } else {
            // TODO: What value should select?
            return 0.0;
        }
    }

    fn get_mean(&self, minbp: usize, maxbp: usize, offset: usize) -> f64 {
        let mut count: usize = 0;
        let mut total: f64 = 0.0;

        for i in (minbp - 1)..maxbp {
            if self.quality_counts[i].total_counts() > 0 {
                count += 1;
                total += self.quality_counts[i].get_mean(offset);
            }
        }

        if count > 0 {
            return total / count as f64;
        }

        return 0.0;
    }

    pub fn add_quality_counts(&mut self, quality_counts: &Vec<QualityCount>) {
        for i in 0..self.quality_counts.len() {
            self.quality_counts[i].add_quality_count(&quality_counts[i]);
        }
    }

    /// Some data values should be calculated after all sequences have been processed
    fn finish(&mut self, offset: usize) {
        // calculate statistics data from all quality scores
        self.get_percentages(offset);
    }

    pub fn merge(&mut self, other: &PerBaseSeqQuality) {
        self.add_quality_counts(&other.quality_counts);
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
    gc_percentage: f64,
    lowest_char: usize,
    highest_char: usize,
    file_type: String,
    min_length: usize,
    max_length: usize,
    phred: PhredEncoding,
}

impl BasicStats {
    pub fn new() -> BasicStats {
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

    pub fn process_sequence(&mut self, record: &impl Record) {
        let mut seq_len = 0;
        for base in record.seq() {
            let base_char = *base as char;
            match base_char {
                // match char::from(base.clone()).to_uppercase().to_string().as_str() {
                'T' => {
                    self.t_count += 1;
                    seq_len += 1;
                }
                'C' => {
                    self.c_count += 1;
                    seq_len += 1;
                }
                'G' => {
                    self.g_count += 1;
                    seq_len += 1;
                }
                'A' => {
                    self.a_count += 1;
                    seq_len += 1;
                }
                'N' => {
                    self.n_count += 1;
                    seq_len += 1;
                }
                _ => {}
            }
        }

        self.update_highest_lowest_char(&record.qual());
        self.total_bases += seq_len;
        self.total_reads += 1;

        if self.total_reads == 1 {
            self.min_length = seq_len;
            self.max_length = seq_len;
        } else {
            if seq_len < self.min_length {
                self.min_length = seq_len;
            }
            if seq_len > self.max_length {
                self.max_length = seq_len;
            }
        }
    }

    pub fn update_highest_lowest_char(&mut self, qual: &[u8]) {
        for c in qual {
            let num = c.clone() as usize;
            if self.lowest_char > num {
                self.lowest_char = num;
            }

            if self.highest_char < num {
                self.highest_char = num;
            }
        }
    }

    pub fn update_name(mut self, filename: &str) -> BasicStats {
        self.name = filename.to_string();
        self
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
        self.gc_percentage = (self.g_count + self.c_count) as f64 / self.total_bases as f64;
    }

    /// Some data values should be calculated after all sequences have been processed
    fn finish(&mut self) {
        self.set_phred();
        self.set_gc_percentage();
    }

    pub fn merge(&mut self, other: &BasicStats) {
        self.total_reads += other.total_reads;
        self.total_bases += other.total_bases;
        self.t_count += other.t_count;
        self.c_count += other.c_count;
        self.a_count += other.a_count;
        self.g_count += other.g_count;
        self.n_count += other.n_count;
        self.min_length = self.min_length.min(other.min_length);
        self.max_length = self.max_length.max(other.max_length);
        self.lowest_char = self.lowest_char.min(other.lowest_char);
        self.highest_char = self.highest_char.max(other.highest_char);
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct PerSeqQualityScore {
    average_score_counts: HashMap<usize, usize>,
    y_category_count: Vec<usize>,
    x_category_quality: Vec<usize>,
    max_counts: usize,
    most_frequent_score: usize,
    lowest_char: usize,
}

impl PerSeqQualityScore {
    fn new() -> PerSeqQualityScore {
        return PerSeqQualityScore {
            average_score_counts: HashMap::new(),
            y_category_count: vec![],
            x_category_quality: vec![],
            max_counts: 0,
            most_frequent_score: 0,
            lowest_char: 126,
        };
    }

    // analysis the average quality scores for a sequence and update the average_score_counts
    fn process_sequence(&mut self, record: &impl Record) {
        let mut average_quality = 0;
        for c in record.qual() {
            let num = *c as usize;
            if num < self.lowest_char {
                self.lowest_char = num;
            }
            average_quality += num;
        }

        if record.qual().len() > 0 {
            average_quality = average_quality / record.qual().len();

            if self.average_score_counts.contains_key(&average_quality) {
                let mut current_count = self.average_score_counts[&average_quality];
                current_count += 1;
                self.average_score_counts
                    .insert(average_quality, current_count);
            } else {
                self.average_score_counts.insert(average_quality, 1);
            }
        }
    }

    fn calculate_distribution(&mut self) {
        let encoding = PhredEncoding::get_fastq_encoding_offset(self.lowest_char);

        let mut raw_scores = self
            .average_score_counts
            .keys()
            .copied()
            .collect::<Vec<_>>();
        raw_scores.sort();

        self.y_category_count = vec![0; raw_scores[raw_scores.len() - 1] - raw_scores[0] + 1];
        self.x_category_quality = vec![0; self.y_category_count.len()];

        for i in 0..self.y_category_count.len() {
            self.x_category_quality[i] = (raw_scores[0] + i) - encoding.offset();
            if self.average_score_counts.contains_key(&(raw_scores[0] + i)) {
                self.y_category_count[i] = self.average_score_counts[&(raw_scores[0] + i)];
            }
        }

        for i in 0..raw_scores.len() {
            if self.y_category_count[i] > self.max_counts {
                self.max_counts = self.y_category_count[i];
                self.most_frequent_score = self.x_category_quality[i];
            }
        }
    }

    /// Some data values should be calculated after all sequences have been processed
    fn finish(&mut self) {
        self.calculate_distribution();
    }

    pub fn merge(&mut self, other: &PerSeqQualityScore) {
        self.lowest_char = self.lowest_char.min(other.lowest_char);

        // merge average_score_counts
        for (score, count) in other.average_score_counts.clone() {
            if self.average_score_counts.contains_key(&score) {
                let mut current_count = self.average_score_counts[&score];
                current_count += count;
                self.average_score_counts.insert(score, current_count);
            } else {
                self.average_score_counts.insert(score, count);
            }
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct PerBaseSeqContent {
    g_counts: Vec<usize>,
    c_counts: Vec<usize>,
    a_counts: Vec<usize>,
    t_counts: Vec<usize>,
    percentages: Vec<Vec<f64>>,
    x_category: Vec<String>,
}

impl PerBaseSeqContent {
    fn new() -> PerBaseSeqContent {
        return PerBaseSeqContent {
            g_counts: vec![],
            c_counts: vec![],
            a_counts: vec![],
            t_counts: vec![],
            percentages: vec![],
            x_category: vec![],
        };
    }

    fn get_percentages(&mut self) {
        let groups: Vec<BaseGroup> = BaseGroup::make_base_groups(self.g_counts.len());
        let length = groups.len();

        let mut g_percent = vec![0.0; length];
        let mut a_percent = vec![0.0; length];
        let mut t_percent = vec![0.0; length];
        let mut c_percent = vec![0.0; length];

        let mut g_count = 0;
        let mut a_count = 0;
        let mut t_count = 0;
        let mut c_count = 0;
        let mut total = 0;
        for i in 0..length {
            self.x_category.push(groups[i].name());
            g_count = 0;
            a_count = 0;
            t_count = 0;
            c_count = 0;
            total = 0;
            let current_group = &groups[i];
            for bp in current_group.lower_count() - 1..current_group.upper_count() {
                total += self.g_counts[bp];
                total += self.c_counts[bp];
                total += self.a_counts[bp];
                total += self.t_counts[bp];

                g_count += self.g_counts[bp];
                c_count += self.c_counts[bp];
                a_count += self.a_counts[bp];
                t_count += self.t_counts[bp];
            }

            g_percent[i] = (g_count as f64 / total as f64) * 100 as f64;
            a_percent[i] = (a_count as f64 / total as f64) * 100 as f64;
            t_percent[i] = (t_count as f64 / total as f64) * 100 as f64;
            c_percent[i] = (c_count as f64 / total as f64) * 100 as f64;
        }

        self.percentages = vec![t_percent, c_percent, a_percent, g_percent];
    }

    fn process_sequence(&mut self, record: &impl Record) {
        let seq = record.seq();
        let seq_len = seq.len();
        let g_counts_len = self.g_counts.len();

        if g_counts_len < seq_len {
            for _ in g_counts_len..seq_len {
                self.g_counts.push(0);
                self.a_counts.push(0);
                self.c_counts.push(0);
                self.t_counts.push(0);
            }
        }

        for i in 0..seq_len {
            let base_char = seq[i] as char;
            match base_char {
                // match char::from(base.clone()).to_uppercase().to_string().as_str() {
                'T' => {
                    self.t_counts[i] += 1;
                }
                'C' => {
                    self.c_counts[i] += 1;
                }
                'G' => {
                    self.g_counts[i] += 1;
                }
                'A' => {
                    self.a_counts[i] += 1;
                }
                _ => {}
            }
        }
    }

    /// Some data values should be calculated after all sequences have been processed
    fn finish(&mut self) {
        self.get_percentages();
    }

    pub fn merge(&mut self, other: &PerBaseSeqContent) {
        let this_len = self.g_counts.len();
        let other_len = other.g_counts.len();
        if this_len < other_len {
            for _ in this_len..other_len {
                self.g_counts.push(0);
                self.a_counts.push(0);
                self.c_counts.push(0);
                self.t_counts.push(0);
            }
        }

        for i in 0..other_len {
            self.t_counts[i] += other.t_counts[i];
            self.a_counts[i] += other.a_counts[i];
            self.c_counts[i] += other.c_counts[i];
            self.g_counts[i] += other.g_counts[i];
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct GCModelValue {
    percentage: usize,
    increment: f64,
}

impl GCModelValue {
    fn new(_percentage: usize, _increment: f64) -> GCModelValue {
        return GCModelValue {
            percentage: _percentage,
            increment: _increment,
        };
    }

    pub fn percentage(&self) -> usize {
        return self.percentage;
    }

    pub fn increment(&self) -> f64 {
        return self.increment;
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct GCModel {
    read_length: usize,
    models: Vec<Vec<GCModelValue>>,
}

impl GCModel {
    fn new() -> GCModel {
        return GCModel {
            read_length: 0,
            models: vec![vec![]],
        };
    }

    fn new_by_len(length: usize) -> GCModel {
        let mut claim_counts = vec![0; 101];
        let read_length = length;
        let mut models = vec![vec![]; length + 1];

        for pos in 0..length + 1 {
            let mut low_count = (pos as f64 - 0.5) as f64;
            let mut high_count = (pos as f64 + 0.5) as f64;

            if low_count < 0.0 {
                low_count = 0.0;
            }
            if high_count < 0.0 {
                high_count = 0.0;
            }
            if high_count > length as f64 {
                high_count = length as f64;
            }
            if low_count > length as f64 {
                low_count = length as f64;
            }
            let low_percent = (low_count * 100 as f64 / length as f64).round() as usize;
            let high_percent = (high_count * 100 as f64 / length as f64).round() as usize;
            for p in low_percent..high_percent + 1 {
                claim_counts[p] += 1;
            }
        }

        // We now do a second pass to make up the model using the weightings
        // we calculated previously.

        for pos in 0..length + 1 {
            let mut low_count = (pos as f64 - 0.5) as f64;
            let mut high_count = (pos as f64 + 0.5) as f64;
            if low_count < 0.0 {
                low_count = 0.0;
            }
            if high_count < 0.0 {
                high_count = 0.0;
            }
            if high_count > length as f64 {
                high_count = length as f64;
            }
            if low_count > length as f64 {
                low_count = length as f64;
            }

            let low_percent = (low_count * 100 as f64 / length as f64).round() as usize;
            let high_percent = (high_count * 100 as f64 / length as f64).round() as usize;

            let mut model_values: Vec<GCModelValue> =
                Vec::with_capacity(high_percent - low_percent + 1);
            for p in low_percent..high_percent + 1 {
                model_values.insert(
                    p - low_percent,
                    GCModelValue::new(p, 1 as f64 / claim_counts[p] as f64),
                );
            }
            models[pos] = model_values;
        }
        return GCModel {
            read_length: read_length,
            models: models,
        };
    }

    pub fn get_model_values(&self, gc_count: usize) -> &Vec<GCModelValue> {
        return self.models.get(gc_count).unwrap();
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct NormalDistribution {
    mean: f64,
    stdev: f64,
}

impl NormalDistribution {
    pub fn new(_mean: f64, _stdev: f64) -> NormalDistribution {
        return NormalDistribution {
            mean: _mean,
            stdev: _stdev,
        };
    }

    pub fn get_zscore_for_values(&self, value: f64) -> f64 {
        let _stdev = self.stdev;
        let lhs = 1.0 / (2.0 * _stdev * _stdev * PI).sqrt();
        let rhs = E.powf(0.0 - (value - self.mean).powf(2.0) / (2.0 * _stdev * _stdev));

        return lhs * rhs;
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct PerSeqGCContent {
    x_category: Vec<usize>,
    y_gc_distribution: Vec<f64>,
    y_theo_distribution: Vec<f64>,
    max: f64,
    deviation_percent: f64,
    #[serde(skip_serializing)]
    cached_models: Vec<GCModel>,
}

impl PerSeqGCContent {
    fn new() -> PerSeqGCContent {
        return PerSeqGCContent {
            max: 0.0,
            deviation_percent: 0.0,
            x_category: vec![],
            y_gc_distribution: vec![0.0; 101],
            y_theo_distribution: vec![0.0; 101],
            cached_models: Vec::with_capacity(200),
        };
    }

    fn process_sequence(&mut self, record: &impl Record) {
        // Because we keep a model around for every possible sequence length we
        // encounter we need to reduce the number of models.  We can do this by
        // rounding off the sequence once we get above a certain size

        let seq = self.truncate_sequence(record);
        let this_seq_length = seq.len();
        if this_seq_length == 0 {
            return;
        }

        let mut this_seq_gc_count = 0;
        for i in 0..this_seq_length {
            let base_char = seq[i] as char;
            if base_char == 'G' || base_char == 'C' {
                this_seq_gc_count += 1;
            }
        }

        let cached_models_len = self.cached_models.len();
        if cached_models_len <= this_seq_length {
            for _ in cached_models_len..this_seq_length {
                self.cached_models.push(GCModel::new());
            }
        }

        match self.cached_models.get(this_seq_length) {
            None => {
                self.cached_models
                    .push(GCModel::new_by_len(this_seq_length));
            }
            _ => {}
        }

        let values: &Vec<GCModelValue> =
            self.cached_models[this_seq_length].get_model_values(this_seq_gc_count);

        for i in 0..values.len() {
            self.y_gc_distribution[values[i].percentage()] += values[i].increment();
        }
    }

    fn truncate_sequence<'a>(&'a mut self, record: &'a impl Record) -> &[u8] {
        let _seq = record.seq();
        let seq_len = _seq.len();
        if seq_len > 1000 {
            let length = (seq_len / 1000) * 1000;
            return &_seq[0..length];
        }

        if seq_len > 100 {
            let length = (seq_len / 100) * 100;
            return &_seq[0..length];
        }
        return _seq;
    }

    fn calculate_distribution(&mut self) {
        self.max = 0.0;
        self.x_category = vec![0; self.y_gc_distribution.len()];
        let mut total_count: f64 = 0.0;
        // We use the mode to calculate the theoretical distribution
        // so that we cope better with skewed distributions.
        let mut first_mode = 0;
        let mut mode_count: f64 = 0.0;

        for i in 0..self.y_gc_distribution.len() {
            self.x_category[i] = i;
            total_count += self.y_gc_distribution[i];

            if self.y_gc_distribution[i] > mode_count {
                mode_count = self.y_gc_distribution[i];
                first_mode = i;
            }
            if self.y_gc_distribution[i] > self.max {
                self.max = self.y_gc_distribution[i];
            }
        }

        // The mode might not be a very good measure of the centre
        // of the distribution either due to duplicated vales or
        // several very similar values next to each other.  We therefore
        // average over adjacent points which stay above 95% of the modal
        // value

        let mut mode: f64 = 0.0;
        let mut mode_duplicate = 0;
        let mut fell_off_top = true;

        for i in first_mode..self.y_gc_distribution.len() {
            if self.y_gc_distribution[i]
                > self.y_gc_distribution[first_mode]
                    - (self.y_gc_distribution[first_mode] / 10 as f64)
            {
                mode += i as f64;
                mode_duplicate += 1;
            } else {
                fell_off_top = false;
                break;
            }
        }

        let mut fell_off_bottom = true;

        for i in (0..first_mode).rev() {
            if self.y_gc_distribution[i]
                > self.y_gc_distribution[first_mode]
                    - (self.y_gc_distribution[first_mode] / 10 as f64)
            {
                mode += i as f64;
                mode_duplicate += 1;
            } else {
                fell_off_bottom = false;
                break;
            }
        }

        if fell_off_bottom || fell_off_top {
            // If the distribution is so skewed that 95% of the mode
            // is off the 0-100% scale then we keep the mode as the
            // centre of the model
            mode = first_mode as f64;
        } else {
            mode /= mode_duplicate as f64;
        }

        // We can now work out a theoretical distribution
        let mut stdev: f64 = 0.0;

        for i in 0..self.y_gc_distribution.len() {
            stdev += (i as f64 - mode).powi(2) * self.y_gc_distribution[i];
        }

        stdev /= total_count - 1.0;
        stdev = stdev.sqrt();

        let nd = NormalDistribution::new(mode, stdev);

        self.deviation_percent = 0.0;
        for i in 0..self.y_theo_distribution.len() {
            let probability = nd.get_zscore_for_values(i as f64);
            self.y_theo_distribution[i] = probability * total_count;

            if self.y_theo_distribution[i] > self.max {
                self.max = self.y_theo_distribution[i];
            }

            self.deviation_percent +=
                (self.y_theo_distribution[i] - self.y_gc_distribution[i]).abs();
        }

        self.deviation_percent /= total_count;
        self.deviation_percent *= 100.0;
    }

    /// Some data values should be calculated after all sequences have been processed
    fn finish(&mut self) {
        self.calculate_distribution();
    }

    pub fn merge(&mut self, other: &PerSeqGCContent) {}
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct PerBaseNContent {
    n_counts: Vec<usize>,
    not_n_counts: Vec<usize>,
    percentages: Vec<f64>,
    x_categories: Vec<String>,
}

impl PerBaseNContent {
    pub fn new() -> PerBaseNContent {
        return PerBaseNContent {
            n_counts: vec![],
            not_n_counts: vec![],
            percentages: vec![],
            x_categories: vec![],
        };
    }

    pub fn process_sequence(&mut self, record: &impl Record) {
        let seq = record.seq();
        let seq_len = seq.len();
        let n_counts_len = self.n_counts.len();
        if n_counts_len < seq_len {
            // We need to expand the size of the data structures
            for _ in n_counts_len..seq_len {
                self.n_counts.push(0);
                self.not_n_counts.push(0);
            }
        }

        for i in 0..seq_len {
            let base_char = seq[i] as char;
            if base_char == 'N' {
                self.n_counts[i] += 1;
            } else {
                self.not_n_counts[i] += 1;
            }
        }
    }

    fn get_percentages(&mut self) {
        let groups: Vec<BaseGroup> = BaseGroup::make_base_groups(self.n_counts.len());
        let groups_len = groups.len();

        self.x_categories = vec!["".to_string(); groups_len];
        self.percentages = vec![0.0; groups_len];

        let mut total: usize;
        let mut n_count: usize;

        for i in 0..groups_len {
            self.x_categories[i] = groups[i].name();

            n_count = 0;
            total = 0;

            for bp in (groups[i].lower_count() - 1)..groups[i].upper_count() {
                n_count += self.n_counts[bp];
                total += self.n_counts[bp];
                total += self.not_n_counts[bp];
            }

            self.percentages[i] = (n_count as f64) / (total as f64) * 100.0;
        }
    }

    /// Some data values should be calculated after all sequences have been processed
    fn finish(&mut self) {
        self.get_percentages();
    }

    pub fn merge(&mut self, other: &PerBaseNContent) {
        let this_len = self.n_counts.len();
        let other_len = other.n_counts.len();

        if this_len < other_len {
            for i in this_len..other_len {
                self.n_counts.push(0);
                self.not_n_counts.push(0);
            }
        }

        for i in 0..other_len {
            self.n_counts[i] += other.n_counts[i];
            self.not_n_counts[i] += other.not_n_counts[i];
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct SeqLenDistribution {
    len_counts: Vec<usize>,
    graph_counts: Vec<f64>,
    x_categories: Vec<String>,
    max: usize,
}

impl SeqLenDistribution {
    pub fn new() -> SeqLenDistribution {
        return SeqLenDistribution {
            len_counts: vec![],
            graph_counts: vec![],
            x_categories: vec![],
            max: 0,
        };
    }

    pub fn process_sequence(&mut self, record: &impl Record) {
        let seq_len = record.seq().len();
        if seq_len + 2 > self.len_counts.len() {
            for _ in self.len_counts.len()..seq_len + 2 {
                self.len_counts.push(0);
            }
        }
        self.len_counts[seq_len] += 1;
    }

    fn get_size_distribution(&mut self, min: usize, max: usize) -> Vec<usize> {
        // We won't group if they've asked us not to
        // some codes haven't completed

        let mut base = 1;
        while base > (max - min) {
            base /= 10;
        }

        let mut interval: usize = 1;
        let starting: usize;
        let divisions = vec![1, 2, 5];

        'outer: loop {
            for d in 0..divisions.len() {
                let tester = base * divisions[d];
                if (max - min) / tester <= 50 {
                    interval = tester;
                    break 'outer;
                }
            }
            base *= 10;
        }

        let basic_division = (min as f64 / interval as f64).round();
        let test_start = basic_division as usize * interval;
        starting = test_start;

        return vec![starting, interval];
    }

    fn calculate_distribution(&mut self) {
        let mut max_len: isize = 0;
        let mut min_len: isize = -1;
        self.max = 0;

        // Find the min and max lengths
        for i in 0..self.len_counts.len() {
            if self.len_counts[i] > 0 {
                if min_len < 0 {
                    min_len = i as isize;
                }
                max_len = i as isize;
            }
        }

        // We can get a -1 value for min if there aren't any valid sequences
        // at all.
        if min_len < 0 {
            min_len = 0;
        }

        // We put one extra category either side of the actual size
        if min_len > 0 {
            min_len -= 1;
        }
        max_len += 1;

        let start_and_interval = self.get_size_distribution(min_len as usize, max_len as usize);

        // Work out how many categories we need
        let mut categories_counts: usize = 0;
        let mut current_value = start_and_interval[0];
        while current_value <= max_len as usize {
            categories_counts += 1;
            current_value += start_and_interval[1];
        }

        self.graph_counts = vec![0.0; categories_counts];
        self.x_categories = vec!["".to_string(); categories_counts];

        for i in 0..self.graph_counts.len() {
            let min_val = start_and_interval[0] + (start_and_interval[1] * i);
            let mut max_val = start_and_interval[0] + start_and_interval[1] * (i + 1) - 1;

            if max_val > max_len as usize {
                max_val = max_len as usize;
            }

            for bp in min_val..max_val + 1 {
                if bp < self.len_counts.len() {
                    self.graph_counts[i] += self.len_counts[bp] as f64;
                }
            }

            self.x_categories[i] = if start_and_interval[1] == 1 {
                format!("{}", min_val)
            } else {
                format!("{}-{}", min_val, max_val)
            };

            if self.graph_counts[i] as usize > self.max {
                self.max = self.graph_counts[i] as usize;
            }
        }
    }

    /// Some data values should be calculated after all sequences have been processed
    fn finish(&mut self) {
        self.calculate_distribution();
    }

    pub fn merge(&mut self, other: &SeqLenDistribution) {
        let this_len = self.len_counts.len();
        let other_len = other.len_counts.len();

        if this_len < other_len {
            for _ in this_len..other_len {
                self.len_counts.push(0);
            }
        }

        for i in 0..other_len {
            self.len_counts[i] += other.len_counts[i];
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Contaminant {
    name: String,
    forward: Vec<u8>,
    reverse: Vec<u8>,
}

impl Contaminant {
    pub fn new(_name: String, _sequence: String) -> Contaminant {
        let sequence = _sequence.to_uppercase();
        let _forward = _sequence.into_bytes();
        let mut _reverse: Vec<u8> = vec![0; sequence.len()];
        for i in 0.._forward.len() {
            let rev_pos = _reverse.len() - 1 - i;
            let base_char = _forward[i] as char;
            match base_char {
                'G' => {
                    _reverse[rev_pos] = 'C' as u8;
                }
                'A' => {
                    _reverse[rev_pos] = 'T' as u8;
                }
                'T' => {
                    _reverse[rev_pos] = 'A' as u8;
                }
                'C' => {
                    _reverse[rev_pos] = 'G' as u8;
                }
                _ => {
                    panic!(
                        "Contaminant contained the illegal character '{}'",
                        base_char
                    );
                }
            }
        }
        return Contaminant {
            name: _name,
            forward: _forward,
            reverse: _reverse,
        };
    }

    pub fn name(&self) -> String {
        return self.name.clone();
    }

    pub fn find_match(&self, query: &String) -> Option<ContaminantHit> {
        let query = query.to_uppercase();
        let length = query.len();

        // We have a special case for queries between 8 - 20bp where we will allow a hit
        // if it's an exact substring of this contaminant
        if length < 20 && length >= 8 {
            let forward_string = from_utf8(&self.forward).unwrap();
            let reverse_string = from_utf8(&self.reverse).unwrap();
            if forward_string.contains(&query) {
                return Some(ContaminantHit::new(
                    self.clone(),
                    FORWARD_TYPE,
                    query.len(),
                    100,
                ));
            }
            if reverse_string.contains(&query) {
                return Some(ContaminantHit::new(
                    self.clone(),
                    REVERSE_TYPE,
                    query.len(),
                    100,
                ));
            }
        }

        let mut best_hit_option: Option<ContaminantHit> = None;

        // We're going to allow only one mismatch and will require
        // a match of at least 20bp to consider this a match at all

        for offset in (0 - (self.forward.len() as isize - 20))..(query.len() as isize - 20) {
            let this_hit_option: Option<ContaminantHit> = self.sub_find_match(
                &self.forward,
                &query.as_bytes().to_vec(),
                offset,
                FORWARD_TYPE,
            );
            if this_hit_option.clone().is_none() {
                continue;
            }
            if best_hit_option.clone().is_none()
                || this_hit_option.clone().unwrap().length()
                    > best_hit_option.clone().unwrap().length()
            {
                best_hit_option = this_hit_option;
            }
        }

        for offset in (0 - (self.forward.len() as isize - 20))..(query.len() as isize - 20) {
            let this_hit_option: Option<ContaminantHit> = self.sub_find_match(
                &self.forward,
                &query.as_bytes().to_vec(),
                offset,
                REVERSE_TYPE,
            );
            if this_hit_option.clone().is_none() {
                continue;
            }
            if best_hit_option.clone().is_none()
                || this_hit_option.clone().unwrap().length()
                    > best_hit_option.clone().unwrap().length()
            {
                best_hit_option = this_hit_option;
            }
        }

        return best_hit_option;
    }

    pub fn sub_find_match(
        &self,
        ca: &Vec<u8>,
        cb: &Vec<u8>,
        offset: isize,
        direction: usize,
    ) -> Option<ContaminantHit> {
        let mut best_hit_option: Option<ContaminantHit> = None;

        let mut mismatch_count = 0;
        let mut start = 0;
        let mut end = 0;

        for i in 0..ca.len() as isize {
            if i + offset < 0 {
                start = i + 1;
                continue;
            }

            if (i + offset) as usize >= cb.len() {
                break;
            }

            if ca[i as usize] == cb[(i + offset) as usize] {
                end = i;
            } else {
                mismatch_count += 1;
                if mismatch_count > 1 {
                    if 1 + end - start > 20 {
                        // That's the end of this match, see if it's worth recording
                        let id = (((1 + (end - start)) - (mismatch_count - 1)) * 100)
                            / (1 + (end - start));
                        if best_hit_option.clone().is_none()
                            || best_hit_option.clone().unwrap().length()
                                < (1 + end as usize - start as usize)
                            || (best_hit_option.clone().unwrap().length()
                                == (1 + (end - start) as usize))
                                && (best_hit_option.clone().unwrap().percent_id() < id as usize)
                        {
                            best_hit_option = Some(ContaminantHit::new(
                                self.clone(),
                                direction,
                                1 + end as usize - start as usize,
                                id as usize,
                            ));
                        }
                    }
                    start = i + 1;
                    end = i + 1;
                    mismatch_count = 0;
                }
            }
        }

        // See if we ended with a match.
        if 1 + (end - start) > 20 {
            let id = (((1 + (end - start)) - mismatch_count) * 100) / (1 + (end - start));
            if best_hit_option.clone().is_none()
                || best_hit_option.clone().unwrap().length() < (1 + end as usize - start as usize)
                || (best_hit_option.clone().unwrap().length() == (1 + (end - start) as usize))
                    && (best_hit_option.clone().unwrap().percent_id() < id as usize)
            {
                best_hit_option = Some(ContaminantHit::new(
                    self.clone(),
                    direction,
                    1 + end as usize - start as usize,
                    id as usize,
                ));
            }
        }
        return best_hit_option;
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct ContaminantHit {
    contaminant: Contaminant,
    direction: usize,
    length: usize,
    percent_id: usize,
}

impl ContaminantHit {
    pub fn new(
        _contaminant: Contaminant,
        _direction: usize,
        _length: usize,
        _percent_id: usize,
    ) -> ContaminantHit {
        if _direction != FORWARD_TYPE && _direction != REVERSE_TYPE {
            panic!("Direction of hit must be FORWARD or REVERSE");
        }
        return ContaminantHit {
            contaminant: _contaminant,
            direction: _direction,
            length: _length,
            percent_id: _percent_id,
        };
    }

    pub fn contaminant(&self) -> Contaminant {
        return self.contaminant.clone();
    }

    pub fn direction(&self) -> usize {
        return self.direction;
    }

    pub fn length(&self) -> usize {
        return self.length;
    }

    pub fn percent_id(&self) -> usize {
        return self.percent_id;
    }

    pub fn value_string(&self) -> String {
        return format!(
            "{} ({}% over {}bp)",
            self.contaminant.name(),
            self.percent_id,
            self.length
        );
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct OverRepresentedSeq {
    seq: String,
    count: usize,
    percentage: f64,
    contaminant_hit: Option<ContaminantHit>,
}

impl OverRepresentedSeq {
    pub fn new(
        _seq: String,
        _count: usize,
        _percentage: f64,
        contaminants: &Vec<Contaminant>,
    ) -> OverRepresentedSeq {
        return OverRepresentedSeq {
            seq: _seq.clone(),
            count: _count,
            percentage: _percentage,
            contaminant_hit: OverRepresentedSeq::find_contaminants_hit(_seq, contaminants),
        };
    }

    pub fn find_contaminants_hit(
        sequences: String,
        contaminants: &Vec<Contaminant>,
    ) -> Option<ContaminantHit> {
        if contaminants.is_empty() {
            panic!("No any contaminants");
        }

        let mut best_hit: Option<ContaminantHit> = None;

        for c in 0..contaminants.len() {
            let this_hit = (contaminants[c]).find_match(&sequences);
            if this_hit.is_none() {
                continue;
            }
            if best_hit.clone().is_none()
                || this_hit.clone().unwrap().length() > best_hit.clone().unwrap().length()
                || (this_hit.clone().unwrap().length() == best_hit.clone().unwrap().length()
                    && this_hit.clone().unwrap().percent_id()
                        > best_hit.clone().unwrap().percent_id())
            {
                best_hit = this_hit;
            }
        }
        return best_hit;
    }

    pub fn seq(&self) -> String {
        return self.seq.clone();
    }

    pub fn count(&self) -> usize {
        return self.count;
    }

    pub fn percentage(&self) -> f64 {
        return self.percentage;
    }

    pub fn contaminant_hit(&self) -> String {
        if self.contaminant_hit.is_none() {
            return "No Hit".to_string();
        } else {
            return self.contaminant_hit.clone().unwrap().value_string();
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct OverRepresentedSeqs {
    sequences: HashMap<String, usize>,
    count: usize,
    overrepresented_seqs: Vec<OverRepresentedSeq>,
    frozen: bool,
    #[serde(skip_serializing)]
    duplication_module: Option<Box<SeqDuplicationLevel>>,
    observation_cut_off: usize,
    unique_seq_count: usize,
    count_at_unique_limit: usize,
    #[serde(skip_serializing)]
    contaminants: Vec<Contaminant>,
}

impl OverRepresentedSeqs {
    pub fn new(contaminants: &String) -> OverRepresentedSeqs {
        return OverRepresentedSeqs {
            sequences: HashMap::new(),
            count: 0,
            overrepresented_seqs: vec![],
            frozen: false,
            duplication_module: None,
            observation_cut_off: 100000,
            unique_seq_count: 0,
            count_at_unique_limit: 0,
            contaminants: OverRepresentedSeqs::make_contaminants_list(contaminants),
        };
    }

    pub fn read_contaminants_file(contaminants_file: &str) -> String {
        let f = match std::fs::File::open(contaminants_file) {
            Ok(f) => f,
            Err(msg) => panic!("Cannot open {} - {}", contaminants_file, msg),
        };

        return OverRepresentedSeqs::read_contaminants_list(f);
    }

    pub fn read_contaminants_list<R: Read>(mut reader: R) -> String {
        let mut contents = String::new();
        reader.read_to_string(&mut contents).unwrap();
        return contents;
    }

    pub fn make_contaminants_list(contaminants: &String) -> Vec<Contaminant> {
        let mut contaminants_vec: Vec<Contaminant> = vec![];
        for s in contaminants.lines() {
            if s.starts_with('#') {
                continue;
            }
            if s.trim().len() == 0 {
                continue;
            }

            let r = Regex::new("\\t+").unwrap();
            let sections: Vec<&str> = r.split(s).collect();
            contaminants_vec.push(Contaminant::new(
                sections[0].trim().to_string(),
                sections[1].trim().to_string(),
            ));
        }

        return contaminants_vec;
    }

    pub fn duplication_level_module(&mut self) -> Option<Box<SeqDuplicationLevel>> {
        return self.duplication_module.clone();
    }

    pub fn seq(&self) -> HashMap<String, usize> {
        return self.sequences.clone();
    }

    pub fn count_at_unique_limit(&self) -> usize {
        return self.count_at_unique_limit;
    }

    pub fn unique_seq_count(&self) -> usize {
        return self.unique_seq_count;
    }
    pub fn count(&self) -> usize {
        return self.count;
    }

    fn calculate_overrepresented_seq(&mut self) {
        // If the duplication module hasn't already done
        // its calculation it needs to do it now before
        // we stomp all over the data
        let new_duplication_module = self.duplication_module.clone();
        let mut new_duplication_obj = *new_duplication_module.unwrap();
        new_duplication_obj.calculate_levels();
        self.duplication_module = Some(Box::new(new_duplication_obj));

        for (seq, seq_count) in self.sequences.clone() {
            let percentage: f64 = seq_count as f64 / self.count as f64 * 100.0;
            if percentage > INDICATOR_CONFIG_OVERREPESENTED_WARN {
                let os: OverRepresentedSeq =
                    OverRepresentedSeq::new(seq.clone(), seq_count, percentage, &self.contaminants);
                self.overrepresented_seqs.push(os);
            }
        }

        // sort self.overrepresented_seqs by count in OverRepresentedSeq
        self.overrepresented_seqs
            .sort_by(|a, b| b.count.cmp(&a.count));
        self.sequences.clear();
    }

    pub fn process_sequence(&mut self, record: &impl Record) {
        self.count += 1;
        let mut seq = record.seq();

        if FASTQC_CONFIG_DUP_LENGTH != 0 {
            seq = &seq[0..FASTQC_CONFIG_DUP_LENGTH];
        } else if seq.len() > 75 {
            seq = &seq[0..50];
        }

        let seq_string: String = from_utf8(seq).unwrap().to_string();
        if self.sequences.contains_key(&seq_string) {
            self.sequences
                .insert(seq_string.clone(), self.sequences[&seq_string] + 1);

            if !self.frozen {
                self.count_at_unique_limit = self.count;
            }
        } else {
            if !self.frozen {
                self.sequences.insert(seq_string.clone(), 1);
                self.unique_seq_count += 1;
                self.count_at_unique_limit = self.count;
                if self.unique_seq_count == self.observation_cut_off {
                    self.frozen = true;
                }
            }
        }
    }

    /// Some data values should be calculated after all sequences have been processed
    fn finish(&mut self) {
        // update duplication_module when self has been updated
        self.duplication_module = Some(Box::new(SeqDuplicationLevel::new(self)));

        self.calculate_overrepresented_seq();
    }

    pub fn merge(&mut self, other: &OverRepresentedSeqs) {
        for (seq, count) in other.sequences.clone() {
            self.count += count;
            if self.sequences.contains_key(&seq) {
                let current_count = self.sequences[&seq] + count;
                self.sequences.insert(seq, current_count);

                if !self.frozen {
                    self.count_at_unique_limit = self.count;
                }
            } else {
                if !self.frozen {
                    self.sequences.insert(seq, count);
                    self.unique_seq_count += 1;
                    self.count_at_unique_limit = self.count;
                    if self.unique_seq_count == self.observation_cut_off {
                        self.frozen = true;
                    }
                }
            }
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct SeqDuplicationLevel {
    #[serde(skip_serializing)]
    overrepresented_module: Box<OverRepresentedSeqs>,
    dedup_percentages: Vec<f64>,
    total_percentages: Vec<f64>,
    #[serde(skip_serializing)]
    max_count: f64,
    percent_diff_seq: f64,
    labels: Vec<String>,
}

impl SeqDuplicationLevel {
    pub fn new(_overrepresented_module: &OverRepresentedSeqs) -> SeqDuplicationLevel {
        return SeqDuplicationLevel {
            overrepresented_module: Box::new((*_overrepresented_module).clone()),
            dedup_percentages: vec![],
            total_percentages: vec![],
            max_count: 0.0,
            percent_diff_seq: 0.0,
            labels: vec![],
        };
    }

    fn get_corrected_count(
        &mut self,
        count_at_limit: usize,
        total_count: usize,
        duplicate_level: usize,
        number_of_observe: usize,
    ) -> f64 {
        if (count_at_limit == total_count) || (total_count - number_of_observe < number_of_observe)
        {
            return number_of_observe as f64;
        }

        let mut pnot_see_at_limit: f64 = 1.0;
        let limit_of_care = 1.0 - (number_of_observe as f64 / (number_of_observe as f64 + 0.01));

        for i in 0..count_at_limit {
            pnot_see_at_limit *=
                ((total_count - i) - duplicate_level) as f64 / (total_count - i) as f64;

            if pnot_see_at_limit < limit_of_care as f64 {
                pnot_see_at_limit = 0.0;
                break;
            }
        }

        let p_see_at_limit = 1.0 - pnot_see_at_limit;
        let true_count = number_of_observe as f64 / p_see_at_limit;
        return true_count;
    }

    fn calculate_levels(&mut self) {
        if self.dedup_percentages.len() != 0 {
            return;
        }

        self.dedup_percentages = vec![0.0; 16];
        self.total_percentages = vec![0.0; 16];

        let mut collated_counts: HashMap<usize, usize> = HashMap::new();

        for (_, this_count) in self.overrepresented_module.seq() {
            if collated_counts.contains_key(&this_count) {
                collated_counts.insert(this_count, collated_counts[&this_count] + 1);
            } else {
                collated_counts.insert(this_count, 1);
            }
        }

        // Now we can correct each of these
        let mut corrected_counts: HashMap<usize, f64> = HashMap::new();
        for (dup_level, count) in collated_counts {
            corrected_counts.insert(
                dup_level,
                self.get_corrected_count(
                    self.overrepresented_module.count_at_unique_limit(),
                    self.overrepresented_module.count(),
                    dup_level,
                    count,
                ),
            );
        }

        // From the corrected counts we can now work out the raw and deduplicated proportions
        let mut dedup_total: f64 = 0.0;
        let mut raw_total: f64 = 0.0;

        for (dup_level, corrected_count) in corrected_counts {
            dedup_total += corrected_count;
            raw_total += corrected_count * dup_level as f64;

            let mut dup_slot = dup_level - 1;

            if dup_slot > 9999 {
                dup_slot = 15;
            } else if dup_slot > 4999 {
                dup_slot = 14;
            } else if dup_slot > 999 {
                dup_slot = 13;
            } else if dup_slot > 499 {
                dup_slot = 12;
            } else if dup_slot > 99 {
                dup_slot = 11;
            } else if dup_slot > 49 {
                dup_slot = 10;
            } else if dup_slot > 9 {
                dup_slot = 9;
            }

            self.dedup_percentages[dup_slot] += corrected_count;
            self.total_percentages[dup_slot] += corrected_count * dup_level as f64;
        }

        self.labels = vec!["".to_string(); 16];
        for i in 0..self.dedup_percentages.len() {
            if i < 9 {
                self.labels[i] = (i + 1).to_string();
            } else if i == 9 {
                self.labels[i] = ">10".to_string();
            } else if i == 10 {
                self.labels[i] = ">50".to_string();
            } else if i == 11 {
                self.labels[i] = ">100".to_string();
            } else if i == 12 {
                self.labels[i] = ">500".to_string();
            } else if i == 13 {
                self.labels[i] = ">1k".to_string();
            } else if i == 14 {
                self.labels[i] = ">5k".to_string();
            } else if i == 15 {
                self.labels[i] = ">10k".to_string();
            }

            self.dedup_percentages[i] /= dedup_total;
            self.total_percentages[i] /= raw_total;
            self.dedup_percentages[i] *= 100.0;
            self.total_percentages[i] *= 100.0;
        }

        self.percent_diff_seq = (dedup_total / raw_total) * 100.0;
        if raw_total == 0.0 as f64 {
            self.percent_diff_seq = 100.0;
        }
    }

    pub fn process_sequence(&mut self, record: &impl Record) {
        // We don't need to do anything since we use
        // the data structure from the overrepresented sequences
        // module.
    }

    /// Some data values should be calculated after all sequences have been processed
    fn finish(&mut self) {
        // it has been finished in over-represented modules
    }

    pub fn merge(&mut self, other: &SeqDuplicationLevel) {
        // We don't need to do anything since we use
        // the data structure from the overrepresented sequences
        // module.
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Adapter {
    name: String,
    sequence: String,
    positions: Vec<usize>,
}

impl Adapter {
    pub fn new(name: String, sequence: String) -> Adapter {
        return Adapter {
            name: name,
            sequence: sequence,
            positions: vec![0],
        };
    }

    pub fn merge(&mut self, other: &Adapter) {
        let this_len = self.positions.len();
        let other_len = other.positions.len();

        if this_len < other_len {
            self.expand_length_to(other_len);
        }

        for i in 0..other_len {
            self.positions[i] += other.positions[i];
        }
    }

    pub fn increment_count(&mut self, position: usize) {
        self.positions[position] += 1;
    }

    pub fn increment_num_count(&mut self, position: usize, number: usize) {
        self.positions[position] += number;
    }

    pub fn get_count_by_position(&self, position: usize) -> usize {
        return self.positions[position];
    }

    pub fn expand_length_to(&mut self, new_length: usize) {
        let old_len = self.positions.len();
        if new_length > old_len {
            for _ in old_len..new_length {
                self.positions.push(self.positions[old_len - 1]);
            }
        }
    }

    pub fn positions(&mut self) -> Vec<usize> {
        return self.positions.clone();
    }

    pub fn sequence(&mut self) -> String {
        return self.sequence.clone();
    }

    pub fn name(&mut self) -> String {
        return self.name.clone();
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct AdapterContent {
    #[serde(skip_serializing)]
    longest_sequence: usize,
    #[serde(skip_serializing)]
    longest_adapter: usize,
    #[serde(skip_serializing)]
    total_count: usize,

    // This is the full set of Kmers to be reported
    #[serde(skip_serializing)]
    adapters: Vec<Adapter>,

    // This is the data for the Kmers which are going to be placed on the graph
    enrichments: Vec<Vec<f64>>,
    labels: Vec<String>,
    x_labels: Vec<String>,
    #[serde(skip_serializing)]
    groups: Vec<BaseGroup>,
}

impl AdapterContent {
    pub fn new(adapters_content: &String) -> AdapterContent {
        let mut adapters: Vec<Adapter> = Vec::new();
        let mut labels: Vec<String> = Vec::new();
        let mut longest_adapter = 0;

        for s in adapters_content.lines() {
            if s.starts_with('#') {
                continue;
            }
            if s.trim().len() == 0 {
                continue;
            }

            let r = Regex::new("\\t+").unwrap();
            let sections: Vec<&str> = r.split(s).collect();
            let mut adapter = Adapter::new(
                sections[0].trim().to_string(),
                sections[1].trim().to_string(),
            );
            adapters.push(adapter.clone());
            labels.push(adapter.name());
            if adapter.sequence().len() > longest_adapter {
                longest_adapter = adapter.sequence().len();
            }
        }

        return AdapterContent {
            longest_sequence: 0,
            longest_adapter: longest_adapter,
            total_count: 0,
            adapters: adapters,
            enrichments: vec![vec![0.0]],
            labels: labels,
            x_labels: vec![],
            groups: vec![],
        };
    }

    pub fn read_adapter_file(adapter_filepath: &str) -> String {
        let f = match std::fs::File::open(adapter_filepath) {
            Ok(f) => f,
            Err(msg) => panic!("Cannot open {} - {}", adapter_filepath, msg),
        };

        return AdapterContent::read_adapter_list(f);
    }

    pub fn read_adapter_list<R: Read>(mut reader: R) -> String {
        let mut contents = String::new();
        reader.read_to_string(&mut contents).unwrap();
        return contents;
    }

    pub fn process_sequence(&mut self, record: &impl Record) {
        self.total_count += 1;
        // We need to be careful about making sure that a sequence is not only longer
        // than we've seen before, but also that the last position we could find a hit
        // is a positive position.

        // If the sequence is longer than it was then we need to expand the storage in
        // all of the adapter objects to account for this.

        let seq_len = record.seq().len();
        if seq_len > self.longest_sequence && seq_len > self.longest_adapter {
            self.longest_sequence = seq_len;
            for a in 0..self.adapters.len() {
                self.adapters[a].expand_length_to(self.longest_sequence - self.longest_adapter + 1);
            }
        }

        // Now we go through all of the Adapters to see where they occur

        for a in 0..self.adapters.len() {
            let index_option = from_utf8(record.seq())
                .unwrap()
                .find(&self.adapters[a].sequence());
            match index_option {
                Some(index) => {
                    for i in index..(self.longest_sequence - self.longest_adapter + 1) {
                        self.adapters[a].increment_count(i);
                    }
                }

                None => {}
            }
        }
    }

    pub fn calculate_enrichment(&mut self) {
        let mut max_len = 0;
        for a in 0..self.adapters.len() {
            if self.adapters[a].positions().len() > max_len {
                max_len = self.adapters[a].positions().len();
            }
        }

        // We'll be grouping together positions later so make up the groups now
        self.groups = BaseGroup::make_base_groups(max_len);

        self.x_labels = vec!["".to_string(); self.groups.len()];
        for i in 0..self.groups.len() {
            self.x_labels[i] = self.groups[i].name();
        }

        self.enrichments = vec![vec![0.0; self.groups.len()]; self.adapters.len()];

        for a in 0..self.adapters.len() {
            let positions = self.adapters[a].positions();

            for g in 0..self.groups.len() {
                let mut p = self.groups[g].lower_count() - 1;
                while p < self.groups[g].upper_count() && p < positions.len() {
                    self.enrichments[a][g] +=
                        (positions[p] as f64 * 100.0) / self.total_count as f64;
                    p += 1;
                }
                self.enrichments[a][g] /= (self.groups[g].upper_count() as f64
                    - self.groups[g].lower_count() as f64)
                    + 1.0;
            }
        }
    }

    /// Some data values should be calculated after all sequences have been processed
    fn finish(&mut self) {
        self.calculate_enrichment();
    }

    pub fn merge(&mut self, other: &AdapterContent) {
        self.total_count += other.total_count;

        let this_longest_len = self.longest_sequence;
        let other_longest_len = other.longest_sequence;

        if other_longest_len > this_longest_len {
            self.longest_sequence = other_longest_len;
            let new_len = self.longest_sequence - self.longest_adapter + 1;
            for a in 0..self.adapters.len() {
                self.adapters[a].expand_length_to(new_len);
                self.adapters[a].merge(&other.adapters[a]);
            }
        } else if other_longest_len < this_longest_len {
            for a in 0..self.adapters.len() {
                self.adapters[a].merge(&other.adapters[a]);
            }
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Kmer {
    sequence: String,
    count: usize,
    lowest_pvalue: f64,
    #[serde(skip_serializing)]
    obs_exp_position: Vec<f64>,
    #[serde(skip_serializing)]
    positions: Vec<usize>,
}

impl Kmer {
    pub fn new(sequence: String, position: usize, seq_len: usize) -> Kmer {
        let mut positions = vec![0; seq_len];
        positions[position] += 1;

        return Kmer {
            sequence: sequence,
            count: 1,
            lowest_pvalue: 0.0,
            obs_exp_position: vec![],
            positions: positions,
        };
    }

    pub fn merge(&mut self, other: &Kmer) {
        if self.positions.len() != other.positions.len() {
            return;
        }

        for i in 0..self.positions.len() {
            self.positions[i] += other.positions[i];
        }

        self.count += other.count;
    }

    pub fn sequence(&mut self) -> String {
        return self.sequence.clone();
    }

    pub fn count(&mut self) -> usize {
        return self.count;
    }

    pub fn positions(&mut self) -> Vec<usize> {
        return self.positions.clone();
    }

    pub fn increment_count(&mut self, position: usize) {
        self.count += 1;
        if position > self.positions.len() {
            for _ in self.positions.len()..(position + 1) {
                self.positions.push(0);
            }
        }
        self.positions[position] += 1;
    }

    pub fn max_obs_exp(&self) -> f64 {
        let mut max: f64 = 0.0;
        for i in 0..self.obs_exp_position.len() {
            if self.obs_exp_position[i] > max {
                max = self.obs_exp_position[i];
            }
        }
        return max;
    }

    pub fn max_position(&self) -> usize {
        let mut max: f64 = 0.0;
        let mut position: usize = 0;
        for i in 0..self.obs_exp_position.len() {
            if self.obs_exp_position[i] > max {
                max = self.obs_exp_position[i];
                position = i + 1;
            }
        }

        if position == 0 {
            print!("No value > 0 for {}", self.sequence);
            position = 1;
        }
        return position;
    }

    pub fn set_obs_exp_positions(&mut self, oe_positions: Vec<f64>) {
        self.obs_exp_position = oe_positions;
    }

    pub fn set_lowest_pvalue(&mut self, pvalue: f64) {
        self.lowest_pvalue = pvalue;
    }

    pub fn obs_exp_position(&self) -> Vec<f64> {
        return self.obs_exp_position.clone();
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct KmerContent {
    #[serde(skip_serializing)]
    kmers: HashMap<String, Kmer>,
    longest_sequence: usize,
    #[serde(skip_serializing)]
    total_kmer_counts: Vec<Vec<usize>>,
    skip_count: usize,
    min_kmer_size: usize,
    max_kmer_size: usize,
    // This is the full set of Kmers to be reported
    enriched_kmers: Vec<Kmer>,
    // This is the data for the Kmers which are going to be placed on the graph
    enrichments: Vec<Vec<f64>>,
    // For the graph we also need to know the scale we need to use on the axes.
    min_gragh_value: f64,
    max_gragh_value: f64,

    x_categories: Vec<String>,
    x_labels: Vec<String>,
    #[serde(skip_serializing)]
    groups: Vec<BaseGroup>,
}

impl KmerContent {
    pub fn new() -> KmerContent {
        let mut min_kmer_size = 7;
        let mut max_kmer_size = 7;
        if FASTQC_CONFIG_KMER_SIZE != 0 {
            min_kmer_size = FASTQC_CONFIG_KMER_SIZE;
            max_kmer_size = FASTQC_CONFIG_KMER_SIZE;
        }
        return KmerContent {
            kmers: HashMap::new(),
            longest_sequence: 0,
            total_kmer_counts: vec![vec![0; max_kmer_size]],
            skip_count: 0,
            max_kmer_size: max_kmer_size,
            min_kmer_size: min_kmer_size,
            enriched_kmers: vec![],
            enrichments: vec![],
            min_gragh_value: 0.0,
            max_gragh_value: 0.0,
            x_categories: vec![],
            x_labels: vec![],
            groups: vec![],
        };
    }

    fn add_kmer_count(&mut self, position: usize, kemer_len: usize, kmer: &String) {
        let total_kmer_counts_len = self.total_kmer_counts.len();
        if position >= total_kmer_counts_len {
            for _ in total_kmer_counts_len..(position + 1) {
                self.total_kmer_counts.push(vec![0; self.max_kmer_size]);
            }
        }

        if kmer.contains('N') {
            return;
        }

        self.total_kmer_counts[position][kemer_len - 1] += 1;
    }

    fn calculate_enrichment(&mut self) {
        /*
         * For each Kmer we work out whether there is a statistically
         * significant deviation in its coverage at any given position
         * compared to its average coverage over all positions.
         */

        self.groups = BaseGroup::make_base_groups(self.longest_sequence - self.min_kmer_size + 1);

        let mut uneven_kemers: Vec<Kmer> = vec![];

        for (_, kmer) in self.kmers.clone() {
            let mut k = kmer.clone();
            let seq: String = k.sequence();
            let count = k.count();

            let mut total_kmer_count: usize = 0;
            // This gets us the total number of Kmers of this type in the whole
            // dataset.
            for i in 0..self.total_kmer_counts.len() {
                total_kmer_count += self.total_kmer_counts[i][seq.len() - 1];
            }

            // This is the expected proportion of all Kmers which have this
            // specific Kmer sequence.  We no longer make any attempt to judge
            // overall enrichment or depletion of this sequence since once you
            // get to longer lengths the distribution isn't flat anyway

            let expected_proportions: f64 = count as f64 / total_kmer_count as f64;

            let mut obs_exp_positions: Vec<f64> = vec![0.0; self.groups.len()];
            let mut binomial_pvalues: Vec<f64> = vec![0.0; self.groups.len()];
            let position_counts = k.positions();

            for g in 0..self.groups.len() {
                // This is a summation of the number of Kmers of this length which
                // fall into this base group
                let mut total_group_count = 0;

                // This is a summation of the number of hit Kmers which fall within
                // this base group.
                let mut total_group_hits = 0;

                let mut p = self.groups[g].lower_count() - 1;
                while p < self.groups[g].upper_count() && p < position_counts.len() {
                    total_group_count += self.total_kmer_counts[p][seq.len() - 1];
                    total_group_hits += position_counts[p];
                    p += 1;
                }

                let predicted = expected_proportions * total_group_count as f64;
                obs_exp_positions[g] = total_group_hits as f64 / predicted;

                // Now we can run a binomial test to see if there is a significant
                // deviation from what we expect given the number of observations we've
                // made

                let bd: Binomial = Binomial::new(total_group_count, expected_proportions as f64);
                if total_group_hits as f64 > predicted {
                    binomial_pvalues[g] = (1.0 - bd.distribution(total_group_hits as f64)) as f64
                        * 4_i32.pow(k.sequence().len() as u32) as f64;
                } else {
                    binomial_pvalues[g] = 1.0;
                }
            }

            k.set_obs_exp_positions(obs_exp_positions.clone());

            // To keep this we need a p-value below 0.01 and an obs/exp above 5 (actual values are log2 transformed)
            let mut lowest_value: f64 = 1.0;

            for i in 0..binomial_pvalues.len() {
                if binomial_pvalues[i] < 0.01 && obs_exp_positions[i] > 5.0 {
                    if binomial_pvalues[i] < lowest_value {
                        lowest_value = binomial_pvalues[i];
                    }
                }
            }

            if lowest_value < 0.01 {
                k.set_lowest_pvalue(lowest_value);
                uneven_kemers.push(k);
            }
        }

        let mut final_kmers = uneven_kemers;

        // We sort by the highest degree of enrichment over the average
        final_kmers.sort_by(|a, b| {
            if b.max_obs_exp() < a.max_obs_exp() {
                Ordering::Less
            } else if b.max_obs_exp() > a.max_obs_exp() {
                Ordering::Greater
            } else {
                Ordering::Equal
            }
        });

        // So we don't end up with stupidly long lists of Kmers in the
        // report we'll only report the top 20
        if final_kmers.len() > 20 {
            let mut short_kmers: Vec<Kmer> = Vec::new();
            for i in 0..20 {
                short_kmers.push(final_kmers[i].clone());
            }
            final_kmers = short_kmers;
        }

        // Now we take the enrichment positions for the top 6 hits and
        // record these so we can plot them on a line graph
        self.enrichments = vec![vec![]; cmp::min(6, final_kmers.len())];
        self.x_labels = vec!["".to_string(); self.enrichments.len()];

        self.x_categories = vec!["".to_string(); self.groups.len()];

        for i in 0..self.x_categories.len() {
            self.x_categories[i] = self.groups[i].name();
        }

        for k in 0..self.enrichments.len() {
            self.enrichments[k] = vec![0.0; self.groups.len()];

            let obs_exp_pos = final_kmers[k].obs_exp_position();
            for g in 0..self.groups.len() {
                self.enrichments[k][g] = obs_exp_pos[g];
                if obs_exp_pos[g] > self.max_gragh_value {
                    self.max_gragh_value = obs_exp_pos[g];
                }
                if obs_exp_pos[g] < self.min_gragh_value {
                    self.min_gragh_value = obs_exp_pos[g];
                }
            }

            self.x_labels[k] = final_kmers[k].sequence();
        }

        self.min_gragh_value = 0.0;
        self.enriched_kmers = final_kmers;

        // Delete the initial data structure so we don't suck up more memory
        // than we have to.
        self.kmers.clear();
    }

    pub fn process_sequence(&mut self, record: &impl Record) {
        /*
         * The processing done by this module is quite intensive so to speed things
         * up we don't look at every sequence.  Instead we take only 2% of the
         * submitted sequences and extrapolate from these to the full set in the file.
         */
        self.skip_count += 1;

        if self.skip_count % 50 != 0 {
            return;
        }

        let seq: String;
        if record.seq().len() > 500 {
            seq = from_utf8(&record.seq()[0..500]).unwrap().to_string();
        } else {
            seq = from_utf8(&record.seq()).unwrap().to_string();
        }

        if seq.len() > self.longest_sequence {
            self.longest_sequence = seq.len();
        }

        // Now we go through all of the Kmers to count these
        for kmer_size in self.min_kmer_size..(self.max_kmer_size + 1) {
            for i in 0..(seq.len() - kmer_size + 1) {
                let kmer: String = seq[i..(i + kmer_size)].to_string();

                if kmer.len() != kmer_size {
                    panic!(
                        "String length {} wasn't the same as the kmer length {}",
                        kmer.len(),
                        kmer_size
                    );
                }
                // Add to the counts before skipping Kmers containing Ns (see
                // explanation in addKmerCount for the reasoning).
                self.add_kmer_count(i, kmer_size, &kmer);

                // Skip Kmers containing N
                if kmer.contains('N') {
                    continue;
                }

                if self.kmers.contains_key(&kmer) {
                    let mut tt: Kmer = self.kmers[&kmer].clone();
                    tt.increment_count(i);
                    self.kmers.insert(kmer, tt);
                } else {
                    self.kmers
                        .insert(kmer.clone(), Kmer::new(kmer, i, seq.len() - kmer_size + 1));
                }
            }
        }
    }

    /// Some data values should be calculated after all sequences have been processed
    fn finish(&mut self) {
        self.calculate_enrichment();
    }

    pub fn merge(&mut self, other: &KmerContent) {
        // merge longest_sequence
        self.longest_sequence = self.longest_sequence.max(other.longest_sequence);

        // merge skip_count
        self.skip_count += other.skip_count;

        // merge kmers
        for (seq, kmer) in other.kmers.clone() {
            if self.kmers.contains_key(&seq) {
                let mut tt: Kmer = self.kmers[&seq].clone();
                tt.merge(&kmer);
                self.kmers.insert(seq.clone(), tt);
            } else {
                self.kmers.insert(seq.clone(), kmer.clone());
            }
        }

        //merge total_kmer_counts
        let this_total_kmer_len = self.total_kmer_counts.len();
        let other_total_kmer_len = other.total_kmer_counts.len();

        if other_total_kmer_len > this_total_kmer_len {
            for _ in this_total_kmer_len..other_total_kmer_len {
                self.total_kmer_counts.push(vec![0; self.max_kmer_size]);
            }
        }

        for kmer_size in self.min_kmer_size..(self.max_kmer_size + 1) {
            for pos in 0..other_total_kmer_len {
                self.total_kmer_counts[pos][kmer_size - 1] +=
                    other.total_kmer_counts[pos][kmer_size - 1];
            }
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct PerTileQualityScore {
    #[serde(skip_serializing)]
    per_tile_quality_counts: HashMap<usize, Vec<QualityCount>>,
    #[serde(skip_serializing)]
    current_length: usize,
    means: Vec<Vec<f64>>,
    x_labels: Vec<String>,
    tiles: Vec<usize>,
    // #[serde(skip_serializing)]
    high: usize,
    // #[serde(skip_serializing)]
    total_count: usize,
    // #[serde(skip_serializing)]
    split_position: isize,
    // #[serde(skip_serializing)]
    max_deviation: f64,
    #[serde(skip_serializing)]
    ignore_in_report: bool,
}

impl PerTileQualityScore {
    pub fn new() -> PerTileQualityScore {
        return PerTileQualityScore {
            per_tile_quality_counts: HashMap::new(),
            current_length: 0,
            means: vec![],
            x_labels: vec![],
            tiles: vec![],
            high: 0,
            total_count: 0,
            split_position: -1,
            max_deviation: 0.0,
            ignore_in_report: false,
        };
    }

    fn get_mean(&self, tile: usize, min_bp: usize, max_bp: usize, offset: usize) -> f64 {
        let mut count: usize = 0;
        let mut total: f64 = 0.0;
        let quality_counts: Vec<QualityCount> = self.per_tile_quality_counts[&tile].clone();

        for i in min_bp - 1..max_bp {
            if quality_counts[i].total_counts() > 0 {
                count += 1;
                total += quality_counts[i].get_mean(offset);
            }
        }

        if count > 0 {
            return total / count as f64;
        }

        return 0.0;
    }

    pub fn process_sequence(&mut self, record: &impl Record) {
        // Check if we can skip counting because the module is being ignored anyway
        if self.total_count == 0 {
            if INDICATOR_CONFIG_TILE_IGNORE > 0 {
                self.ignore_in_report = true;
            }
        }

        // Don't waste time calculating this if we're not going to use it anyway
        if self.ignore_in_report {
            return;
        }

        // Don't bother with sequences with zero length as they don't have any
        // quality information anyway.
        if record.qual().len() == 0 {
            return;
        }

        self.total_count += 1;
        if self.total_count > 10000 && self.total_count % 10 != 0 {
            return;
        }

        // First try to split the id by :
        let mut tile: usize = 0;
        let id_string = from_utf8(record.head()).unwrap().to_string();
        let split_id_array: Vec<&str> = id_string.split(":").collect();

        // analysis tile id
        if self.split_position >= 0 {
            if split_id_array.len() <= self.split_position as usize {
                println!("Can't extract a number - not enough data");
                self.ignore_in_report = true;
                return;
            }
            tile = split_id_array[self.split_position as usize]
                .parse::<usize>()
                .unwrap();
        } else if split_id_array.len() >= 7 {
            self.split_position = 4;
            tile = split_id_array[4].parse::<usize>().unwrap();
        } else if split_id_array.len() >= 5 {
            self.split_position = 2;
            tile = split_id_array[2].parse::<usize>().unwrap();
        } else {
            // We're not going to get a tile out of this
            self.ignore_in_report = true;
            return;
        }

        let qual = record.qual();

        if self.current_length < qual.len() {
            for (this_tile, quality_count) in self.per_tile_quality_counts.clone() {
                let mut quality_count_new = quality_count.clone();
                for _ in quality_count.len()..qual.len() {
                    quality_count_new.push(QualityCount::new());
                }
                self.per_tile_quality_counts
                    .insert(this_tile, quality_count_new);
            }

            self.current_length = qual.len();
        }

        if !self.per_tile_quality_counts.contains_key(&tile) {
            if self.per_tile_quality_counts.len() > 1000 {
                println!("Too many tiles (>1000) so giving up trying to do per-tile qualities since we're probably parsing the file wrongly");
                self.ignore_in_report = true;
                self.per_tile_quality_counts.clear();
                return;
            }
            let quality_counts: Vec<QualityCount> = vec![QualityCount::new(); self.current_length];
            self.per_tile_quality_counts.insert(tile, quality_counts);
        }

        let mut quality_count: Vec<QualityCount> = self.per_tile_quality_counts[&tile].clone();

        for i in 0..qual.len() {
            quality_count[i].add_value(qual[i] as usize);
        }

        // I guess author forgot the steps as follows:
        self.per_tile_quality_counts.insert(tile, quality_count);
    }

    fn calculate_offset(&self) -> Vec<u8> {
        // Works out from the set of chars what is the most
        // likely encoding scale for this file.
        let mut min_char: u8 = 0;
        let mut max_char: u8 = 0;

        // Iterate through the tiles to check them all in case
        // we're dealing with unrepresentative data in the first one.

        for (tile, quality_count) in self.per_tile_quality_counts.clone() {
            println!("{}", tile);
            for q in 0..quality_count.len() {
                if min_char == 0 {
                    max_char = quality_count[q].clone().get_max_char() as u8;
                    min_char = quality_count[q].clone().get_min_char() as u8;
                } else {
                    if (quality_count[q].get_min_char() as u8) < min_char {
                        min_char = quality_count[q].get_min_char() as u8;
                    }
                    if (quality_count[q].get_max_char() as u8) > max_char {
                        max_char = quality_count[q].get_max_char() as u8;
                    }
                }
            }
        }
        let result: Vec<u8> = vec![min_char, max_char];
        return result;
    }

    fn get_percentages(&mut self, offset: usize) {
        let range = self.calculate_offset();
        self.high = range[1] as usize - offset;

        if self.high < 35 {
            self.high = 35;
        }

        let groups = BaseGroup::make_base_groups(self.current_length);

        let mut tile_numbers = self
            .per_tile_quality_counts
            .keys()
            .copied()
            .collect::<Vec<_>>();
        tile_numbers.sort();

        self.tiles = vec![0; tile_numbers.len()];

        for i in 0..tile_numbers.len() {
            self.tiles[i] = tile_numbers[i];
        }

        self.means = vec![vec![0.0; groups.len()]; tile_numbers.len()];
        self.x_labels = vec!["".to_string(); groups.len()];

        for t in 0..tile_numbers.len() {
            for i in 0..groups.len() {
                if t == 0 {
                    self.x_labels[i] = groups[i].name();
                }

                let min_base = groups[i].lower_count();
                let max_base = groups[i].upper_count();
                self.means[t][i] = self.get_mean(tile_numbers[t], min_base, max_base, offset);
            }
        }

        // Now we normalise across each column to see if there are any tiles with unusually
        // high or low quality.
        let mut max_deviation: f64 = 0.0;

        let mut average_qualities_per_group: Vec<f64> = vec![0.0; groups.len()];

        for t in 0..tile_numbers.len() {
            for i in 0..groups.len() {
                average_qualities_per_group[i] += self.means[t][i];
            }
        }

        for i in 0..average_qualities_per_group.len() {
            average_qualities_per_group[i] /= tile_numbers.len() as f64;
        }

        for i in 0..groups.len() {
            for t in 0..tile_numbers.len() {
                self.means[t][i] -= average_qualities_per_group[i];
                if self.means[t][i].abs() > max_deviation {
                    max_deviation = self.means[t][i].abs();
                }
            }
        }

        self.max_deviation = max_deviation;
    }

    /// Some data values should be calculated after all sequences have been processed
    fn finish(&mut self, offset: usize) {
        self.get_percentages(offset);
    }

    pub fn merge(&mut self, other: &PerTileQualityScore) {
        self.total_count += other.total_count;
        if self.total_count > 10000 && self.total_count % 10 != 0 {
            return;
        }

        self.ignore_in_report = self.ignore_in_report && other.ignore_in_report;

        // judge if expand the per_tile_quality_counts
        if self.current_length < other.current_length {
            for (this_tile, quality_count) in self.per_tile_quality_counts.clone() {
                let mut quality_count_new = quality_count.clone();
                for _ in quality_count.len()..other.current_length {
                    quality_count_new.push(QualityCount::new());
                }
                self.per_tile_quality_counts
                    .insert(this_tile, quality_count_new);
            }

            self.current_length = other.current_length;
        }

        // update the value of per_tile_quality_counts
        for (other_tile, other_quality_count) in other.per_tile_quality_counts.clone() {
            if !self.per_tile_quality_counts.contains_key(&other_tile) {
                if self.per_tile_quality_counts.len() > 1000 {
                    self.ignore_in_report = true;
                    self.per_tile_quality_counts.clear();
                    return;
                }
                self.per_tile_quality_counts
                    .insert(other_tile.clone(), other_quality_count.clone());
            } else {
                let mut current_quality_counts = self.per_tile_quality_counts[&other_tile].clone();
                for pos in 0..self.current_length {
                    current_quality_counts[pos].merge(&other_quality_count[pos]);
                }
                self.per_tile_quality_counts
                    .insert(other_tile.clone(), current_quality_counts);
            }
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct FastQC {
    pub basic_stats: BasicStats,
    pub per_base_seq_quality: PerBaseSeqQuality,
    pub per_seq_quality_score: PerSeqQualityScore,
    pub per_base_seq_content: PerBaseSeqContent,
    pub per_seq_gc_content: PerSeqGCContent,
    pub per_base_n_content: PerBaseNContent,
    pub seq_len_distribution: SeqLenDistribution,
    pub overrepresented_seqs: OverRepresentedSeqs,
    pub seq_duplication_level: Option<SeqDuplicationLevel>,
    pub adapter_content: AdapterContent,
    pub kmer_content: KmerContent,
    pub per_tile_quality_score: PerTileQualityScore,
}

impl FastQC {
    pub fn new(contaminants: &String, adapters: &String) -> FastQC {
        return FastQC {
            basic_stats: BasicStats::new(),
            per_base_seq_quality: PerBaseSeqQuality::new(),
            per_seq_quality_score: PerSeqQualityScore::new(),
            per_base_seq_content: PerBaseSeqContent::new(),
            per_seq_gc_content: PerSeqGCContent::new(),
            per_base_n_content: PerBaseNContent::new(),
            seq_len_distribution: SeqLenDistribution::new(),
            overrepresented_seqs: OverRepresentedSeqs::new(contaminants),
            seq_duplication_level: None,
            adapter_content: AdapterContent::new(adapters),
            kmer_content: KmerContent::new(),
            per_tile_quality_score: PerTileQualityScore::new(),
        };
    }

    pub fn update_name(mut self, filename: &str) -> FastQC {
        self.basic_stats.name = filename.to_string();
        self
    }

    /// Finish method is crucial, don't forget it.
    pub fn finish(&mut self) {
        self.basic_stats.finish();

        self.per_base_seq_quality
            .finish(self.basic_stats.phred.offset);

        self.per_seq_quality_score.finish();

        self.per_base_seq_content.finish();

        self.per_seq_gc_content.finish();

        self.per_base_n_content.finish();

        self.seq_len_distribution.finish();

        self.overrepresented_seqs.finish();

        let mut seq_duplication_level = *self
            .overrepresented_seqs
            .duplication_level_module()
            .unwrap();
        self.seq_duplication_level = Some(seq_duplication_level);

        self.adapter_content.finish();

        self.kmer_content.finish();

        // self.per_tile_quality_score
        //     .finish(self.basic_stats.phred.offset);
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
    pub fn process_sequence(&mut self, record: &impl Record) {
        self.basic_stats.process_sequence(record);

        self.per_base_seq_quality.process_sequence(record);

        self.per_seq_quality_score.process_sequence(record);

        self.per_base_seq_content.process_sequence(record);

        self.per_seq_gc_content.process_sequence(record);

        self.per_base_n_content.process_sequence(record);

        self.seq_len_distribution.process_sequence(record);

        self.overrepresented_seqs.process_sequence(record);

        // For seq_duplication_level sharing the same data structure with overrepresented_seqs
        // we don't call seq_duplication_level' function of process_sequence which is empty.

        self.adapter_content.process_sequence(record);

        self.kmer_content.process_sequence(record);

        // self.per_tile_quality_score.process_sequence(record);
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
            self.basic_stats.merge(&i.basic_stats);
            self.per_base_seq_quality.merge(&i.per_base_seq_quality);
            self.per_seq_quality_score.merge(&i.per_seq_quality_score);
            self.per_base_seq_content.merge(&i.per_base_seq_content);
            self.per_seq_gc_content.merge(&i.per_seq_gc_content);
            self.per_base_n_content.merge(&i.per_base_n_content);
            self.seq_len_distribution.merge(&i.seq_len_distribution);
            self.overrepresented_seqs.merge(&i.overrepresented_seqs);
            // For seq_duplication_level sharing the same data structure with overrepresented_seqs
            // we don't call seq_duplication_level' function of merge which is empty.
            self.adapter_content.merge(&i.adapter_content);
            self.kmer_content.merge(&i.kmer_content);
            self.per_tile_quality_score.merge(&i.per_tile_quality_score);
        }

        // Finish method is crucial, don't forget it.
        // self.finish();
    }
}

pub type FilteredFastQC = FastQC;
