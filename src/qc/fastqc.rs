use fastq::{OwnedRecord, Record};
use serde::{Deserialize, Serialize};
use std::{collections::HashMap, vec, sync::PoisonError, f32::consts::{PI, E}};

const SANGER_ENCODING_OFFSET: usize = 32;
const ILLUMINA_1_3_ENCODING_OFFSET: usize = 64;
const SANGER_ILLUMINA_1_9: &str = "Sanger / Illumina 1.9";
const ILLUMINA_1_3: &str = "Illumina 1.3";
const ILLUMINA_1_5: &str = "Illumina 1.5";

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
pub struct PerSeqQualityScore {
    average_score_counts : HashMap<usize, usize>,
    y_category_count : Vec<usize>,
    x_category_quality : Vec<usize>,
    max_counts: usize,
    most_frequent_score : usize,
    lowest_char : usize,
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
    fn process_sequence(&mut self, record: &OwnedRecord) {
        let mut average_quality = 0;
        for c in record.qual.clone() {
            let num = c.clone() as usize;
            if  num < self.lowest_char {
                self.lowest_char = num;
            }
            average_quality += c as usize;
        }

        if record.qual.len() > 0 {
            average_quality = average_quality / record.qual.len() ;

            if self.average_score_counts.contains_key(&average_quality) {
                let mut current_count = self.average_score_counts[&average_quality];
                current_count += 1;
                self.average_score_counts.insert(average_quality, current_count);
            }
            else {
                self.average_score_counts.insert(average_quality, 1);
            }
        }
    }

    fn calculate_distribution(&mut self) {
        let encoding = PhredEncoding::get_fastq_encoding_offset(self.lowest_char);

        let mut raw_scores = self.average_score_counts.keys().copied().collect::<Vec<_>>();
        raw_scores.sort();
        
        self.y_category_count = vec![0;raw_scores[raw_scores.len()-1]-raw_scores[0]+1];
        self.x_category_quality = vec![0;self.y_category_count.len()];

        for i in 0..self.y_category_count.len() {
            self.x_category_quality[i] = (raw_scores[0]+i)-encoding.offset();
            if self.average_score_counts.contains_key(&(raw_scores[0]+i)) {
                self.y_category_count[i] = self.average_score_counts[&(raw_scores[0]+i)];
            }
        }

        for i in 0..raw_scores.len() {
            if self.y_category_count[i] >self.max_counts {
                self.max_counts = self.y_category_count[i];
                self.most_frequent_score = self.x_category_quality[i];
            }
        }

        
    }
}


#[cfg(test)]
mod per_seq_qua_score_tests {
    use super::*;

    #[test]
    fn test_phred_encoding() {
        let read1 = OwnedRecord {
          head: b"some_name".to_vec(),
          seq: b"GTCGCACTGATCTGGGTTAGGCGCGGAGCCGAGGGTTGCACCATTTTTCATTATTGAATGCCAAGATA".to_vec(),
          qual: b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII".to_vec(),
          sep: None,
        };
        let mut tt = PerSeqQualityScore::new();
        tt.process_sequence(&read1);
        println!("{:?}", tt);
        // assert_eq!(phred.name, "Illumina 1.3".to_string());
        // assert_eq!(phred.offset, 33);
    }

}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct PerBaseSeqContent {
    g_counts: Vec<usize>,
    c_counts: Vec<usize>,
    a_counts: Vec<usize>,
    t_counts: Vec<usize>,
    percentages: Vec<Vec<f32>>,
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
        }
    }

    fn get_percentages(&mut self, offset: usize) {
        let groups: Vec<BaseGroup> = BaseGroup::make_base_groups(self.g_counts.len());
        let length = groups.len();

        let mut g_percent = vec![0.0 ;length];
        let mut a_percent = vec![0.0;length];
        let mut t_percent = vec![0.0;length];
        let mut c_percent = vec![0.0;length];
        
        let mut g_count = 0;
        let mut a_count = 0;
        let mut t_count = 0;
        let mut c_count = 0;
        let mut total = 0;
        for i in 0..length {
            self.x_category[i] = groups[i].name();
            g_count = 0;
            a_count = 0;
            t_count = 0;
            c_count = 0;
            total = 0;
            let current_group = &groups[i];
            for bp in current_group.lower_count()-1 .. current_group.upper_count() {
                total += self.g_counts[bp];
                total += self.c_counts[bp];
                total += self.a_counts[bp];
                total += self.t_counts[bp];

                g_count += self.g_counts[bp];
                c_count += self.c_counts[bp];
                a_count += self.a_counts[bp];
                t_count += self.t_counts[bp];
            }

            g_percent[i] = (g_count as f32/total as f32) * 100 as f32;
            a_percent[i] = (a_count as f32/total as f32) * 100 as f32;
            t_percent[i] = (t_count as f32/total as f32) * 100 as f32;
            c_percent[i] = (c_count as f32/total as f32) * 100 as f32;
        }
        
        self.percentages = vec![t_percent, c_percent, a_percent, g_percent];

    }

    fn process_sequence(&mut self, record: &OwnedRecord) {
        let seq = record.seq();
        let seq_len = seq.len();
        let g_counts_len = self.g_counts.len();

        if g_counts_len <  seq_len {
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
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct GCModelValue {
    percentage: usize,
    increment:f32,
}

impl GCModelValue {
    fn new(_percentage:usize,_increment:f32) ->GCModelValue {
        return GCModelValue {
            percentage:_percentage,
            increment:_increment,
        }
    }

    pub fn percentage(&self) -> usize {
        return self.percentage;
    }

    pub fn increment(&self) -> f32 {
        return self.increment;
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct GCModel {
    read_length: usize,
    models:Vec<Vec<GCModelValue>>,
}

impl GCModel {
    fn new () -> GCModel {
        return GCModel {
            read_length: 0,
            models: vec![vec![]],
        }
    }

    fn new_by_len(length:usize) ->GCModel {
        let mut claim_counts =  vec![0;101];
        let read_length = length;
        let mut models = vec![vec![];length+1];

        for pos in 0..length+1 {
            let mut low_count = (pos as f32 - 0.5) as f32;
            let mut high_count = (pos as f32 + 0.5) as f32;

            if low_count < 0.0 {
                low_count = 0.0;   
            }   
            if high_count < 0.0 {
                high_count = 0.0;   
            }
            if high_count > length as f32{
                high_count = length as f32;
            }
            if low_count > length as f32 {
                low_count = length as f32;
            }
            let low_percent = (low_count*100 as f32 /length as f32).round() as usize;
            let high_percent = (high_count*100 as f32 /length as f32).round() as usize;
            for p in low_percent..high_percent+1 {
                claim_counts[p] += 1;
            }
        }

        // We now do a second pass to make up the model using the weightings
		// we calculated previously.

        for pos in 0..length+1 {
            let mut low_count = (pos as f32 - 0.5) as f32;
            let mut high_count = (pos as f32 + 0.5) as f32;
            if low_count < 0.0 {
                low_count = 0.0;   
            }
            if high_count < 0.0 {
                high_count = 0.0;   
            }
            if high_count > length as f32{
                high_count = length as f32;
            }
            if low_count > length as f32 {
                low_count = length as f32;
            }

            let low_percent = (low_count*100 as f32 /length as f32).round() as usize;
            let high_percent = (high_count*100 as f32 /length as f32).round() as usize;

            let mut model_values:Vec<GCModelValue> = Vec::with_capacity(high_percent-low_percent+1);
            for p in low_percent..high_percent+1 {
                model_values.insert(p-low_percent, GCModelValue::new(p
                    , 1 as f32/claim_counts[p] as f32));
            }
            models[pos] = model_values;
        }
        return GCModel { 
            read_length: read_length, 
            models: models, 
        }
    }

    pub fn get_model_values(&self, gc_count:usize) -> &Vec<GCModelValue> {
        return self.models.get(gc_count).unwrap();
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct NormalDistribution {
    mean:f32,
    stdev:f32,
}

impl NormalDistribution {
    pub fn new(_mean:f32, _stdev:f32) -> NormalDistribution {
        return NormalDistribution {
            mean:_mean,
            stdev:_stdev,
        }
    }

    pub fn get_zscore_for_values(&self, value:f32)->f32 {
        let _stdev = self.stdev;
        let lhs = 1.0/(2.0*_stdev*_stdev*PI).sqrt();
        let rhs = E.powf(0.0-(value-self.mean).powf(2.0) / (2.0*_stdev*_stdev));

        return lhs*rhs;
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct PerSeqGCContent {
    x_category: Vec<usize>,
    y_gc_distribution: Vec<f32>,
    y_theo_distribution: Vec<f32>,
    max:f32,
    deviation_percent:f32,
    cached_models:Vec<GCModel>
}

impl PerSeqGCContent {
    fn new() -> PerSeqGCContent {
        return PerSeqGCContent { 
            max: 0.0,
            deviation_percent: 0.0,
            x_category: vec![], 
            y_gc_distribution: vec![0.0;101],
            y_theo_distribution: vec![0.0;101],
            cached_models:Vec::with_capacity(200),
        }
    }

    fn process_sequence(&mut self, record: &OwnedRecord) {
        let seq = self.truncate_sequence(record);
        let this_seq_length = seq.len();
        if this_seq_length ==0 {
            return;
        }

        let mut this_seq_gc_count = 0;
        for i in 0..this_seq_length {
            let base_char = seq[i] as char;
            if base_char=='G' || base_char=='C' {
                this_seq_gc_count += 1;
            }
        }

        let cached_models_len = self.cached_models.len();
        if  cached_models_len <= this_seq_length { 
            for _ in cached_models_len .. this_seq_length {
                self.cached_models.push(GCModel::new());
            }

            match self.cached_models.get(this_seq_length) {
                None => {
                    self.cached_models.push(GCModel::new_by_len(this_seq_length));
                }
                _ =>{}
            }

            let values:&Vec<GCModelValue> =self.cached_models[this_seq_length].get_model_values(this_seq_gc_count);

            for i in 0..values.len() {
                self.y_gc_distribution[values[i].percentage()] += values[i].increment();
            }
        }
    }

    fn truncate_sequence<'a>(&'a mut self, record: &'a OwnedRecord) -> &[u8]{
        let _seq= record.seq();
        let seq_len = _seq.len();
        if seq_len > 1000 {
            let length = (seq_len / 1000) * 1000;
            return &_seq[0..length];
        }
        else if seq_len > 100 {
            let length = (seq_len / 100) * 100;
            return &_seq[0..length];
        }
        return _seq;
    }

    fn calculate_distribution(&mut self) {
        self.max = 0.0;
        self.x_category = vec![0;self.y_gc_distribution.len()];
        let mut total_count :f32  =0.0;
        // We use the mode to calculate the theoretical distribution
		// so that we cope better with skewed distributions.
        let mut first_mode = 0;
        let mut mode_count:f32 = 0.0;

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

        let mut mode:f32 =0.0;
        let mut mode_duplicate = 0;
        let mut fell_off_top = true;

        for i in first_mode..self.y_gc_distribution.len() {
            if self.y_gc_distribution[i] > self.y_gc_distribution[first_mode] - (self.y_gc_distribution[first_mode]/10 as f32) {
                mode += i as f32;
                mode_duplicate += 1;
            }
            else {
                fell_off_top = false;
                break;
            }
        }

        let mut fell_off_bottom = true;

        for i in (0..first_mode).rev() {
            if self.y_gc_distribution[i] > self.y_gc_distribution[first_mode] - (self.y_gc_distribution[first_mode]/10 as f32) {
                mode += i as f32;
                mode_duplicate += 1;
            }
            else {
				fell_off_bottom = false;
				break;
			}
        }

        if fell_off_bottom || fell_off_top {
			// If the distribution is so skewed that 95% of the mode
			// is off the 0-100% scale then we keep the mode as the 
			// centre of the model
			mode = first_mode as f32;
		}
        else {
            mode /= mode_duplicate as f32;
        }

        // We can now work out a theoretical distribution
        let mut stdev:f32 = 0.0;

        for i in 0..self.y_gc_distribution.len() {
            stdev +=  (i-mode as usize).pow(2) as f32 * self.y_gc_distribution[i] as f32;
        }

        stdev /= total_count-1.0;
        stdev = stdev.sqrt();

        let nd = NormalDistribution::new(mode,stdev);

        self.deviation_percent = 0.0;
        for i in 0..self.y_theo_distribution.len() {
            let probability =nd.get_zscore_for_values(i as f32);
            self.y_theo_distribution[i] = probability*total_count;

            if self.y_theo_distribution[i] > self.max {
                self.max = self.y_theo_distribution[i];
            }

            self.deviation_percent += (self.y_theo_distribution[i] - self.y_gc_distribution[i]).abs();
        }

        self.deviation_percent /= total_count;
        self.deviation_percent *= 100.0;
    }

}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct PerBaseNContent {
    n_counts:Vec<usize>,
    not_n_counts:Vec<usize>,
    percentages:Vec<f32>,
    x_categories:Vec<String>,
}

impl PerBaseNContent {
    pub fn new() -> PerBaseNContent {
        return PerBaseNContent {
            n_counts: vec![],
            not_n_counts: vec![],
            percentages: vec![],
            x_categories: vec![],
        }
    }

    pub fn process_sequence(&mut self, record: &OwnedRecord) {
        let seq = record.seq();
        let seq_len = seq.len();
        let n_counts_len = self.n_counts.len();
        if n_counts_len < seq_len {
            // We need to expand the size of the data structures
            for _ in n_counts_len .. seq_len {
                self.n_counts.push(0);
                self.not_n_counts.push(0);
            }
        }

        for i in 0..seq_len {
            let base_char = seq[i] as char;
            if base_char =='N' {
                self.n_counts[i] += 1;
            }
            else {
                self.not_n_counts[i] += 1;
            }
        }
    }

    fn get_percentages(&mut self) {
        let groups: Vec<BaseGroup> = BaseGroup::make_base_groups(self.n_counts.len());
        let groups_len = groups.len();

        self.x_categories = vec!["".to_string();groups_len];
        self.percentages = vec![0.0;groups_len];

        let mut total:usize;
        let mut n_count:usize;

        for i in 0..groups_len {
            self.x_categories[i] = groups[i].name();

            n_count = 0;
            total = 0;

            for bp in (groups[i].lower_count()-1)..groups[i].upper_count() {
                n_count += self.n_counts[bp];
                total += self.n_counts[bp];
                total += self.not_n_counts[bp];
            }

            self.percentages[i] = (n_count as f32)/(total as f32) *100.0;
        }

    }


}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct SeqLenDistribution {
    len_counts:Vec<usize>,
    graph_counts:Vec<f32>,
    x_categories: Vec<String>,
    max:usize,
}

impl SeqLenDistribution {
    pub fn new() -> SeqLenDistribution {
        return SeqLenDistribution {
            len_counts: vec![],
            graph_counts: vec![],
            x_categories: vec![],
            max:0,
        }
    }

    pub fn process_sequence(&mut self,  record: &OwnedRecord) {
        let seq_len = record.seq().len();
        if seq_len+2 > self.len_counts.len() {
            for _ in self.len_counts.len() .. seq_len+2 {
                self.len_counts.push(0);
            }
        }
        self.len_counts[seq_len] += 1;
    }

    fn get_size_distribution(&mut self,min:usize, max:usize) -> Vec<usize> {
        // We won't group if they've asked us not to
        // some codes haven't completed

        let mut base = 1;
        while base > (max-min) {
            base /= 10;
        }

        let mut interval:usize = 1;
        let mut starting:usize;
        let divisions = vec![1,2,5];

        'outer: while true {
            for d in  0.. divisions.len() {
                let tester = base * divisions[d];
                if (max-min) / tester <=50 {
                    interval = tester;
                    break 'outer;
                }
            }
            base *= 10;
        }

        let basic_division = (min as f32/ interval as f32).round();
        let test_start = basic_division as usize * interval;
        starting = test_start;

        return vec![starting, interval];
    }

    fn calculate_distribution (&mut self) {
        let mut max_len:isize = 0;
        let mut min_len:isize = -1;
        self.max = 0;

        // Find the min and max lengths
        for i in 0..self.len_counts.len() {
            if self.len_counts[i] > 0 {
                if min_len < 0 {
                    min_len  = i as isize;
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
        let mut categories_counts:usize = 0;
        let mut current_value = start_and_interval[0];
        while current_value <= max_len as usize {
            categories_counts += 1;
            current_value = start_and_interval[1];
        }

        self.graph_counts = vec![0.0;categories_counts];
        self.x_categories = vec!["".to_string();categories_counts];

        for i in 0..self.graph_counts.len() {
            let mut min_val = start_and_interval[0] + (start_and_interval[1]*i);
            let mut max_val = start_and_interval[0]+ start_and_interval[1]*(i+1) -1;

            if max_val > max_len as usize {
                max_val = max_len as usize;
            }

            for bp in min_val..max_val+1 {
                if bp < self.len_counts.len() {
                    self.graph_counts[i] += self.len_counts[bp]  as f32;
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

    
}

pub struct OverRepresentedSeqs {
    sequences:HashMap<String, usize>,
    count:usize,
    x_categories: Vec<String>,
    max:usize,
}

pub struct SeqDuplicationLevel {
    overrepresented_module: OverRepresentedSeqs,
    dedup_percentages:Vec<f32>,
    total_percentages:Vec<f32>,
    percent_diff_seq:f32,
    labels:Vec<String>,
}

impl SeqDuplicationLevel {
    pub fn new(&mut self, _overrepresented_module:OverRepresentedSeqs) -> SeqDuplicationLevel {
        return SeqDuplicationLevel {
            overrepresented_module:_overrepresented_module,
            dedup_percentages:vec![],
            total_percentages:vec![],
            percent_diff_seq:0.0,
            labels:vec![],
        }
    }

    fn calculate_levels(&mut self) {
        if self.dedup_percentages.len() != 0 {
            return;
        }

        self.dedup_percentages = vec![0.0;16];
        self.total_percentages = vec![0.0;16];

        let  mut collated_counts:HashMap<usize, usize> = HashMap::new();
    }

    fn get_corrected_count(count_at_limit:usize,  total_count:usize, duplicate_level:usize, number_of_observe:usize) -> f32 {
        if (count_at_limit == total_count) || (total_count - number_of_observe < number_of_observe) {
            return number_of_observe as f32;
        }

       let mut pnot_see_at_limit:f32 = 1.0;
       let limit_of_care = 1.0 - (number_of_observe as f32 / (number_of_observe as f32 + 0.01));

       for i in 0..count_at_limit {
            pnot_see_at_limit *=  ((total_count - i) - duplicate_level) as f32 / (total_count - i) as f32;

            if pnot_see_at_limit <limit_of_care as f32 {
                pnot_see_at_limit = 0.0;
                break;
            }
        }

        let p_see_at_limit = 1.0 -pnot_see_at_limit;
        let true_count = number_of_observe as f32 / p_see_at_limit;
        return true_count;
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
}

impl FastQC {
    pub fn new() -> FastQC {
        return FastQC {
            basic_stats: BasicStats::new(),
            per_base_seq_quality: PerBaseSeqQuality::new(),
            per_seq_quality_score: PerSeqQualityScore::new(),
            per_base_seq_content: PerBaseSeqContent::new(),
            per_seq_gc_content: PerSeqGCContent::new(),
            per_base_n_content: PerBaseNContent::new(),
            seq_len_distribution: SeqLenDistribution::new(),
        };
    }

    pub fn update_name(mut self, filename: &str) -> FastQC {
        self.basic_stats.name = filename.to_string();
        self
    }

    /// Finish method is crucial, don't forget it.
    pub fn finish(&mut self) -> &FastQC {
        self.basic_stats.finish();
        self.per_base_seq_quality
            .get_percentages(self.basic_stats.phred.offset);
        self
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
        let mut seq_len = 0;
        for base in record.seq() {
            let base_char = *base as char;
            match base_char {
                // match char::from(base.clone()).to_uppercase().to_string().as_str() {
                'T' => {
                    self.basic_stats.t_count += 1;
                    seq_len += 1;
                }
                'C' => {
                    self.basic_stats.c_count += 1;
                    seq_len += 1;
                }
                'G' => {
                    self.basic_stats.g_count += 1;
                    seq_len += 1;
                }
                'A' => {
                    self.basic_stats.a_count += 1;
                    seq_len += 1;
                }
                'N' => {
                    self.basic_stats.n_count += 1;
                    seq_len += 1;
                }
                _ => {}
            }
        }

        self.basic_stats.total_bases += seq_len;
        self.basic_stats.set_min_len(seq_len);
        self.basic_stats.set_max_len(seq_len);

        self.set_highest_lowest_char(&record.qual);
        self.basic_stats.total_reads += 1;

        self.per_base_seq_quality.process_qual(&record.qual);

        self.per_seq_quality_score.process_sequence(&record);

        self.per_base_seq_content.process_sequence(&record);

        self.per_seq_gc_content.process_sequence(&record);

        self.per_base_n_content.process_sequence(&record);

        self.seq_len_distribution.process_sequence(&record);
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
            self.per_base_seq_quality
                .add_quality_counts(&i.per_base_seq_quality.quality_counts);
        }

        // Finish method is crucial, don't forget it.
        self.finish();
    }
}

pub type FilteredFastQC = FastQC;
