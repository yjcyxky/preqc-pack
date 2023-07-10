use core::panic;
use rust_htslib::bam::Record;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::fmt::format;
use std::ops::Bound::*;
use stderrlog::new;
use store_interval_tree::Interval;
use store_interval_tree::IntervalTree;

#[derive(Debug, Clone)]
pub struct BamRecord {
    record: Record,
    in_region_flag: bool,
}

impl BamRecord {
    pub fn new(_record: Record) -> Self {
        Self {
            record: _record,
            in_region_flag: true,
        }
    }

    pub fn record(&self) -> &Record {
        &self.record
    }

    pub fn in_region_flag(&self) -> bool {
        self.in_region_flag
    }

    pub fn set_in_region_flag(&mut self, flag: bool) {
        self.in_region_flag = flag;
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct OverlapResult {
    is_overlap_flag: bool,
    is_correct_strand_flag: bool,
}

impl OverlapResult {
    pub fn new(_is_overlap_flag: bool, _is_correct_strand_flag: bool) -> Self {
        Self {
            is_overlap_flag: _is_overlap_flag,
            is_correct_strand_flag: _is_correct_strand_flag,
        }
    }

    pub fn is_strand_matches(&self) -> bool {
        self.is_correct_strand_flag
    }

    pub fn is_interval_overlap(&self) -> bool {
        self.is_overlap_flag
    }
}

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq, Eq, Hash)]
pub struct Feature {
    name: String,
    is_positive_strand_flag: bool,
}

impl Feature {
    pub fn new(_name: String, _is_positive_strand_flag: bool) -> Self {
        Self {
            name: _name,
            is_positive_strand_flag: _is_positive_strand_flag,
        }
    }

    pub fn get_name(&self) -> &str {
        &self.name
    }

    pub fn get_is_positive_strand_flag(&self) -> bool {
        self.is_positive_strand_flag
    }
}

#[derive(Debug, Clone)]
pub struct OverlapDetector {
    lhs_buffer: usize,
    rhs_buffer: usize,
    cache: HashMap<String, IntervalTree<usize, HashSet<Feature>>>,
    unique_intervals: HashMap<Interval<usize>, bool>,
}

impl OverlapDetector {
    pub fn new(_lhs_buffer: usize, _rhs_buffer: usize) -> Self {
        Self {
            lhs_buffer: _lhs_buffer,
            rhs_buffer: _rhs_buffer,
            cache: HashMap::new(),
            unique_intervals: HashMap::new(),
        }
    }

    pub fn add_lhs(&mut self, feature: Feature, start_pos: usize, end_pos: usize, seq_name: &str) {
        let seq_id = seq_name.to_string();

        self.cache.entry(seq_id.clone()).or_insert_with(|| {
            let interval_tree = IntervalTree::<usize, HashSet<Feature>>::new();
            interval_tree
        });

        let tree = self.cache.get_mut(&seq_id).unwrap();
        let start = start_pos + self.lhs_buffer;
        let end = end_pos - self.lhs_buffer;
        let mut features = HashSet::new();
        features.insert(feature.clone());

        if start <= end {
            let interval = Interval::new(Included(start), Included(end));
            if self.unique_intervals.contains_key(&interval) {
                for mut entry in tree.query_mut(&interval) {
                    let interval_key = entry.interval().clone();
                    if interval_key == interval {
                        let features_value = entry.value();
                        features_value.insert(feature.clone());
                    } else {
                        continue;
                    }
                }
            } else {
                self.unique_intervals.insert(interval.clone(), true);
                tree.insert(interval, features);
            }
        }
    }

    pub fn get_overlaps(&self, seq_name: &str, read_start: usize, read_end: usize) -> Vec<Feature> {
        let mut matches: Vec<Feature> = vec![];
        let seq_id = seq_name.to_string();

        if !self.cache.contains_key(&seq_id) {
            return matches;
        }
        let tree = self.cache.get(&seq_id).unwrap();
        let start = read_start + self.lhs_buffer;
        let end = read_end - self.lhs_buffer;

        if start <= end {
            let interval = Interval::new(Included(start), Included(end));
            for entry in tree.query(&interval) {
                let features_value = entry.value();
                for feature in features_value.iter() {
                    matches.push(feature.to_owned());
                }
            }
        }

        matches
    }

    pub fn is_overlap_flag(&self, seq_name: &str, read_start: usize, read_end: usize) -> bool {
        let mut is_overlap_flag = false;
        let seq_id = seq_name.to_string();

        if !self.cache.contains_key(&seq_id) {
            return is_overlap_flag;
        }
        let tree = self.cache.get(&seq_id).unwrap();
        let start = read_start + self.lhs_buffer;
        let end = read_end - self.lhs_buffer;

        if start <= end {
            let interval = Interval::new(Included(start), Included(end));
            is_overlap_flag = tree.overlaps(&interval);
        }

        is_overlap_flag
    }
}

#[derive(Debug, Clone)]
pub struct RegionOverlapLookupTable {
    overlap_detector: OverlapDetector,
    sequence_names: HashSet<String>,
    region_count: usize,
}
impl RegionOverlapLookupTable {
    pub fn new() -> Self {
        Self {
            overlap_detector: OverlapDetector::new(0, 0),
            sequence_names: HashSet::new(),
            region_count: 0,
        }
    }

    pub fn put_region(
        &mut self,
        start_pos: usize,
        end_pos: usize,
        seq_name: &str,
        is_positive_strand_flag: bool,
    ) {
        // unique feature name
        let feature_name = format!("region{}", self.region_count + 1);
        self.region_count += 1;

        self.overlap_detector.add_lhs(
            Feature::new(feature_name, is_positive_strand_flag),
            start_pos,
            end_pos,
            seq_name,
        );
        self.sequence_names.insert(seq_name.to_string());
    }

    pub fn overlaps(&self, read_start: usize, read_end: usize, seq_name: &str) -> bool {
        return self
            .overlap_detector
            .is_overlap_flag(seq_name, read_start, read_end);
    }

    pub fn overlaps_with_strand(
        &self,
        read_start: usize,
        read_end: usize,
        seq_name: &str,
        is_positive_strand_flag: bool,
    ) -> OverlapResult {
        let overlaps = self
            .overlap_detector
            .get_overlaps(seq_name, read_start, read_end);

        let mut num_strand_matches = 0;
        let mut num_total_matches = 0;

        for feature in &overlaps {
            num_total_matches += 1;
            if feature.is_positive_strand_flag == is_positive_strand_flag {
                num_strand_matches += 1;
            }
        }
        OverlapResult::new(
            num_total_matches > 0,
            num_total_matches == num_strand_matches && num_strand_matches > 0,
        )
    }
}
