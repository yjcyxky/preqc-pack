use crate::qc::config::bam_config::Constants;
use bit_vec::BitVec;
use chashmap_serde::CHashMap;
use executor_service::{ExecutorService, Executors, Future};
use exitcode::OK;
use futures::select_biased;
use lazy_static::lazy_static;
use linear_map::LinearMap;
use math::stats;
use regex::Regex;
use rust_htslib::bam::record::{Aux, AuxArray, Cigar};
use rust_htslib::bam::{Header, Read, Reader, Record};
use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use std::{
    collections::HashMap,
    f64::consts::{E, PI},
    str::from_utf8,
    vec,
};
use time::{Date, Duration, Instant};

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub enum SkipDuplicatesMode {
    NONE,
    ONLY_MARKED_DUPLICATES,
    ONLY_DETECTED_DUPLICATES,
    BOTH,
}

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub enum LibraryProtocol {
    NonStrandSpecific,
    StrandSpecificForward,
    StrandSpecificReverse,
}

impl LibraryProtocol {
    pub const PROTOCOL_NON_STRAND_SPECIFIC: &'static str = "non-strand-specific";
    pub const PROTOCOL_FORWARD_STRAND: &'static str = "strand-specific-forward";
    pub const PROTOCOL_REVERSE_STRAND: &'static str = "strand-specific-reverse";
    pub const PROTOCOL_UNKNOWN: &'static str = "unknown";

    pub fn to_string(&self) -> &str {
        match *self {
            LibraryProtocol::NonStrandSpecific => LibraryProtocol::PROTOCOL_NON_STRAND_SPECIFIC,
            LibraryProtocol::StrandSpecificForward => LibraryProtocol::PROTOCOL_FORWARD_STRAND,
            LibraryProtocol::StrandSpecificReverse => LibraryProtocol::PROTOCOL_REVERSE_STRAND,
        }
    }

    pub fn get_protocol_by_name(protocol_name: &str) -> LibraryProtocol {
        match protocol_name {
            LibraryProtocol::PROTOCOL_FORWARD_STRAND => LibraryProtocol::StrandSpecificForward,
            LibraryProtocol::PROTOCOL_NON_STRAND_SPECIFIC => LibraryProtocol::NonStrandSpecific,
            LibraryProtocol::PROTOCOL_REVERSE_STRAND => LibraryProtocol::StrandSpecificReverse,
            _ => panic!(
                "Unknown library protocol name: {}\nSupported protocols: {:?}",
                protocol_name,
                LibraryProtocol::get_protocol_names()
            ),
        }
    }

    pub fn get_protocol_names() -> Vec<&'static str> {
        vec![
            LibraryProtocol::PROTOCOL_NON_STRAND_SPECIFIC,
            LibraryProtocol::PROTOCOL_FORWARD_STRAND,
            LibraryProtocol::PROTOCOL_REVERSE_STRAND,
        ]
    }

    pub fn get_protocol_names_string() -> String {
        format!(
            "{}, {} or {}",
            LibraryProtocol::PROTOCOL_FORWARD_STRAND,
            LibraryProtocol::PROTOCOL_REVERSE_STRAND,
            LibraryProtocol::PROTOCOL_NON_STRAND_SPECIFIC
        )
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
struct XYItem {
    x: f64,
    y: f64,
}

impl XYItem {
    fn new(x: f64, y: f64) -> XYItem {
        XYItem { x, y }
    }

    fn get_x(&self) -> f64 {
        self.x
    }

    fn set_x(&mut self, x: f64) {
        self.x = x;
    }

    fn get_y(&self) -> f64 {
        self.y
    }

    fn set_y(&mut self, y: f64) {
        self.y = y;
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
struct XYIntervalItem {
    #[serde(skip_serializing)]
    x_up_deviation: f64,
    #[serde(skip_serializing)]
    x_down_deviation: f64,
    #[serde(skip_serializing)]
    y_up_deviation: f64,
    #[serde(skip_serializing)]
    y_down_deviation: f64,

    xy_item: XYItem,
}

impl XYIntervalItem {
    fn new(
        x: f64,
        x_down_deviation: f64,
        x_up_deviation: f64,
        y: f64,
        y_down_deviation: f64,
        y_up_deviation: f64,
    ) -> Self {
        let xy_item = XYItem::new(x, y);
        XYIntervalItem {
            x_up_deviation,
            x_down_deviation,
            y_up_deviation,
            y_down_deviation,
            xy_item,
        }
    }

    pub fn get_x_up_deviation(&self) -> f64 {
        self.x_up_deviation
    }

    pub fn set_x_up_deviation(&mut self, x_up_deviation: f64) {
        self.x_up_deviation = x_up_deviation;
    }

    pub fn get_x_down_deviation(&self) -> f64 {
        self.x_down_deviation
    }

    pub fn set_x_down_deviation(&mut self, x_down_deviation: f64) {
        self.x_down_deviation = x_down_deviation;
    }

    pub fn get_y_up_deviation(&self) -> f64 {
        self.y_up_deviation
    }

    pub fn set_y_up_deviation(&mut self, y_up_deviation: f64) {
        self.y_up_deviation = y_up_deviation;
    }

    pub fn get_y_down_deviation(&self) -> f64 {
        self.y_down_deviation
    }

    pub fn set_y_down_deviation(&mut self, y_down_deviation: f64) {
        self.y_down_deviation = y_down_deviation;
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
struct XYVector {
    max_value: f64,
    items: Vec<XYIntervalItem>,
}

impl XYVector {
    fn new() -> XYVector {
        XYVector {
            items: vec![],
            max_value: 0.0,
        }
    }

    fn from_arrays(x: &[f64], y: &[f64]) -> XYVector {
        let mut items = Vec::new();
        let mut max_value = 0.0;
        for i in 0..x.len() {
            items.push(XYIntervalItem::new(x[i], 0.0, 0.0, y[i], 0.0, 0.0));
            if y[i] > max_value {
                max_value = y[i];
            }
        }
        XYVector { items, max_value }
    }

    fn from_lists(x: &[f64], y: &[f64]) -> XYVector {
        let mut items = Vec::new();
        let mut max_value = 0.0;
        for i in 0..x.len() {
            items.push(XYIntervalItem::new(x[i], 0.0, 0.0, y[i], 0.0, 0.0));
            if y[i] > max_value {
                max_value = y[i];
            }
        }
        XYVector { items, max_value }
    }

    fn from_lists_with_deviation(
        x: &[f64],
        y: &[f64],
        deviation: &[f64],
        is_deviation: bool,
    ) -> XYVector {
        let mut items = Vec::new();
        let mut max_value = 0.0;
        for i in 0..x.len() {
            if is_deviation {
                items.push(XYIntervalItem::new(
                    x[i],
                    x[i],
                    x[i],
                    y[i],
                    y[i] - deviation[i],
                    y[i] + deviation[i],
                ));
            } else {
                let (up, down) = if y[i] > deviation[i] {
                    (y[i], deviation[i])
                } else {
                    (deviation[i], y[i])
                };
                items.push(XYIntervalItem::new(x[i], x[i], x[i], y[i], down, up));
            }
            if y[i] > max_value {
                max_value = y[i];
            }
        }
        XYVector { items, max_value }
    }

    fn add_item(&mut self, item: XYItem) {
        self.items
            .push(XYIntervalItem::new(item.x, 0.0, 0.0, item.y, 0.0, 0.0));
        if item.y > self.max_value {
            self.max_value = item.y;
        }
    }

    fn get_x_vector(&self) -> Vec<f64> {
        self.items.iter().map(|item| item.xy_item.x).collect()
    }

    fn get_y_vector(&self) -> Vec<f64> {
        self.items.iter().map(|item| item.xy_item.y).collect()
    }

    fn get_size(&self) -> usize {
        self.items.len()
    }

    fn get(&self, index: usize) -> Option<&XYItem> {
        match self.items.get(index) {
            None => None,
            Some(item) => Some(&item.xy_item),
        }
    }

    fn get_max_value(&self) -> f64 {
        self.max_value
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
struct ReadStartsHistogram {
    current_read_start_position: i64,
    read_start_counter: i32,
    read_starts_histogram: Vec<i64>,
}

impl ReadStartsHistogram {
    pub const MAX_READ_STARTS_PER_POSITION: usize = 50;
    fn new() -> Self {
        Self {
            current_read_start_position: -1,
            read_start_counter: 1,
            read_starts_histogram: vec![0; Self::MAX_READ_STARTS_PER_POSITION + 1],
        }
    }

    pub fn update(&mut self, position: i64) -> bool {
        if position == self.current_read_start_position {
            self.read_start_counter += 1;
        } else {
            let hist_pos = if self.read_start_counter < Self::MAX_READ_STARTS_PER_POSITION as i32 {
                self.read_start_counter
            } else {
                Self::MAX_READ_STARTS_PER_POSITION as i32
            };
            self.read_starts_histogram[hist_pos as usize] += 1;
            self.read_start_counter = 1;
            self.current_read_start_position = position;
        }
        self.read_start_counter > 1
    }

    pub fn get_histogram(&self) -> &[i64] {
        &self.read_starts_histogram
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct ContigRecord {
    name: String,
    position: i64,
    relative: i64,
    size: i32,
}

impl ContigRecord {
    pub fn new(name: String, size: i32) -> ContigRecord {
        ContigRecord {
            name: name,
            position: 1,
            relative: 1,
            size: size,
        }
    }

    pub fn new_with_position(name: String, position: i64, size: i32) -> ContigRecord {
        ContigRecord {
            name: name,
            position: position,
            relative: 1,
            size: size,
        }
    }

    pub fn new_with_relative(
        name: String,
        position: i64,
        relative: i64,
        size: i32,
    ) -> ContigRecord {
        ContigRecord {
            name: name,
            position: position,
            relative: relative,
            size: size,
        }
    }

    pub fn start(&self) -> i64 {
        self.position
    }

    pub fn end(&self) -> i64 {
        self.position + i64::from(self.size) - 1
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn set_name(&mut self, name: String) {
        self.name = name;
    }

    pub fn size(&self) -> i32 {
        self.size
    }

    pub fn set_size(&mut self, size: i32) {
        self.size = size;
    }

    pub fn position(&self) -> i64 {
        self.position
    }

    pub fn set_position(&mut self, position: i64) {
        self.position = position;
    }

    pub fn relative(&self) -> i64 {
        self.relative
    }

    pub fn set_relative(&mut self, relative: i64) {
        self.relative = relative;
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct GenomeLocator {
    contigs: Vec<ContigRecord>,
    positions: HashMap<String, i64>,
    total_size: i64,
}

impl GenomeLocator {
    pub fn new() -> Self {
        Self {
            contigs: vec![],
            /// The absolute start position in all contigs for each contig
            positions: HashMap::new(),
            total_size: 0,
        }
    }

    pub fn add_contig(&mut self, name: String, size: i32) {
        let contig = ContigRecord::new_with_position(name.clone(), self.total_size + 1, size);
        self.contigs.push(contig);
        self.positions.insert(name, self.total_size + 1);
        self.total_size += size as i64;
    }

    pub fn get_absolute_coordinates(&self, name: &str, relative: i32) -> i64 {
        if let Some(position) = self.positions.get(name) {
            *position + i64::from(relative) - 1
        } else {
            -1
        }
    }

    pub fn get_contig(&self, contig_id: usize) -> Option<&ContigRecord> {
        return self.contigs.get(contig_id);
    }

    // pub fn get_contig_coordinates(&self, absolute: i64) -> Option<&ContigRecord> {
    // empty contig list
    // if self.contigs.is_empty() {
    //     None
    // }
    // // mega contig
    // else if self.contigs.len() == 1 {
    //     let mut contig = self.contigs[0];
    //     contig.set_relative(absolute);
    //     Some(&contig)
    // }
    // // search contig
    // else {
    //     let mut last = &self.contigs[0];
    //     for contig in self.contigs.iter().skip(1) {
    //         if contig.position() > absolute {
    //             break;
    //         }
    //         last = contig;
    //     }
    //     let mut contig = last.clone();
    //     contig.set_relative(absolute - last.position() + 1);
    //     Some(&contig)
    // }
    // }

    pub fn get_contigs(&self) -> &Vec<ContigRecord> {
        &self.contigs
    }

    pub fn get_total_size(&self) -> i64 {
        self.total_size
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct ChromosomeInfo {
    name: String,
    length: i64,
    num_bases: i64,
    cov_mean: f64,
    cov_std: f64,
}

impl ChromosomeInfo {
    pub fn new(name: String, length: i64, num_bases: i64, cov_mean: f64, cov_std: f64) -> Self {
        ChromosomeInfo {
            name,
            length,
            num_bases,
            cov_mean,
            cov_std,
        }
    }

    pub fn get_name(&self) -> &String {
        &self.name
    }

    pub fn get_length(&self) -> i64 {
        self.length
    }

    pub fn get_num_bases(&self) -> i64 {
        self.num_bases
    }

    pub fn get_cov_mean(&self) -> f64 {
        self.cov_mean
    }

    pub fn get_cov_std(&self) -> f64 {
        self.cov_std
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
struct BamStats {
    name: String,
    source_file: String,

    // globals
    number_of_mapped_bases: i64,
    number_of_sequenced_bases: i64,
    number_of_aligned_bases: i64,
    number_of_reads: i64,
    number_of_mapped_reads: i64,
    number_of_secondary_alignments: i64,
    number_of_paired_reads: i64,
    number_of_singletons: i64,
    number_of_supp_alignments: i64,
    number_of_mapped_first_of_pair: i64,
    number_of_mapped_second_of_pair: i64,

    // regions related
    number_of_mapped_reads_in_regions: i64,
    number_of_paired_reads_in_regions: i64,
    number_of_singletons_in_regions: i64,
    number_of_mapped_first_of_pair_in_regions: i64,
    number_of_mapped_second_of_pair_in_regions: i64,
    num_correct_strand_reads: i64,

    num_mapped_bases_per_window: Vec<i64>,
    coverage_per_window: Vec<i64>,
    coverage_squared_per_window: Vec<i64>,
    window_lengths: Vec<i64>,

    // ***reference params***
    reference_size: i64,
    number_of_reference_contigs: i64,
    in_region_reference_size: i64,
    num_selected_regions: i32,

    // ***sample params***
    // coverage data
    mean_coverage: f64,
    std_coverage: f64,
    coverage_across_reference: Vec<f64>,
    std_coverage_across_reference: Vec<f64>,
    coverage_histogram_map: HashMap<i64, i64>,
    coverage_histogram_cache: Vec<i64>,
    coverage_histogram: XYVector,
    acum_coverage_histogram: XYVector,
    read_starts_histogram: ReadStartsHistogram,
    unique_read_starts_histogram: XYVector,
    balanced_coverage_histogram: XYVector,
    max_coverage_quota: i32,
    coverage_quotes: XYVector,
    balanced_coverage_bar_names: HashMap<i64, String>,
    duplication_rate: f64,

    // quality
    mean_mapping_quality_per_window: f64,
    mapping_quality_across_reference: Vec<f64>,
    mapping_quality_histogram_map: HashMap<i64, i64>,
    mapping_quality_histogram_cache: Vec<i64>,
    mapping_quality_histogram: XYVector,

    // A content
    number_of_as: i64,
    mean_a_content_per_window: f64,
    mean_a_relative_content_per_window: f64,
    mean_a_relative_content: f64,
    a_content_across_reference: Vec<f64>,
    a_relative_content_across_reference: Vec<f64>,

    // C content
    number_of_cs: i64,
    mean_c_content_per_window: f64,
    mean_c_relative_content_per_window: f64,
    mean_c_relative_content: f64,
    c_content_across_reference: Vec<f64>,
    c_relative_content_across_reference: Vec<f64>,

    // T content
    number_of_ts: i64,
    mean_t_content_per_window: f64,
    mean_t_relative_content_per_window: f64,
    mean_t_relative_content: f64,
    t_content_across_reference: Vec<f64>,
    t_relative_content_across_reference: Vec<f64>,

    // G content
    number_of_gs: i64,
    mean_g_content_per_window: f64,
    mean_g_relative_content_per_window: f64,
    mean_g_relative_content: f64,
    g_content_across_reference: Vec<f64>,
    g_relative_content_across_reference: Vec<f64>,

    // N content
    number_of_ns: i64,
    mean_n_content_per_window: f64,
    mean_n_relative_content_per_window: f64,
    mean_n_relative_content: f64,
    n_content_across_reference: Vec<f64>,
    n_relative_content_across_reference: Vec<f64>,

    // GC content
    mean_gc_content: f64,
    mean_gc_content_per_window: f64,
    mean_gc_relative_content_per_window: f64,
    mean_gc_relative_content: f64,
    gc_content_across_reference: Vec<f64>,
    gc_relative_content_across_reference: Vec<f64>,

    // insert size
    mean_insert_size: f64,
    p25_insert_size: i32,
    median_insert_size: i32,
    p75_insert_size: i32,
    std_insert_size: f64,
    insert_size_across_reference: Vec<f64>,
    insert_size_histogram: XYVector,
    insert_size_array: Vec<i32>,
    insert_size_histogram_map: HashMap<i64, i64>,
    insert_size_histogram_cache: Vec<i64>,

    // reads stats
    read_mean_size: f64,
    read_max_size: i32,
    read_min_size: i32,
    num_clipped_reads: i32,
    num_reads_with_insertion: i32,
    num_reads_with_deletion: i32,
    reads_as_data: Vec<i32>,
    reads_cs_data: Vec<i32>,
    reads_gs_data: Vec<i32>,
    reads_ts_data: Vec<i32>,
    reads_ns_data: Vec<i32>,
    reads_clipping_data: Vec<i32>,
    reads_as_histogram: XYVector,
    reads_cs_histogram: XYVector,
    reads_gs_histogram: XYVector,
    reads_ts_histogram: XYVector,
    reads_ns_histogram: XYVector,
    reads_clipping_profile_histogram: XYVector,
    homopolymer_indels_data: Vec<i32>,
    report_non_zero_coverage_only: bool,

    num_detected_duplicate_reads: u64,
    num_estimated_duplicate_reads: u64,
    num_duplicates_skipped: u64,
    skip_duplicates_mode: SkipDuplicatesMode,
    adapted_mean_coverage: f64,
    report_overlapping_read_pairs: bool,
    num_overlapping_read_pairs: u64,
    num_of_intersecting_mapped_bases: u64,

    // chromosome stats
    chromosome_stats: Vec<ChromosomeInfo>,

    // windows
    number_of_windows: i32,
    number_of_processed_windows: i32,
    number_of_initialized_windows: i32,
    window_sizes: Vec<i64>,
    window_starts: Vec<i64>,
    window_ends: Vec<i64>,
    window_names: Vec<String>,

    // reporting
    active_window_reporting: bool,
    // window_report: Option<std::fs::File>,
    active_coverage_reporting: bool,
    // coverage_report: Option<std::fs::File>,
    sum_coverage_squared: i64,
    sum_coverage: i64,
    warnings: std::collections::HashMap<String, String>,
    locator: GenomeLocator,

    // gc content histogram
    gc_content_histogram: Vec<f64>,
    available_genome_gc_content_data: bool,
    sample_count: i64,
    num_insertions: i32,
    num_deletions: i32,
    num_mismatches: i64,
    acum_edit_distance: i64,
}

impl BamStats {
    pub const CACHE_SIZE: i32 = 2000;
    pub const NUM_BINS: i32 = 1000;
    const SMOOTH_DISTANCE: i32 = 0;

    fn new(
        _name: String,
        _locator: GenomeLocator,
        _reference_size: i64,
        _number_of_windows: usize,
    ) -> Self {
        Self {
            name: _name,
            source_file: "".to_string(),

            // globals
            number_of_mapped_bases: 0,
            number_of_sequenced_bases: 0,
            number_of_aligned_bases: 0,
            number_of_reads: 0,
            number_of_mapped_reads: 0,
            number_of_secondary_alignments: 0,
            number_of_paired_reads: 0,
            number_of_singletons: 0,
            number_of_supp_alignments: 0,
            number_of_mapped_first_of_pair: 0,
            number_of_mapped_second_of_pair: 0,

            // regions related
            number_of_mapped_reads_in_regions: 0,
            number_of_paired_reads_in_regions: 0,
            number_of_singletons_in_regions: 0,
            number_of_mapped_first_of_pair_in_regions: 0,
            number_of_mapped_second_of_pair_in_regions: 0,
            num_correct_strand_reads: 0,

            num_mapped_bases_per_window: vec![],
            coverage_per_window: vec![],
            coverage_squared_per_window: vec![],
            window_lengths: vec![],

            // ***reference params***
            reference_size: _reference_size,
            number_of_reference_contigs: 0,
            in_region_reference_size: 0,
            num_selected_regions: 0,

            // ***sample params***
            // coverage data
            mean_coverage: 0.0,
            std_coverage: 0.0,
            coverage_across_reference: Vec::with_capacity(_number_of_windows),
            std_coverage_across_reference: Vec::with_capacity(_number_of_windows),
            coverage_histogram_map: HashMap::with_capacity(_number_of_windows),
            coverage_histogram_cache: vec![0; BamStats::CACHE_SIZE as usize],
            coverage_histogram: XYVector::new(),
            acum_coverage_histogram: XYVector::new(),
            read_starts_histogram: ReadStartsHistogram::new(),
            unique_read_starts_histogram: XYVector::new(),
            balanced_coverage_histogram: XYVector::new(),
            max_coverage_quota: 50,
            coverage_quotes: XYVector::new(),
            balanced_coverage_bar_names: HashMap::new(),
            duplication_rate: 0.0,

            // quality
            mean_mapping_quality_per_window: 0.0,
            mapping_quality_across_reference: Vec::with_capacity(_number_of_windows),
            mapping_quality_histogram_map: HashMap::with_capacity(_number_of_windows),
            mapping_quality_histogram_cache: vec![0; BamStats::CACHE_SIZE as usize],
            mapping_quality_histogram: XYVector::new(),

            // A content
            number_of_as: 0,
            mean_a_content_per_window: 0.0,
            mean_a_relative_content_per_window: 0.0,
            mean_a_relative_content: 0.0,
            a_content_across_reference: Vec::with_capacity(_number_of_windows),
            a_relative_content_across_reference: Vec::with_capacity(_number_of_windows),

            // C content
            number_of_cs: 0,
            mean_c_content_per_window: 0.0,
            mean_c_relative_content_per_window: 0.0,
            mean_c_relative_content: 0.0,
            c_content_across_reference: Vec::with_capacity(_number_of_windows),
            c_relative_content_across_reference: Vec::with_capacity(_number_of_windows),

            // T content
            number_of_ts: 0,
            mean_t_content_per_window: 0.0,
            mean_t_relative_content_per_window: 0.0,
            mean_t_relative_content: 0.0,
            t_content_across_reference: Vec::with_capacity(_number_of_windows),
            t_relative_content_across_reference: Vec::with_capacity(_number_of_windows),

            // G content
            number_of_gs: 0,
            mean_g_content_per_window: 0.0,
            mean_g_relative_content_per_window: 0.0,
            mean_g_relative_content: 0.0,
            g_content_across_reference: Vec::with_capacity(_number_of_windows),
            g_relative_content_across_reference: Vec::with_capacity(_number_of_windows),

            // N content
            number_of_ns: 0,
            mean_n_content_per_window: 0.0,
            mean_n_relative_content_per_window: 0.0,
            mean_n_relative_content: 0.0,
            n_content_across_reference: Vec::with_capacity(_number_of_windows),
            n_relative_content_across_reference: Vec::with_capacity(_number_of_windows),

            // GC content
            mean_gc_content: 0.0,
            mean_gc_content_per_window: 0.0,
            mean_gc_relative_content_per_window: 0.0,
            mean_gc_relative_content: 0.0,
            gc_content_across_reference: Vec::with_capacity(_number_of_windows),
            gc_relative_content_across_reference: Vec::with_capacity(_number_of_windows),

            // insert size
            mean_insert_size: 0.0,
            p25_insert_size: 0,
            median_insert_size: 0,
            p75_insert_size: 0,
            std_insert_size: 0.0,
            insert_size_histogram: XYVector::new(),
            insert_size_across_reference: Vec::with_capacity(_number_of_windows),
            insert_size_histogram_map: HashMap::with_capacity(_number_of_windows),
            insert_size_histogram_cache: vec![0; BamStats::CACHE_SIZE as usize],
            insert_size_array: vec![],

            // reads stats
            read_mean_size: 0.0,
            read_max_size: 0,
            read_min_size: 0,
            num_clipped_reads: 0,
            num_reads_with_insertion: 0,
            num_reads_with_deletion: 0,
            reads_as_data: vec![],
            reads_cs_data: vec![],
            reads_gs_data: vec![],
            reads_ts_data: vec![],
            reads_ns_data: vec![],
            reads_clipping_data: vec![],
            reads_as_histogram: XYVector::new(),
            reads_cs_histogram: XYVector::new(),
            reads_gs_histogram: XYVector::new(),
            reads_ts_histogram: XYVector::new(),
            reads_ns_histogram: XYVector::new(),
            reads_clipping_profile_histogram: XYVector::new(),
            homopolymer_indels_data: vec![0; 6],
            report_non_zero_coverage_only: true,

            num_detected_duplicate_reads: 0,
            num_estimated_duplicate_reads: 0,
            num_duplicates_skipped: 0,
            skip_duplicates_mode: SkipDuplicatesMode::NONE,
            adapted_mean_coverage: 0.0,
            report_overlapping_read_pairs: true,
            num_overlapping_read_pairs: 0,
            num_of_intersecting_mapped_bases: 0,

            // chromosome stats
            chromosome_stats: vec![],

            // windows
            number_of_windows: _number_of_windows as i32,
            number_of_processed_windows: 0,
            number_of_initialized_windows: 0,
            window_sizes: vec![0; _number_of_windows],
            window_starts: vec![0; _number_of_windows],
            window_ends: vec![0; _number_of_windows],
            window_names: vec!["".to_string(); _number_of_windows],

            // reporting
            active_window_reporting: false,
            // window_report: None,
            active_coverage_reporting: false,
            // coverage_report: None,
            sum_coverage_squared: 0,
            sum_coverage: 0,
            warnings: HashMap::new(),
            locator: _locator,

            // gc content histogram
            gc_content_histogram: vec![0.0; BamStats::NUM_BINS as usize + 1],
            available_genome_gc_content_data: false,
            sample_count: 0,
            num_insertions: 0,
            num_deletions: 0,
            num_mismatches: 0,
            acum_edit_distance: 0,
        }
    }

    fn ensure_list_size(array: &mut Vec<i32>, expected_size: usize) {
        let size = array.len();
        if size < expected_size {
            for _ in 0..expected_size - size + 1 {
                // tothink!() may be it's expected_size-size ,original is wrong
                array.push(0);
            }
        }
    }

    fn inc_processed_windows(&mut self) {
        self.number_of_processed_windows += 1;
    }

    fn set_window_references(&mut self, prefix: &str, window_positions: Vec<i64>) {
        for i in 0..self.number_of_windows as usize {
            self.window_names[i] = format!("{}_{}", prefix, i + 1);
            self.window_starts[i] = window_positions[i];
            if i + 1 == self.number_of_windows as usize {
                self.window_ends[i] = self.reference_size;
            } else {
                self.window_ends[i] = window_positions[i + 1] - 1;
            }
            self.window_sizes[i] = self.window_ends[i] - self.window_starts[i] + 1;
        }
    }

    fn set_number_of_reads(&mut self) {}

    fn set_num_duplicates_skipped(&mut self) {}

    fn set_skip_duplicates_mode(&mut self, skip_marked: bool, skip_detected: bool) {
        if skip_marked && skip_detected {
            self.skip_duplicates_mode = SkipDuplicatesMode::BOTH;
        } else if skip_marked {
            self.skip_duplicates_mode = SkipDuplicatesMode::ONLY_MARKED_DUPLICATES;
        } else if skip_detected {
            self.skip_duplicates_mode = SkipDuplicatesMode::ONLY_DETECTED_DUPLICATES;
        }
    }

    fn set_number_of_secondary_alignments(&mut self) {}
    // inside of regions

    fn set_num_selected_regions(&mut self) {}

    fn set_in_region_reference_size(&mut self) {}

    fn set_number_of_mapped_reads_in_regions(&mut self) {}

    fn set_number_of_mapped_first_of_pair_in_regions(&mut self) {}

    fn set_number_of_mapped_second_of_pair_in_regions(&mut self) {}

    fn set_number_of_singletons_in_regions(&mut self) {}

    fn set_number_of_correct_strand_reads(&mut self) {}

    fn set_number_of_intersecting_read_pairs(
        &mut self,
        num_intersecting_reads_pairs: u64,
        num_intersecting_bases: u64,
    ) {
        self.report_overlapping_read_pairs = true;
        self.num_overlapping_read_pairs = num_intersecting_reads_pairs;
        self.num_of_intersecting_mapped_bases = num_intersecting_bases;
    }

    fn increment_initialized_windows(&mut self) {
        self.number_of_initialized_windows += 1;
    }

    pub fn set_source_file(&mut self, _source_file: String) {
        self.source_file = _source_file;
    }

    fn get_error_rate(&self) -> f64 {
        let mut error_rate = 0.0;
        if self.acum_edit_distance > 0 && self.number_of_mapped_bases > 0 {
            error_rate = self.acum_edit_distance as f64 / self.number_of_mapped_bases as f64;
        }
        error_rate
    }

    fn get_homopolymer_indels_fraction(&self) -> f64 {
        1.0 - self.homopolymer_indels_data[5] as f64
            / (self.num_insertions + self.num_deletions) as f64
    }

    fn get_number_of_processed_windows(&self) -> i32 {
        self.number_of_processed_windows
    }

    fn get_number_of_windows(&self) -> i32 {
        self.number_of_windows
    }

    fn get_number_of_initialized_windows(&self) -> i32 {
        self.number_of_initialized_windows
    }

    fn get_window_start(&self, index: usize) -> i64 {
        self.window_starts[index]
    }

    fn get_window_end(&self, index: usize) -> i64 {
        self.window_ends[index]
    }

    fn get_window_name(&self, index: usize) -> String {
        self.window_names[index as usize].clone()
    }

    fn get_current_window_start(&self) -> i64 {
        self.window_starts[self.number_of_processed_windows as usize]
    }

    fn get_current_window_end(&self) -> i64 {
        self.window_ends[self.number_of_processed_windows as usize]
    }

    fn get_current_window_name(&self) -> String {
        self.window_names[self.number_of_processed_windows as usize].clone()
    }

    fn get_reference_size(&self) -> i64 {
        self.reference_size
    }

    fn get_num_indels(&self) -> usize {
        return (self.num_insertions + self.num_deletions) as usize;
    }

    fn get_num_insertions(&self) -> usize {
        return self.num_insertions as usize;
    }

    fn get_num_deletions(&self) -> usize {
        return self.num_deletions as usize;
    }

    fn get_gc_content_histogram(&self) -> XYVector {
        let mut result = XYVector::new();

        let iter_count = Self::NUM_BINS / 100;
        let mut counter = 0;
        let mut acum = 0.0;
        let mut index = 1.0;

        for i in 1..Self::NUM_BINS + 1 {
            counter += 1;
            acum += self.gc_content_histogram[i as usize];

            if counter == iter_count {
                result.add_item(XYItem::new(index, acum));
                counter = 0;
                acum = 0.0;
                index += 1.0;
            }
        }
        result
    }

    fn update_read_start_histogram_and_judge_is_detected_dup(&mut self, position: i64) -> bool {
        let duplicate = self.read_starts_histogram.update(position);
        if duplicate {
            self.num_estimated_duplicate_reads += 1;
        }
        duplicate
    }

    fn update_insert_size_histogram(&mut self, insert_size: i32) {
        if insert_size > 0 {
            self.insert_size_array.push(insert_size);
        }
    }

    fn update_histograms(&mut self, window: &BamGenomeWindow) {
        for i in 0..window.get_coverage_across_reference().len() {
            if window.selected_regions_available_flag && !window.get_region().get(i).unwrap() {
                continue;
            }

            // coverageData
            self.update_coverage_histogram_value(window.get_coverage_across_reference()[i] as i64);
            let quality = window.get_mapping_quality_across_reference()[i];

            if quality != -1 {
                self.update_mapping_histogram_value(quality);
            }
        }
    }

    fn update_coverage_histogram_value(&mut self, key: i64) {
        if key < Self::CACHE_SIZE as i64 {
            self.coverage_histogram_cache[key as usize] += 1;
        } else if !self.coverage_histogram_map.contains_key(&key) {
            self.coverage_histogram_map.insert(key, 1);
        } else {
            let mut value = self.coverage_histogram_map.get(&key).unwrap();
            let update_value = *value + 1;
            self.coverage_histogram_map.insert(key, update_value);
        }
    }

    fn update_mapping_histogram_value(&mut self, key: i64) {
        if key < Self::CACHE_SIZE as i64 {
            self.mapping_quality_histogram_cache[key as usize] += 1;
        } else if !self.coverage_histogram_map.contains_key(&key) {
            self.coverage_histogram_map.insert(key, 1);
        } else {
            let mut value = self.coverage_histogram_map.get(&key).unwrap();
            let update_value = *value + 1;
            self.coverage_histogram_map.insert(key, update_value);
        }
    }

    fn update_insert_size_histogram_value(&mut self, key: i64) {
        if key < Self::CACHE_SIZE as i64 {
            self.insert_size_histogram_cache[key as usize] += 1;
        } else if !self.insert_size_histogram_map.contains_key(&key) {
            self.insert_size_histogram_map.insert(key, 1);
        } else {
            let mut value = self.insert_size_histogram_map.get(&key).unwrap();
            let update_value = *value + 1;
            self.insert_size_histogram_map.insert(key, update_value);
        }
    }

    fn add_read_stats_data(&mut self, read_stats_collector: &ReadStatsCollector) {
        self.add_reads_as_data(read_stats_collector.get_reads_as_content());
        self.add_reads_ts_data(read_stats_collector.get_reads_ts_content());
        self.add_reads_cs_data(read_stats_collector.get_reads_cs_content());
        self.add_reads_gs_data(read_stats_collector.get_reads_gs_content());
        self.add_reads_ns_data(read_stats_collector.get_reads_ns_content());
        self.add_reads_clipping_info(read_stats_collector.get_reads_clipping_info());

        let reads_gc_content = read_stats_collector.get_reads_gc_content();

        for i in 0..reads_gc_content.len() {
            let index = (reads_gc_content[i] * Self::NUM_BINS as f32).floor();
            self.gc_content_histogram[index as usize] += 1.0;
        }

        self.sample_count += reads_gc_content.len() as i64;
        self.num_clipped_reads += read_stats_collector.get_num_clipped_reads();
        self.num_reads_with_insertion += read_stats_collector.get_num_reads_with_insertion();
        self.num_reads_with_deletion += read_stats_collector.get_num_reads_with_deletion();

        let homopolymer_indels = read_stats_collector.get_homopolymer_indels();
        let mut num_homopolymer_indels = 0;
        for i in 0..5 {
            self.homopolymer_indels_data[i] += homopolymer_indels[i];
            num_homopolymer_indels += homopolymer_indels[i];
        }
        self.num_insertions += read_stats_collector.get_num_insertions();
        self.num_deletions += read_stats_collector.get_num_deletions();
        self.num_mismatches += read_stats_collector.get_num_mismatches() as i64;
        self.acum_edit_distance += read_stats_collector.get_edit_distance() as i64;

        self.homopolymer_indels_data[5] +=
            read_stats_collector.get_num_indels() - num_homopolymer_indels;
    }

    fn add_reads_clipping_info(&mut self, reads_clipping_content: &Vec<i32>) {
        Self::ensure_list_size(&mut self.reads_clipping_data, reads_clipping_content.len());
        for i in 0..reads_clipping_content.len() {
            self.reads_clipping_data[i] = self.reads_clipping_data[i] + reads_clipping_content[i];
        }
    }
    fn add_reads_as_data(&mut self, reads_a_content: &Vec<i32>) {
        Self::ensure_list_size(&mut self.reads_as_data, reads_a_content.len());
        for i in 0..reads_a_content.len() {
            self.reads_as_data[i] = self.reads_as_data[i] + reads_a_content[i];
        }
    }
    fn add_reads_ts_data(&mut self, reads_t_content: &Vec<i32>) {
        Self::ensure_list_size(&mut self.reads_ts_data, reads_t_content.len());
        for i in 0..reads_t_content.len() {
            self.reads_ts_data[i] = self.reads_ts_data[i] + reads_t_content[i];
        }
    }
    fn add_reads_cs_data(&mut self, reads_c_content: &Vec<i32>) {
        Self::ensure_list_size(&mut self.reads_cs_data, reads_c_content.len());
        for i in 0..reads_c_content.len() {
            self.reads_cs_data[i] = self.reads_cs_data[i] + reads_c_content[i];
        }
    }
    fn add_reads_gs_data(&mut self, reads_g_content: &Vec<i32>) {
        Self::ensure_list_size(&mut self.reads_gs_data, reads_g_content.len());
        for i in 0..reads_g_content.len() {
            self.reads_gs_data[i] = self.reads_gs_data[i] + reads_g_content[i];
        }
    }
    fn add_reads_ns_data(&mut self, reads_n_content: &Vec<i32>) {
        Self::ensure_list_size(&mut self.reads_ns_data, reads_n_content.len());
        for i in 0..reads_n_content.len() {
            self.reads_ns_data[i] = self.reads_ns_data[i] + reads_n_content[i];
        }
    }

    fn add_window_information(&mut self, window: &BamGenomeWindow) {
        // global
        self.number_of_mapped_bases += window.get_number_of_mapped_bases();
        self.number_of_sequenced_bases += window.get_number_of_sequenced_bases();
        self.number_of_aligned_bases += window.get_number_of_aligned_bases();

        self.window_lengths
            .push(window.get_effective_window_length());
        self.num_mapped_bases_per_window
            .push(window.get_num_mapped_bases());

        /*
         * Sample
         */

        // coverageData across reference
        self.coverage_across_reference
            .push(window.get_mean_coverage());
        self.std_coverage_across_reference
            .push(window.get_std_coverage());
        if window.detailed_flag {
            self.coverage_squared_per_window
                .push(window.get_sum_coverage_sequared());
            self.coverage_per_window.push(window.get_sum_coverage());
            self.sum_coverage_squared += window.get_sum_coverage_sequared();
            self.sum_coverage += window.get_sum_coverage();
            self.update_histograms(window);
        }

        // quality
        self.mapping_quality_across_reference
            .push(window.get_mean_mapping_quality());

        // A
        self.number_of_as += window.get_number_of_as();
        self.a_content_across_reference
            .push(window.get_mean_a_content());
        self.a_relative_content_across_reference
            .push(window.get_mean_a_relative_content());

        // T
        self.number_of_ts += window.get_number_of_ts();
        self.t_content_across_reference
            .push(window.get_mean_t_content());
        self.t_relative_content_across_reference
            .push(window.get_mean_t_relative_content());

        // C
        self.number_of_cs += window.get_number_of_cs();
        self.c_content_across_reference
            .push(window.get_mean_c_content());
        self.c_relative_content_across_reference
            .push(window.get_mean_c_relative_content());

        // G
        self.number_of_gs += window.get_number_of_gs();
        self.g_content_across_reference
            .push(window.get_mean_g_content());
        self.g_relative_content_across_reference
            .push(window.get_mean_g_relative_content());

        // N
        self.number_of_ns += window.get_number_of_ns();
        self.n_content_across_reference
            .push(window.get_mean_n_content());
        self.n_relative_content_across_reference
            .push(window.get_mean_n_relative_content());

        // GC
        self.gc_content_across_reference
            .push(window.get_mean_gc_content());
        self.gc_relative_content_across_reference
            .push(window.get_mean_gc_relative_content());

        // insert size
        self.insert_size_across_reference
            .push(window.get_mean_insert_size());

        // todo!() reporting
        // if(activeWindowReporting) reportWindow(window);
        // if(isInstanceOfBamGenomeWindow){
        // 	if(activeCoverageReporting) {
        //         reportCoverage((BamDetailedGenomeWindow)window);
        //     }
        // }
    }

    pub fn compute_descriptors(&mut self) {
        //TODO: this can be parallel!
        println!("numberOfMappedBases: {}", self.number_of_mapped_bases);
        println!("referenceSize: {}", self.reference_size);
        println!("numberOfSequencedBases: {}", self.number_of_sequenced_bases);
        println!("numberOfAs: {}", self.number_of_as);

        let effective_ref_size = if self.num_selected_regions > 0 {
            self.in_region_reference_size
        } else {
            self.reference_size
        };
        self.mean_coverage = self.number_of_mapped_bases as f64 / effective_ref_size as f64;

        if self.num_of_intersecting_mapped_bases > 0 {
            self.adapted_mean_coverage =
                (self.number_of_mapped_bases - self.num_of_intersecting_mapped_bases as i64) as f64
                    / effective_ref_size as f64;
        }

        // compute average
        if self.number_of_sequenced_bases != 0 {
            let mean_coverage_squared = self.mean_coverage * self.mean_coverage;
            let std_coverage_squared: f64 = self.sum_coverage_squared as f64
                - 2.0 * self.mean_coverage * self.sum_coverage as f64
                + mean_coverage_squared * effective_ref_size as f64;
            self.std_coverage = (std_coverage_squared / effective_ref_size as f64).sqrt();

            // quality
            self.mean_mapping_quality_per_window =
                Self::compute_average(&self.mapping_quality_across_reference);

            // A
            self.mean_a_content_per_window =
                Self::compute_average(&self.a_content_across_reference);
            self.mean_a_relative_content_per_window =
                Self::compute_average(&self.a_relative_content_across_reference);
            self.mean_a_relative_content =
                self.number_of_as as f64 / self.number_of_sequenced_bases as f64 * 100.0;

            // T
            self.mean_t_content_per_window =
                Self::compute_average(&self.t_content_across_reference);
            self.mean_t_relative_content_per_window =
                Self::compute_average(&self.t_relative_content_across_reference);
            self.mean_t_relative_content =
                self.number_of_ts as f64 / self.number_of_sequenced_bases as f64 * 100.0;

            // C
            self.mean_c_content_per_window =
                Self::compute_average(&self.c_content_across_reference);
            self.mean_c_relative_content_per_window =
                Self::compute_average(&self.c_relative_content_across_reference);
            self.mean_c_relative_content =
                self.number_of_cs as f64 / self.number_of_sequenced_bases as f64 * 100.0;

            // G
            self.mean_g_content_per_window =
                Self::compute_average(&self.g_content_across_reference);
            self.mean_g_relative_content_per_window =
                Self::compute_average(&self.g_relative_content_across_reference);
            self.mean_g_relative_content =
                self.number_of_gs as f64 / self.number_of_sequenced_bases as f64 * 100.0;

            // N
            self.mean_n_content_per_window =
                Self::compute_average(&self.n_content_across_reference);
            self.mean_n_relative_content_per_window =
                Self::compute_average(&self.n_relative_content_across_reference);
            self.mean_n_relative_content =
                self.number_of_ns as f64 / self.number_of_sequenced_bases as f64 * 100.0;

            // GC
            self.mean_gc_content_per_window =
                Self::compute_average(&self.gc_content_across_reference);
            self.mean_gc_relative_content_per_window =
                Self::compute_average(&self.gc_relative_content_across_reference);
            self.mean_gc_relative_content = (self.number_of_gs + self.number_of_cs) as f64
                / self.number_of_sequenced_bases as f64
                * 100.0;
        }

        // todo!()
        // // reporting
        // if(activeWindowReporting) {
        //     closeWindowReporting();
        // }

        // if(activeCoverageReporting) {
        //     closeCoverageReporting();
        // }
    }

    pub fn compute_chromosome_stats(
        &mut self,
        locator: &GenomeLocator,
        chromosome_window_index: &Vec<usize>,
    ) {
        let chromosome_count = chromosome_window_index.len();
        let contig_records = locator.get_contigs();

        self.chromosome_stats = Vec::with_capacity(chromosome_count);

        for k in 0..chromosome_count {
            let first_window_index = chromosome_window_index[k];
            let last_window_index = if k + 1 < chromosome_count {
                chromosome_window_index[k + 1] - 1
            } else {
                self.number_of_windows as usize - 1
            };

            // Computing mean
            let mut num_bases = 0;
            let mut length = 0;
            let mut sum_cov = 0;
            let mut sum_cov_squared = 0;
            for i in first_window_index..last_window_index + 1 {
                num_bases += self.num_mapped_bases_per_window[i];
                sum_cov += self.coverage_per_window[i];
                sum_cov_squared += self.coverage_squared_per_window[i];
                length += self.window_lengths[i];
            }

            let contig = &contig_records[k];
            let mut info = ChromosomeInfo::new(contig.name.clone(), 0, 0, 0.0, 0.0);
            if length != 0 {
                info.length = length;
                info.num_bases = num_bases;
                info.cov_mean = num_bases as f64 / length as f64;
                let mean_coverage_squared = info.cov_mean * info.cov_mean;
                let std_coverage_squared: f64 = sum_cov_squared as f64
                    - 2.0 * info.cov_mean * sum_cov as f64
                    + mean_coverage_squared * length as f64;
                info.cov_std = (std_coverage_squared / length as f64).sqrt();
            }

            self.chromosome_stats.push(info);
        }
    }

    pub fn compute_histograms(&mut self) {
        for i in 0..Self::CACHE_SIZE as usize {
            let val = self.mapping_quality_histogram_cache[i];
            if val > 0 {
                self.mapping_quality_histogram_map.insert(i as i64, val);
            }
        }

        self.mapping_quality_histogram =
            Self::compute_vector_histogram(&self.mapping_quality_histogram_map);
        if self.num_selected_regions > 0 {
            self.mean_mapping_quality_per_window =
                Self::compute_mean_val_from_histogram(&self.mapping_quality_histogram);
        }

        for i in 0..Self::CACHE_SIZE as usize {
            let val = self.coverage_histogram_cache[i];
            if val > 0 {
                self.coverage_histogram_map.insert(i as i64, val);
            }
        }

        self.compute_insert_size_histogram();
        self.compute_coverage_histogram();
        self.compute_unique_read_starts_histogram();
        self.compute_gc_content_histogram();
        self.compute_reads_content_histogram();
        self.compute_reads_clipping_profile_histogram();
    }

    pub fn compute_average(array: &Vec<f64>) -> f64 {
        let mut sum = 0.0;
        for i in array {
            sum += *i;
        }
        sum / array.len() as f64
    }

    pub fn compute_vector_histogram(map: &HashMap<i64, i64>) -> XYVector {
        let mut coverages = vec![0.0; map.len()];
        let mut freqs = vec![0.0; map.len()];

        // read keys
        let mut index = 0;
        for (key, value) in map.iter() {
            coverages[index] = *key as f64;
            index += 1;
        }
        // sort coverageData
        coverages.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for i in 0..coverages.len() {
            let key = coverages[i] as i64;
            freqs[i] = map[&key] as f64;
        }

        // fill output
        let mut vector_histogram = XYVector::new();
        for i in 0..coverages.len() {
            vector_histogram.add_item(XYItem::new(coverages[i], freqs[i]));
        }

        vector_histogram
    }

    pub fn compute_mean_val_from_histogram(histogram: &XYVector) -> f64 {
        let mut mean_val = 0.0;
        let mut global = 0.0;
        let mut total = 0.0;
        for i in 0..histogram.get_size() {
            let val = histogram.get(i).unwrap();
            global += val.get_x() * val.get_y();
            total += val.get_y();
        }

        if total > 0.0 {
            mean_val = global / total;
        }
        mean_val
    }

    pub fn compute_insert_size_histogram(&mut self) {
        if self.insert_size_array.is_empty() {
            return;
        }

        self.insert_size_array.sort();
        let size = self.insert_size_array.len();
        let median_index = size / 2;
        let percentile_25_index = size / 4;
        let percentile_75_index = percentile_25_index * 3;

        self.p25_insert_size = self.insert_size_array[percentile_25_index];
        self.median_insert_size = self.insert_size_array[median_index];
        self.p75_insert_size = self.insert_size_array[percentile_75_index];

        self.mean_insert_size = stats::mean(self.insert_size_array.iter());
        self.std_insert_size = stats::standard_deviation(self.insert_size_array.iter(), 0);

        let border = self.insert_size_array[percentile_75_index] as f64 * 2.0;

        for val in self.insert_size_array.clone() {
            if val as f64 <= border {
                self.update_insert_size_histogram_value(val as i64);
            }
        }

        for i in 0..Self::CACHE_SIZE as usize {
            let val = self.insert_size_histogram_cache[i];
            if val > 0 {
                self.insert_size_histogram_map.insert(i as i64, val);
            }
        }

        self.insert_size_histogram =
            Self::compute_vector_histogram(&self.insert_size_histogram_map);
    }

    pub fn compute_coverage_histogram(&mut self) {
        let map = &self.coverage_histogram_map;
        let mut coverages = vec![0.0; map.len()];
        let mut freqs = vec![0.0; map.len()];

        // read keys
        let mut index = 0;
        for (key, value) in map.iter() {
            coverages[index] = *key as f64;
            index += 1;
        }
        // sort coverageData
        coverages.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for i in 0..coverages.len() {
            let key = coverages[i] as i64;
            freqs[i] = map[&key] as f64;
        }

        // fill output
        for i in 0..coverages.len() {
            self.coverage_histogram
                .add_item(XYItem::new(coverages[i], freqs[i]));
        }

        //TODO: what if coverage is constant, for example 1 everywhere? This has code has to be checked

        // compute balanced coverage histogram
        // self.balanced_coverage_histogram
        let max_coverage = coverages[coverages.len() - 1];
        let min_coverage = coverages[0];
        let bin_count = 30;
        let n = (max_coverage - min_coverage).powf(1.0 / bin_count as f64);

        let mut border = min_coverage as i64;
        let mut balanced_coverages = vec![];
        balanced_coverages.push(border);

        let mut k = if coverages.len() > 1 { 1 } else { 0 };

        for i in 0..bin_count + 1 {
            let new_border = (n.powf(i as f64).round() + min_coverage) as i64;
            if new_border > border && new_border >= coverages[k as usize] as i64 {
                balanced_coverages.push(new_border);
                border = new_border;
                while (k as usize) < coverages.len() && new_border >= coverages[k as usize] as i64 {
                    k += 1;
                }
            }
        }

        let mut border_index = 0;
        if coverages.len() == 1 {
            let bar_name = format!("{}", coverages[0]);
            self.balanced_coverage_bar_names.insert(0, bar_name);
            self.balanced_coverage_histogram
                .add_item(XYItem::new(0.0, freqs[0]));
        } else {
            for i in 0..balanced_coverages.len() {
                let coverage = balanced_coverages[i];
                let mut sum = 0.0;
                let prev_index = border_index;
                while border_index < coverages.len() && coverages[border_index] <= coverage as f64 {
                    sum += freqs[border_index];
                    border_index += 1;
                }
                let bar_name = if border_index == prev_index + 1 {
                    format!("{}", coverages[prev_index] as i64)
                } else {
                    format!(
                        "{} - {}",
                        coverages[prev_index] as i64,
                        coverages[border_index - 1] as i64
                    )
                };
                self.balanced_coverage_bar_names.insert(i as i64, bar_name);
                self.balanced_coverage_histogram
                    .add_item(XYItem::new(i as f64, sum));
            }
        }

        // compute acum histogram
        let mut acum = 0.0;
        for i in 0..coverages.len() {
            acum += freqs[i];
            self.acum_coverage_histogram
                .add_item(XYItem::new(freqs[i], acum));
        }

        // coverageData quotes
        let total = acum;
        acum = 0.0;
        for i in 0..self.max_coverage_quota as i64 + 1 {
            if self.coverage_histogram_map.contains_key(&i) {
                acum += self.coverage_histogram_map[&i] as f64 / total * 100.0;
            }
            self.coverage_quotes
                .add_item(XYItem::new(i as f64 + 1.0, (100.0 - acum).max(0.0)));
        }
    }

    pub fn compute_unique_read_starts_histogram(&mut self) {
        let unique_read_starts_array = self.read_starts_histogram.get_histogram();

        let mut num_positions = 0.0;

        for i in 1..ReadStartsHistogram::MAX_READ_STARTS_PER_POSITION + 1 {
            let val = unique_read_starts_array[i];
            num_positions += val as f64;
            self.unique_read_starts_histogram
                .add_item(XYItem::new(i as f64, val as f64));
        }

        self.duplication_rate = (1.0 - unique_read_starts_array[1] as f64 / num_positions) * 100.0;
    }

    pub fn compute_gc_content_histogram(&mut self) {
        //normalize
        //double normalizer = sampleCount - gcContentHistogram[0];
        for i in 0..Self::NUM_BINS + 1 {
            self.gc_content_histogram[i as usize] /= self.sample_count as f64;
        }

        //smooth
        for i in Self::SMOOTH_DISTANCE..Self::NUM_BINS - Self::SMOOTH_DISTANCE {
            let mut res = 0.0;
            for j in 0..2 * Self::SMOOTH_DISTANCE + 1 {
                let index = (i + j - Self::SMOOTH_DISTANCE) as usize;
                res += self.gc_content_histogram[index];
            }
            self.gc_content_histogram[i as usize] = res / (Self::SMOOTH_DISTANCE * 2 + 1) as f64;
        }
    }

    pub fn compute_reads_content_histogram(&mut self) {
        let total_size = self.reads_as_data.len()
            + self.reads_ts_data.len()
            + self.reads_cs_data.len()
            + self.reads_gs_data.len()
            + self.reads_ns_data.len();

        // make sure that we have enough data
        Self::ensure_list_size(&mut self.reads_as_data, self.read_max_size as usize);
        Self::ensure_list_size(&mut self.reads_ts_data, self.read_max_size as usize);
        Self::ensure_list_size(&mut self.reads_cs_data, self.read_max_size as usize);
        Self::ensure_list_size(&mut self.reads_gs_data, self.read_max_size as usize);
        Self::ensure_list_size(&mut self.reads_ns_data, self.read_max_size as usize);

        for i in 0..self.read_max_size {
            let num_a = self.reads_as_data[i as usize] as f64;
            let num_t = self.reads_ts_data[i as usize] as f64;
            let num_c = self.reads_cs_data[i as usize] as f64;
            let num_g = self.reads_gs_data[i as usize] as f64;
            let num_n = self.reads_ns_data[i as usize] as f64;
            let sum = num_a + num_t + num_c + num_g + num_n;

            self.reads_as_histogram
                .add_item(XYItem::new(i as f64, (num_a / sum) * 100.0));
            self.reads_ts_histogram
                .add_item(XYItem::new(i as f64, (num_t / sum) * 100.0));
            self.reads_cs_histogram
                .add_item(XYItem::new(i as f64, (num_c / sum) * 100.0));
            self.reads_gs_histogram
                .add_item(XYItem::new(i as f64, (num_g / sum) * 100.0));
            self.reads_ns_histogram
                .add_item(XYItem::new(i as f64, (num_n / sum) * 100.0));
        }
    }

    pub fn compute_reads_clipping_profile_histogram(&mut self) {
        if self.read_max_size == 0 || !self.clipping_is_present() {
            return;
        }
        Self::ensure_list_size(&mut self.reads_clipping_data, self.read_max_size as usize);

        let mut total_bases_clipped = 0.0;
        for val in &self.reads_clipping_data {
            total_bases_clipped += *val as f64;
        }

        for pos in 0..self.read_max_size {
            let val = (self.reads_clipping_data[pos as usize] as f64 / total_bases_clipped) * 100.0;
            self.reads_clipping_profile_histogram
                .add_item(XYItem::new(pos as f64, val));
        }
    }

    pub fn clipping_is_present(&self) -> bool {
        for val in &self.reads_clipping_data {
            if *val > 0 {
                return true;
            }
        }
        return false;
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct BamGenomeWindow {
    detailed_flag: bool,
    name: String,
    // window params
    start: i64,
    end: i64,
    window_size: i64,
    // general
    number_of_mapped_bases: i64,
    number_of_sequenced_bases: i64,
    number_of_aligned_bases: i64,
    // gff-like
    selected_regions_available_flag: bool,
    #[serde(skip_serializing)]
    #[serde(skip_deserializing)]
    selected_regions: BitVec,
    /*
     *
     * Reference params
     *
     */
    reference_available: bool,
    // A content
    number_of_as_in_reference: i64,
    a_relative_content_in_reference: f64,
    // C content
    number_of_cs_in_reference: i64,
    c_relative_content_in_reference: f64,
    // T content
    number_of_ts_in_reference: i64,
    t_relative_content_in_reference: f64,
    // G content
    number_of_gs_in_reference: i64,
    g_relative_content_in_reference: f64,
    // N content
    number_of_ns_in_reference: i64,
    n_relative_content_in_reference: f64,
    // GC content
    number_of_gcs_in_reference: i64,
    gc_relative_content_in_reference: f64,
    // AT content
    number_of_ats_in_reference: i64,
    at_relative_content_in_reference: f64,

    /*
     *
     * Sample params
     *
     */
    // coverageData
    mean_coverage: f64,
    std_coverage: f64,

    // quality
    acum_mapping_quality: f64,
    mean_mapping_quality: f64,

    // A content
    number_of_as: i64,
    mean_a_content: f64,
    mean_a_relative_content: f64,

    // C content
    number_of_cs: i64,
    mean_c_content: f64,
    mean_c_relative_content: f64,

    // G content
    number_of_ts: i64,
    mean_g_content: f64,
    mean_g_relative_content: f64,

    // T content
    number_of_gs: i64,
    mean_t_content: f64,
    mean_t_relative_content: f64,

    // N content
    number_of_ns: i64,
    mean_n_content: f64,
    mean_n_relative_content: f64,

    // GC content
    mean_gc_content: f64,
    mean_gc_relative_content: f64,

    correct_insert_sizes: i32,
    acum_insert_size: f64,
    mean_insert_size: f64,
    effective_window_length: i64,

    /// Details
    // coverageData
    coverage_across_reference: Vec<i32>,
    // quality
    mapping_quality_across_reference: Vec<i64>,
    // required for calculation of mean and std of coverage
    sum_coverage_squared: i64,
    sum_coverage: i64,
}

impl BamGenomeWindow {
    pub fn new(
        _name: String,
        _start: i64,
        _end: i64,
        _reference: &Vec<u8>,
        _detailed_flag: bool,
    ) -> Self {
        let _window_size = _end - _start + 1;
        let mut instance = Self {
            detailed_flag: _detailed_flag,
            name: _name,
            // window params
            start: _start,
            end: _end,
            window_size: _window_size,
            // general
            number_of_mapped_bases: 0,
            number_of_sequenced_bases: 0,
            number_of_aligned_bases: 0,
            // gff-like
            selected_regions_available_flag: false,
            selected_regions: BitVec::new(),
            /*
             *
             * Reference params
             *
             */
            reference_available: false,
            // A content
            number_of_as_in_reference: 0,
            a_relative_content_in_reference: 0.0,
            // C content
            number_of_cs_in_reference: 0,
            c_relative_content_in_reference: 0.0,
            // T content
            number_of_ts_in_reference: 0,
            t_relative_content_in_reference: 0.0,
            // G content
            number_of_gs_in_reference: 0,
            g_relative_content_in_reference: 0.0,
            // N content
            number_of_ns_in_reference: 0,
            n_relative_content_in_reference: 0.0,
            // GC content
            number_of_gcs_in_reference: 0,
            gc_relative_content_in_reference: 0.0,
            // AT content
            number_of_ats_in_reference: 0,
            at_relative_content_in_reference: 0.0,

            /*
             *
             * Sample params
             *
             */
            // coverageData
            mean_coverage: 0.0,
            std_coverage: 0.0,

            // quality
            acum_mapping_quality: 0.0,
            mean_mapping_quality: 0.0,

            // A content
            number_of_as: 0,
            mean_a_content: 0.0,
            mean_a_relative_content: 0.0,

            // C content
            number_of_cs: 0,
            mean_c_content: 0.0,
            mean_c_relative_content: 0.0,

            // G content
            number_of_ts: 0,
            mean_g_content: 0.0,
            mean_g_relative_content: 0.0,

            // T content
            number_of_gs: 0,
            mean_t_content: 0.0,
            mean_t_relative_content: 0.0,

            // N content
            number_of_ns: 0,
            mean_n_content: 0.0,
            mean_n_relative_content: 0.0,

            // GC content
            mean_gc_content: 0.0,
            mean_gc_relative_content: 0.0,

            correct_insert_sizes: 0,
            acum_insert_size: 0.0,
            mean_insert_size: 0.0,
            effective_window_length: 0,

            // coverageData
            coverage_across_reference: vec![0; _window_size as usize],
            // quality
            mapping_quality_across_reference: vec![0; _window_size as usize],
            // required for calculation of mean and std of coverage
            sum_coverage_squared: 0,
            sum_coverage: 0,
        };
        if !_reference.is_empty() {
            instance.process_reference(_reference);
        }
        instance
    }
    fn process_reference(&mut self, reference: &Vec<u8>) {
        for nucleotide in reference {
            match *nucleotide as char {
                'A' => self.number_of_as_in_reference += 1,
                'C' => self.number_of_cs_in_reference += 1,
                'T' => self.number_of_ts_in_reference += 1,
                'G' => self.number_of_gs_in_reference += 1,
                'N' => self.number_of_ns_in_reference += 1,
                _ => (),
            }
        }
        self.number_of_gcs_in_reference =
            self.number_of_gs_in_reference + self.number_of_cs_in_reference;
        self.number_of_ats_in_reference =
            self.number_of_as_in_reference + self.number_of_ts_in_reference;

        // relative contents
        self.a_relative_content_in_reference =
            ((self.number_of_as_in_reference as f64) / (self.window_size as f64)) * 100.0;
        self.c_relative_content_in_reference =
            ((self.number_of_cs_in_reference as f64) / (self.window_size as f64)) * 100.0;
        self.t_relative_content_in_reference =
            ((self.number_of_ts_in_reference as f64) / (self.window_size as f64)) * 100.0;
        self.g_relative_content_in_reference =
            ((self.number_of_gs_in_reference as f64) / (self.window_size as f64)) * 100.0;
        self.n_relative_content_in_reference =
            ((self.number_of_ns_in_reference as f64) / (self.window_size as f64)) * 100.0;
        self.gc_relative_content_in_reference =
            ((self.number_of_gcs_in_reference as f64) / (self.window_size as f64)) * 100.0;
        self.at_relative_content_in_reference =
            ((self.number_of_ats_in_reference as f64) / (self.window_size as f64)) * 100.0;

        self.reference_available = true;
    }

    pub fn acum_base(&mut self, relative: i64) {
        self.number_of_sequenced_bases += 1;
        self.number_of_mapped_bases += 1;
        if self.detailed_flag {
            self.coverage_across_reference[relative as usize] =
                self.coverage_across_reference[relative as usize] + 1;
        }
    }

    pub fn acum_insert_size(&mut self, insert_size: i64) {
        if insert_size > 0 {
            self.correct_insert_sizes += 1;
            self.acum_insert_size += insert_size.abs() as f64;
        }
    }
    /// Calculate the average metrics
    pub fn compute_descriptors(&mut self) {
        self.effective_window_length = self.window_size;
        if self.selected_regions_available_flag {
            // todo!() self.selected_regions.cardinality()
        }

        if self.effective_window_length != 0 {
            self.mean_coverage =
                self.number_of_mapped_bases as f64 / self.effective_window_length as f64;
        } else {
            self.mean_coverage = 0.0;
        }

        if self.correct_insert_sizes > 0 {
            self.mean_insert_size = self.acum_insert_size / self.correct_insert_sizes as f64;
        }

        // ACTG absolute content
        if self.mean_coverage > 0.0 {
            self.mean_mapping_quality =
                self.acum_mapping_quality / self.number_of_mapped_bases as f64;

            self.mean_a_content = self.number_of_as as f64 / self.mean_coverage;
            self.mean_t_content = self.number_of_ts as f64 / self.mean_coverage;
            self.mean_c_content = self.number_of_cs as f64 / self.mean_coverage;
            self.mean_g_content = self.number_of_gs as f64 / self.mean_coverage;
            self.mean_n_content = self.number_of_ns as f64 / self.mean_coverage;
            self.mean_gc_content = self.mean_g_content + self.mean_c_content;

            // ACTG relative content
            let acum_mean_content = self.mean_a_content
                + self.mean_t_content
                + self.mean_c_content
                + self.mean_g_content
                + self.mean_n_content;
            self.mean_a_relative_content = (self.mean_a_content / acum_mean_content) * 100.0;
            self.mean_t_relative_content = (self.mean_t_content / acum_mean_content) * 100.0;
            self.mean_c_relative_content = (self.mean_c_content / acum_mean_content) * 100.0;
            self.mean_g_relative_content = (self.mean_g_content / acum_mean_content) * 100.0;
            self.mean_n_relative_content = (self.mean_n_content / acum_mean_content) * 100.0;
            self.mean_gc_relative_content =
                self.mean_g_relative_content + self.mean_c_relative_content;
        }

        if self.detailed_flag {
            self.sum_coverage = 0;
            for i in 0..self.coverage_across_reference.len() {
                let coverage_at_postion = self.coverage_across_reference[i];
                if coverage_at_postion > 0 {
                    // quality
                    self.mapping_quality_across_reference[i] =
                        self.mapping_quality_across_reference[i] / coverage_at_postion as i64;

                    self.sum_coverage_squared +=
                        coverage_at_postion as i64 * coverage_at_postion as i64;
                    self.sum_coverage += coverage_at_postion as i64;
                } else {
                    // make it invalid for histogram
                    self.mapping_quality_across_reference[i] = -1;
                }
            }

            // compute std coverageData
            let mean_coverage = (self.sum_coverage as f64
                / self.coverage_across_reference.len() as f64)
                .floor() as usize;
            self.std_coverage = (self.sum_coverage_squared as f64
                / self.coverage_across_reference.len() as f64
                - mean_coverage as f64 * mean_coverage as f64)
                .sqrt();
        }
    }

    pub fn add_read_alignment_data(&mut self, read_data: &SingleReadData) {
        self.number_of_mapped_bases += read_data.number_of_mapped_bases;
        self.number_of_sequenced_bases += read_data.number_of_sequenced_bases;

        self.number_of_as += read_data.number_of_as;
        self.number_of_cs += read_data.number_of_cs;
        self.number_of_ts += read_data.number_of_ts;
        self.number_of_gs += read_data.number_of_gs;
        self.number_of_ns += read_data.number_of_ns;

        if self.detailed_flag {
            for pos in read_data.coverage_data.iter() {
                self.coverage_across_reference[*pos as usize] += 1;
            }

            for cell in read_data.mapping_quality_data.iter() {
                let pos = cell.get_position();
                let val = cell.get_value();
                self.mapping_quality_across_reference[pos as usize] += val as i64;
                self.acum_mapping_quality += val as f64;
            }
        }
    }

    fn get_window_size(&self) -> i64 {
        self.window_size
    }

    fn get_window_start(&self) -> i64 {
        self.start
    }

    fn get_window_end(&self) -> i64 {
        self.end
    }

    fn get_number_of_mapped_bases(&self) -> i64 {
        self.number_of_mapped_bases
    }

    fn get_number_of_sequenced_bases(&self) -> i64 {
        self.number_of_sequenced_bases
    }

    fn get_number_of_aligned_bases(&self) -> i64 {
        self.number_of_aligned_bases
    }

    fn get_effective_window_length(&self) -> i64 {
        self.effective_window_length
    }

    fn get_num_mapped_bases(&self) -> i64 {
        self.number_of_mapped_bases
    }

    fn get_mean_coverage(&self) -> f64 {
        self.mean_coverage
    }

    fn get_std_coverage(&self) -> f64 {
        self.std_coverage
    }

    fn get_sum_coverage_sequared(&self) -> i64 {
        self.sum_coverage_squared
    }

    fn get_sum_coverage(&self) -> i64 {
        self.sum_coverage
    }

    fn get_mean_mapping_quality(&self) -> f64 {
        self.mean_mapping_quality
    }

    fn get_number_of_as(&self) -> i64 {
        self.number_of_as
    }

    fn get_mean_a_content(&self) -> f64 {
        self.mean_a_content
    }

    fn get_mean_a_relative_content(&self) -> f64 {
        self.mean_a_relative_content
    }

    fn get_number_of_ts(&self) -> i64 {
        self.number_of_ts
    }

    fn get_mean_t_content(&self) -> f64 {
        self.mean_t_content
    }

    fn get_mean_t_relative_content(&self) -> f64 {
        self.mean_t_relative_content
    }

    fn get_number_of_cs(&self) -> i64 {
        self.number_of_cs
    }

    fn get_mean_c_content(&self) -> f64 {
        self.mean_c_content
    }

    fn get_mean_c_relative_content(&self) -> f64 {
        self.mean_c_relative_content
    }

    fn get_number_of_gs(&self) -> i64 {
        self.number_of_gs
    }

    fn get_mean_g_content(&self) -> f64 {
        self.mean_g_content
    }

    fn get_mean_g_relative_content(&self) -> f64 {
        self.mean_g_relative_content
    }

    fn get_number_of_ns(&self) -> i64 {
        self.number_of_ns
    }

    fn get_mean_n_content(&self) -> f64 {
        self.mean_n_content
    }

    fn get_mean_n_relative_content(&self) -> f64 {
        self.mean_n_relative_content
    }

    fn get_mean_gc_content(&self) -> f64 {
        self.mean_gc_content
    }

    fn get_mean_gc_relative_content(&self) -> f64 {
        self.mean_gc_relative_content
    }

    fn get_mean_insert_size(&self) -> f64 {
        self.mean_insert_size
    }

    fn get_coverage_across_reference(&self) -> &Vec<i32> {
        &self.coverage_across_reference
    }

    fn get_mapping_quality_across_reference(&self) -> &Vec<i64> {
        &self.mapping_quality_across_reference
    }

    fn inverse_regions(&mut self) {
        let mut inverse_selected_regions = BitVec::new();
        for x in self.selected_regions.iter() {
            inverse_selected_regions.push(!x);
        }
        self.selected_regions = inverse_selected_regions;
    }

    fn get_region(&self) -> &BitVec {
        &self.selected_regions
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct RegionOverlapLookupTable {}
impl RegionOverlapLookupTable {
    pub fn new() -> Self {
        Self {}
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
struct ReadAlignmentInfo {
    first_read_end_pos: i32,
    second_read_start_pos: i32,
}

impl ReadAlignmentInfo {
    fn new() -> Self {
        ReadAlignmentInfo {
            first_read_end_pos: 0,
            second_read_start_pos: 0,
        }
    }

    fn first_read_end_pos(&self) -> i32 {
        self.first_read_end_pos
    }

    fn second_read_start_pos(&self) -> i32 {
        self.second_read_start_pos
    }
}
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct BamStatsCollector {
    num_mapped_reads: u64,
    num_paired_reads: u64,
    num_supplementary_alignments: u64,
    num_mapped_first_in_pair: u64,
    num_mapped_second_in_pair: u64,
    num_singletons: u64,
    num_marked_duplicates: u64,
    collect_intersecting_read_pairs_flag: bool,
    num_overlapping_read_pairs: u64,
    num_overlapping_bases: u64,
    pairs_collector: HashMap<String, ReadAlignmentInfo>,
    cur_chromosome_id: i32,
}
impl BamStatsCollector {
    pub fn new() -> Self {
        Self {
            num_mapped_reads: 0,
            num_paired_reads: 0,
            num_supplementary_alignments: 0,
            num_mapped_first_in_pair: 0,
            num_mapped_second_in_pair: 0,
            num_singletons: 0,
            num_marked_duplicates: 0,
            collect_intersecting_read_pairs_flag: false,
            num_overlapping_read_pairs: 0,
            num_overlapping_bases: 0,
            pairs_collector: HashMap::new(),
            cur_chromosome_id: -1,
        }
    }
    pub fn get_num_mapped_reads(&self) -> u64 {
        self.num_mapped_reads
    }

    pub fn get_num_supplementary_alignments(&self) -> u64 {
        self.num_supplementary_alignments
    }

    pub fn get_num_mapped_first_in_pair(&self) -> u64 {
        self.num_mapped_first_in_pair
    }

    pub fn get_num_mapped_second_in_pair(&self) -> u64 {
        self.num_mapped_second_in_pair
    }

    pub fn get_num_singletons(&self) -> u64 {
        self.num_singletons
    }

    pub fn get_num_paired_reads(&self) -> u64 {
        self.num_paired_reads
    }

    pub fn get_num_overlapping_read_pairs(&self) -> u64 {
        self.num_overlapping_read_pairs
    }

    pub fn get_num_overlapping_bases(&self) -> u64 {
        self.num_overlapping_bases
    }

    pub fn get_num_marked_duplicates(&self) -> u64 {
        self.num_marked_duplicates
    }

    pub fn enable_intersecting_reads_collection(&mut self) {
        self.collect_intersecting_read_pairs_flag = true;
    }

    pub fn update_stats_and_judge_is_dup(&mut self, record: &Record) -> bool {
        // let flag_value = record.flags(); flagValue & Constants.SAM_FLAG_SUPP_ALIGNMENT
        if record.is_supplementary() {
            self.num_supplementary_alignments += 1;
        } else {
            self.num_mapped_reads += 1;
            if record.is_paired() {
                self.num_paired_reads += 1;
                if record.is_first_in_template() {
                    // don't sure if it means first in pair
                    self.num_mapped_first_in_pair += 1;
                } else if record.is_last_in_template() {
                    self.num_mapped_second_in_pair += 1;
                }
                if record.is_mate_unmapped() {
                    self.num_singletons += 1;
                }
            }
        }

        if record.is_duplicate() {
            self.num_marked_duplicates += 1;
            true
        } else {
            false
        }
    }

    fn report(&self) -> String {
        let mut result = "".to_string();
        result = result + &format!("Num mapped reads: {}\n", self.num_mapped_reads);
        result = result
            + &format!(
                "Num mapped first of pair: {}\n",
                self.num_mapped_first_in_pair
            );
        result = result
            + &format!(
                "Num mapped second of pair: {}\n",
                self.num_mapped_second_in_pair
            );
        result = result + &format!("Num singletons: {}\n", self.num_singletons);
        result
    }

    pub fn finalize_alignment_info(&mut self) {
        if self.cur_chromosome_id != -1 {
            for (_, info) in self.pairs_collector.iter() {
                if info.second_read_start_pos() > 0 {
                    let intersection_size =
                        info.first_read_end_pos() - info.second_read_start_pos() + 1;
                    // calculate intersection
                    if intersection_size > 0 {
                        self.num_overlapping_bases += intersection_size as u64;
                        self.num_overlapping_read_pairs += 1;
                    }
                }
            }
        }
    }

    pub fn collect_paired_read_info(&mut self, read: &Record) {
        if !read.is_paired() || read.is_mate_unmapped() {
            return;
        }

        let chr_id = read.tid();
        if self.cur_chromosome_id != chr_id {
            self.finalize_alignment_info();
            self.cur_chromosome_id = chr_id;
            self.pairs_collector.clear();
        }

        // TODO: what if there are several alignments of the same read?
        let read_name = String::from_utf8(read.qname().to_vec()).unwrap();
        if let Some(info) = self.pairs_collector.get_mut(&read_name) {
            info.second_read_start_pos = read.pos() as i32;
        } else {
            let info = ReadAlignmentInfo {
                first_read_end_pos: (read.cigar().end_pos() - 1) as i32,
                second_read_start_pos: -1,
            };
            self.pairs_collector.insert(read_name, info);
        }
    }
}
#[derive(Debug, Serialize, Deserialize, Clone)]
struct Cell {
    position: i32,
    value: i32,
}

impl Cell {
    fn new(position: i32, value: i32) -> Self {
        Cell { position, value }
    }

    fn get_position(&self) -> i32 {
        self.position
    }

    fn get_value(&self) -> i32 {
        self.value
    }
}
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct SingleReadData {
    number_of_sequenced_bases: i64,
    number_of_mapped_bases: i64,
    number_of_as: i64,
    number_of_ts: i64,
    number_of_cs: i64,
    number_of_gs: i64,
    number_of_ns: i64,
    // These number denotes how many bases are aligned from sequenced bases
    //pub numberOfAlignedBases: i64;
    coverage_data: Vec<i32>,
    mapping_quality_data: Vec<Cell>,

    window_start: i64,
}

impl SingleReadData {
    pub fn new(_window_start: i64) -> Self {
        Self {
            number_of_sequenced_bases: 0,
            number_of_mapped_bases: 0,
            number_of_as: 0,
            number_of_ts: 0,
            number_of_cs: 0,
            number_of_gs: 0,
            number_of_ns: 0,
            coverage_data: vec![],
            mapping_quality_data: vec![],

            window_start: _window_start,
        }
    }

    pub fn get_window_start(&self) -> i64 {
        self.window_start
    }

    pub fn acum_base(&mut self, relative: u64, base: char) {
        self.number_of_sequenced_bases += 1;

        // ATCG content
        match base {
            'A' => self.acum_a(relative),
            'C' => self.acum_c(relative),
            'T' => self.acum_t(relative),
            'G' => self.acum_g(relative),
            _ => {}
        }

        self.coverage_data.push(relative as i32);
    }

    pub fn acum_a(&mut self, _relative: u64) {
        self.number_of_as += 1;
    }

    pub fn acum_c(&mut self, _relative: u64) {
        self.number_of_cs += 1;
    }

    pub fn acum_t(&mut self, _relative: u64) {
        self.number_of_ts += 1;
    }

    pub fn acum_g(&mut self, _relative: u64) {
        self.number_of_gs += 1;
    }

    pub fn acum_num_of_mapped_base(&mut self) {
        self.number_of_mapped_bases += 1;
    }

    pub fn acum_mapping_quality(&mut self, relative: u64, mapping_quality: i32) {
        if mapping_quality != 0 {
            self.mapping_quality_data
                .push(Cell::new(relative as i32, mapping_quality));
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct ReadStatsCollector {
    reads_a_content: Vec<i32>,
    reads_c_content: Vec<i32>,
    reads_g_content: Vec<i32>,
    reads_t_content: Vec<i32>,
    reads_n_content: Vec<i32>,
    reads_clipping_content: Vec<i32>,
    reads_gc_content: Vec<f32>,
    homopolymer_indels: Vec<i32>,
    num_clipped_reads: i32,
    num_read_with_insertion: i32,
    num_read_with_deletion: i32,
    num_insertions: i32,
    num_deletions: i32,
    num_mismatches: i32,
    edit_distance: i32,
    num_bases: i32,
    num_gc: i32,
    min_homopolymer_size: i32,
    prev_base_inside_indel_region_flag: bool,
    homopolymer_starts_inside_indel_region_flag: bool,
    prev_base: u8,
    homopolymer_size: i32,
}

impl ReadStatsCollector {
    const INITIAL_SIZE: usize = 64;

    fn ensure_array_size(array: &mut Vec<i32>, pos: usize) {
        let old_size = array.len();
        if pos >= old_size {
            let mut new_size = 0;
            if old_size * 2 < pos + 1 {
                new_size = pos + 1;
            } else {
                new_size = old_size * 2;
            }
            for _ in 0..new_size - old_size {
                array.push(0);
            }
        }
    }

    fn new(homopolymer_size: i32) -> Self {
        Self {
            reads_a_content: vec![0; Self::INITIAL_SIZE],
            reads_c_content: vec![0; Self::INITIAL_SIZE],
            reads_g_content: vec![0; Self::INITIAL_SIZE],
            reads_t_content: vec![0; Self::INITIAL_SIZE],
            reads_n_content: vec![0; Self::INITIAL_SIZE],
            reads_clipping_content: vec![0; Self::INITIAL_SIZE],
            reads_gc_content: vec![],
            homopolymer_indels: vec![0; 5],
            num_clipped_reads: 0,
            num_read_with_insertion: 0,
            num_read_with_deletion: 0,
            num_insertions: 0,
            num_deletions: 0,
            num_mismatches: 0,
            edit_distance: 0,
            num_bases: 0,
            num_gc: 0,
            min_homopolymer_size: homopolymer_size,
            prev_base_inside_indel_region_flag: false,
            homopolymer_starts_inside_indel_region_flag: false,
            prev_base: 0,
            homopolymer_size: 0,
        }
    }

    fn save_gc(&mut self) {
        if self.num_gc != 0 {
            let gc_content = self.num_gc as f32 / self.num_bases as f32;
            self.reads_gc_content.push(gc_content);
        }

        self.num_bases = 0;
        self.num_gc = 0;
    }

    pub fn collect_base(&mut self, pos: usize, base: u8, inside_indel_region_flag: bool) {
        // ATCG content
        match base as char {
            'A' => {
                self.inc_as_content(pos);
                self.num_bases += 1;
            }
            'C' => {
                self.inc_cs_content(pos);
                self.num_bases += 1;
                self.num_gc += 1;
            }
            'T' => {
                self.inc_ts_content(pos);
                self.num_bases += 1;
            }
            'G' => {
                self.inc_gs_content(pos);
                self.num_bases += 1;
                self.num_gc += 1;
            }
            'N' => {
                self.inc_ns_content(pos);
            }
            _ => {}
        }

        if self.num_bases >= 1000 {
            self.save_gc();
        }
        self.update_homopolymer_stats(base, inside_indel_region_flag);
    }

    fn update_homopolymer_stats(&mut self, base: u8, inside_indel_region: bool) {
        if base == self.prev_base {
            self.homopolymer_size += 1;
        } else {
            if self.prev_base_inside_indel_region_flag
                || self.homopolymer_starts_inside_indel_region_flag
            {
                if self.homopolymer_size >= self.min_homopolymer_size {
                    self.save_homopolymer_data();
                }
            }
            self.homopolymer_size = 1;
            self.homopolymer_starts_inside_indel_region_flag = inside_indel_region;
        }

        self.prev_base = base;
        self.prev_base_inside_indel_region_flag = inside_indel_region;
    }

    fn save_homopolymer_data(&mut self) {
        match self.prev_base as char {
            'A' => self.homopolymer_indels[0] += 1,
            'C' => self.homopolymer_indels[1] += 1,
            'G' => self.homopolymer_indels[2] += 1,
            'T' => self.homopolymer_indels[3] += 1,
            'N' => self.homopolymer_indels[4] += 1,
            _ => (),
        }
    }

    pub fn collect_deleted_base(&mut self, next_base: u8) {
        if next_base != self.prev_base {
            if self.homopolymer_size + 1 >= self.min_homopolymer_size {
                self.save_homopolymer_data();
                self.homopolymer_starts_inside_indel_region_flag = false;
            } else {
                self.homopolymer_starts_inside_indel_region_flag = true;
            }
            self.prev_base = next_base;
            self.homopolymer_size = 1;
        } else {
            self.homopolymer_size += 1;
        }
    }

    fn inc_as_content(&mut self, pos: usize) {
        Self::ensure_array_size(&mut self.reads_a_content, pos);
        self.reads_a_content[pos] += 1;
    }

    fn inc_cs_content(&mut self, pos: usize) {
        Self::ensure_array_size(&mut self.reads_c_content, pos);
        self.reads_c_content[pos] += 1;
    }

    fn inc_ts_content(&mut self, pos: usize) {
        Self::ensure_array_size(&mut self.reads_t_content, pos);
        self.reads_t_content[pos] += 1;
    }

    fn inc_gs_content(&mut self, pos: usize) {
        Self::ensure_array_size(&mut self.reads_g_content, pos);
        self.reads_g_content[pos] += 1;
    }

    fn inc_ns_content(&mut self, pos: usize) {
        Self::ensure_array_size(&mut self.reads_n_content, pos);
        self.reads_n_content[pos] += 1;
    }

    fn inc_clipping_content(&mut self, pos: usize) {
        Self::ensure_array_size(&mut self.reads_clipping_content, pos);
        self.reads_clipping_content[pos] += 1;
    }

    pub fn inc_num_clipped_reads(&mut self) {
        self.num_clipped_reads += 1;
    }

    pub fn inc_num_insertions(&mut self) {
        self.num_insertions += 1;
    }

    pub fn inc_num_deletions(&mut self) {
        self.num_deletions += 1;
    }

    pub fn inc_num_reads_with_insertion(&mut self) {
        self.num_read_with_insertion += 1;
    }

    pub fn inc_num_reads_with_deletion(&mut self) {
        self.num_read_with_deletion += 1;
    }

    pub fn get_num_clipped_reads(&self) -> i32 {
        self.num_clipped_reads
    }

    pub fn get_num_reads_with_deletion(&self) -> i32 {
        self.num_read_with_deletion
    }

    pub fn get_num_reads_with_insertion(&self) -> i32 {
        self.num_read_with_insertion
    }

    pub fn get_homopolymer_indels(&self) -> &Vec<i32> {
        &self.homopolymer_indels
    }

    pub fn get_num_indels(&self) -> i32 {
        self.num_insertions + self.num_deletions
    }

    pub fn get_num_insertions(&self) -> i32 {
        self.num_insertions
    }

    pub fn get_num_deletions(&self) -> i32 {
        self.num_deletions
    }

    pub fn get_reads_as_content(&self) -> &Vec<i32> {
        &self.reads_a_content
    }

    pub fn get_reads_ts_content(&self) -> &Vec<i32> {
        &self.reads_t_content
    }

    pub fn get_reads_cs_content(&self) -> &Vec<i32> {
        &self.reads_c_content
    }

    pub fn get_reads_gs_content(&self) -> &Vec<i32> {
        &self.reads_g_content
    }

    pub fn get_reads_ns_content(&self) -> &Vec<i32> {
        &self.reads_n_content
    }

    pub fn get_reads_clipping_info(&self) -> &Vec<i32> {
        &self.reads_clipping_content
    }

    pub fn get_reads_gc_content(&self) -> &Vec<f32> {
        &self.reads_gc_content
    }

    pub fn inc_num_mismatches(&mut self, mismatches: i32) {
        self.num_mismatches += mismatches;
    }

    pub fn get_num_mismatches(&self) -> i32 {
        self.num_mismatches
    }

    pub fn inc_edit_distance(&mut self, dist: i32) {
        self.edit_distance += dist;
    }

    pub fn get_edit_distance(&self) -> i32 {
        self.edit_distance
    }
}
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct TaskResult {
    reads_data: Option<Vec<SingleReadData>>,
    out_region_reads_data: Option<Vec<SingleReadData>>,
    read_stats_collector: Option<ReadStatsCollector>,
    out_region_read_stats_collector: Option<ReadStatsCollector>,
}
impl TaskResult {
    pub fn new() -> Self {
        Self {
            reads_data: None,
            out_region_reads_data: None,
            read_stats_collector: None,
            out_region_read_stats_collector: None,
        }
    }

    pub fn get_read_stats_collector(&self) -> &ReadStatsCollector {
        self.read_stats_collector.as_ref().unwrap()
    }

    pub fn get_out_region_read_stats_collector(&self) -> &ReadStatsCollector {
        &self.out_region_read_stats_collector.as_ref().unwrap()
    }

    pub fn set_read_stats_collector(&mut self, read_stats_collector: ReadStatsCollector) {
        self.read_stats_collector = Some(read_stats_collector);
    }

    pub fn set_out_region_read_stats_collector(
        &mut self,
        read_stats_collector: ReadStatsCollector,
    ) {
        self.out_region_read_stats_collector = Some(read_stats_collector);
    }

    pub fn set_global_reads_data(&mut self, reads_data: Vec<SingleReadData>) {
        self.reads_data = Some(reads_data);
    }

    pub fn get_read_alignment_data(&self) -> &Vec<SingleReadData> {
        &self.reads_data.as_ref().unwrap()
    }

    pub fn set_out_of_region_reads_data(&mut self, in_region_reads_data: Vec<SingleReadData>) {
        self.out_region_reads_data = Some(in_region_reads_data);
    }

    pub fn get_out_of_region_reads_data(&self) -> &Vec<SingleReadData> {
        self.out_region_reads_data.as_ref().unwrap()
    }
}
pub struct ProcessBunchOfReadsTask<'a> {
    reads: Vec<Record>,
    ctx: &'a mut BamStatsAnalysis,
    // current_window: BamGenomeWindow,
    compute_insert_size_flag: bool,
    is_paired_data_flag: bool,
    analyze_regions_flag: bool,
    compute_outside_stats_flag: bool,
    analysis_results: HashMap<i64, SingleReadData>,
    out_of_regions_results: HashMap<i64, SingleReadData>,
    reads_gc_content: Vec<f32>,
    read_stats_collector: ReadStatsCollector,
    out_of_regions_reads_gc_content: Vec<f32>,
    out_of_regions_read_stats_collector: ReadStatsCollector,
    locator: GenomeLocator,
}

impl<'a> ProcessBunchOfReadsTask<'a> {
    const CIGAR_M: char = 'M';
    const CIGAR_EQ: char = '=';

    pub fn new(
        _reads: Vec<Record>,
        selected_region_flag: bool,
        compute_outside_stats_flag: bool,
        min_homopolymer_size: i32,
        _locator: GenomeLocator,
        _ctx: &'a mut BamStatsAnalysis,
    ) -> Self {
        Self {
            reads: _reads,
            ctx: _ctx,
            compute_insert_size_flag: true,
            is_paired_data_flag: true,
            analyze_regions_flag: selected_region_flag,
            compute_outside_stats_flag: compute_outside_stats_flag,
            analysis_results: HashMap::new(),
            out_of_regions_results: HashMap::new(),
            reads_gc_content: vec![],
            read_stats_collector: ReadStatsCollector::new(min_homopolymer_size),
            out_of_regions_reads_gc_content: vec![],
            out_of_regions_read_stats_collector: ReadStatsCollector::new(min_homopolymer_size),
            locator: _locator,
        }
    }

    pub fn run(&mut self, time_to_task_run: &mut Duration) -> TaskResult {
        let task_run_start = Instant::now();
        let bunch_clone = self.reads.clone();

        for record in &bunch_clone {
            let contig_id = record.tid();
            let contig_option = self.ctx.get_locator().get_contig(contig_id as usize);
            // compute absolute position
            let mut read_absolute_start = -1;
            if let Some(contig) = contig_option {
                read_absolute_start = self
                    .locator
                    .get_absolute_coordinates(contig.name(), record.pos() as i32 + 1);
            }

            let mut alignment: Vec<char> = vec![];
            // compute alignment
            let read_in_region_flag = true;
            if read_in_region_flag {
                // it costs a lot of time, about 0.5~0.6 seconds
                alignment = self.compute_read_alignment_and_collect(record);
            } else if self.analyze_regions_flag && self.compute_outside_stats_flag {
                alignment = self.compute_outside_read_alignment_and_collect(record);
            } else {
                alignment = self.compute_read_alignment(record);
            }

            if alignment.is_empty() {
                continue;
            }

            // insert size
            if self.compute_insert_size_flag && record.is_proper_pair() {
                let bam_stats = self.ctx.get_bam_stats().unwrap();
                let current_window_start = bam_stats.get_current_window_start();
                let current_window = self
                    .ctx
                    .open_windows
                    .get_mut(&current_window_start)
                    .unwrap();
                // update insert size of current_window
                current_window.acum_insert_size(record.insert_size());
            }

            let mapping_quality = record.mapq();
            let read_absolute_end = read_absolute_start + alignment.len() as i64 - 1;

            let current_window_size = self.get_current_window_size();
            let current_window_start = self.get_current_window_start();
            let current_window_end = self.get_current_window_end();
            let current_selected_region = self.get_current_window_selected_region().clone();

            // acum read

            let out_of_bounds_flag = self.process_read_alignment(
                current_window_size,
                current_window_start,
                current_window_end,
                &current_selected_region,
                &alignment,
                read_absolute_start,
                read_absolute_end,
                mapping_quality,
            );

            if out_of_bounds_flag {
                self.propagate_read(
                    &alignment,
                    read_absolute_start,
                    read_absolute_end,
                    mapping_quality,
                );
            }
        }

        self.read_stats_collector.save_gc();

        let mut task_result = TaskResult::new();
        let mut reads_data = vec![];

        for val in self.analysis_results.values() {
            reads_data.push(val.clone());
        }

        task_result.set_global_reads_data(reads_data);
        task_result.set_read_stats_collector(self.read_stats_collector.clone());

        if self.analyze_regions_flag && self.compute_outside_stats_flag {
            self.out_of_regions_read_stats_collector.save_gc();

            let mut out_of_regionreads_data = vec![];
            for val in self.out_of_regions_results.values() {
                out_of_regionreads_data.push(val.clone());
            }

            task_result.set_global_reads_data(out_of_regionreads_data);
            task_result.set_read_stats_collector(self.out_of_regions_read_stats_collector.clone());
        }
        *time_to_task_run += Instant::now() - task_run_start;
        task_result
    }

    fn compute_read_alignment_and_collect(&mut self, read: &Record) -> Vec<char> {
        let alignment_length = read.cigar().end_pos() - read.pos();
        let read_length = read.seq_len();
        if alignment_length <= 0 || read_length == 0 {
            return vec![];
        }
        // precompute total size of alignment
        let mut total_size = 0;
        let cigar = read.cigar();
        let mut read_is_clipped = false;
        let mut read_has_deletion = false;
        let mut read_has_insertion = false;

        for (__, c) in cigar.iter().enumerate() {
            total_size += c.len();
            if c.char() == 'H' || c.char() == 'S' {
                read_is_clipped = true;
            } else if c.char() == 'I' {
                self.read_stats_collector.inc_num_insertions();
                if !read_has_insertion {
                    self.read_stats_collector.inc_num_reads_with_insertion();
                }
                read_has_insertion = true;
            } else if c.char() == 'D' {
                self.read_stats_collector.inc_num_deletions();
                if !read_has_deletion {
                    self.read_stats_collector.inc_num_reads_with_deletion();
                }
                read_has_deletion = true;
            }
        }

        let mut start_pos = 0;
        let mut end_pos = 0;
        // compute extended cigar
        let mut extended_cigar_vector: Vec<char> = vec![' '; total_size as usize];

        for (_, c) in cigar.iter().enumerate() {
            end_pos = start_pos + c.len();
            for pos in start_pos..end_pos {
                extended_cigar_vector[pos as usize] = c.char();
            }
            start_pos = end_pos;
        }

        let mut read_pos = 0;
        let mut alignment_pos = 0;
        let mut alignment_vector: Vec<char> = vec![0 as char; alignment_length as usize];
        let read_bases = read.seq();

        // function collect base takes 0.15~0.17 seconds, this cycle costs about 0.4 seconds
        for pos in 0..extended_cigar_vector.len() {
            let cigar_char = extended_cigar_vector[pos];

            if cigar_char == 'M' || cigar_char == '=' {
                let base = read_bases[read_pos];
                self.read_stats_collector
                    .collect_base(read_pos, base, false);
                read_pos += 1;
                alignment_vector[alignment_pos as usize] = base as char;
                alignment_pos += 1;
            } else if cigar_char == 'I' {
                let base = read_bases[read_pos];
                self.read_stats_collector.collect_base(read_pos, base, true);
                read_pos += 1;
            } else if cigar_char == 'D' {
                let next_cigar_pos = pos + 1;
                if next_cigar_pos < extended_cigar_vector.len()
                    && extended_cigar_vector[next_cigar_pos] != 'D'
                {
                    let mut next_base: u8 = 0; // origin is -1
                    if read_pos + 1 < read_bases.len() {
                        next_base = read_bases[read_pos + 1];
                    }
                    self.read_stats_collector.collect_deleted_base(next_base);
                }
                alignment_vector[alignment_pos as usize] = '-';
                alignment_pos += 1;
            } else if cigar_char == 'N' {
                alignment_vector[alignment_pos as usize] = 'N';
                alignment_pos += 1;
            } else if cigar_char == 'S' {
                self.read_stats_collector.inc_clipping_content(read_pos);
                read_pos += 1;
            } else if cigar_char == 'H' {
                self.read_stats_collector.inc_clipping_content(read_pos);
            } else if cigar_char == 'P' {
                alignment_vector[alignment_pos as usize] = '-';
                alignment_pos += 1;
            }
        }

        if read_is_clipped {
            self.read_stats_collector.inc_num_clipped_reads();
        }

        let num_mismatches = self.compute_read_mismatches(read);

        self.read_stats_collector
            .inc_num_mismatches(num_mismatches as i32);

        match read.aux(b"NM") {
            Ok(value) => {
                if let Aux::U8(v) = value {
                    self.read_stats_collector.inc_edit_distance(v as i32);
                }
            }
            Err(e) => {}
        }

        alignment_vector
    }

    fn compute_outside_read_alignment_and_collect(&mut self, read: &Record) -> Vec<char> {
        let alignment_length = read.cigar().end_pos() - read.pos();
        let read_length = read.seq_len();
        if alignment_length <= 0 || read_length == 0 {
            return vec![];
        }
        // precompute total size of alignment
        let mut total_size = 0;
        let cigar = read.cigar();
        let mut read_is_clipped = false;
        let mut read_has_deletion = false;
        let mut read_has_insertion = false;
        for (_i, c) in cigar.iter().enumerate() {
            total_size += c.len();
            if c.char() == 'H' || c.char() == 'S' {
                read_is_clipped = true;
            } else if c.char() == 'I' {
                self.out_of_regions_read_stats_collector
                    .inc_num_insertions();
                if !read_has_insertion {
                    self.out_of_regions_read_stats_collector
                        .inc_num_reads_with_insertion();
                }
                read_has_insertion = true;
            } else if c.char() == 'D' {
                self.out_of_regions_read_stats_collector.inc_num_deletions();
                if !read_has_deletion {
                    self.out_of_regions_read_stats_collector
                        .inc_num_reads_with_deletion();
                }
                read_has_deletion = true;
            }
        }

        let mut start_pos = 0;
        let mut end_pos = 0;
        // compute extended cigar
        let mut extended_cigar_vector: Vec<char> = vec![' '; total_size as usize];
        for (_, c) in cigar.iter().enumerate() {
            end_pos = start_pos + c.len();
            for pos in start_pos..end_pos {
                extended_cigar_vector[pos as usize] = c.char();
            }
            start_pos = end_pos;
        }

        let mut alignment_vector: Vec<char> = vec![0 as char; alignment_length as usize];
        let mut read_pos = 0;
        let mut alignment_pos = 0;
        let read_bases = read.seq();

        for pos in 0..extended_cigar_vector.len() {
            let cigar_char = extended_cigar_vector[pos];
            let base = read_bases[read_pos];
            if cigar_char == 'M' || cigar_char == '=' {
                self.out_of_regions_read_stats_collector
                    .collect_base(read_pos, base, false);
                read_pos += 1;
                alignment_vector[alignment_pos as usize] = base as char;
                alignment_pos += 1;
            } else if cigar_char == 'I' {
                self.out_of_regions_read_stats_collector
                    .collect_base(read_pos, base, true);
                read_pos += 1;
            } else if cigar_char == 'D' {
                let next_cigar_pos = pos + 1;
                if next_cigar_pos < extended_cigar_vector.len()
                    && extended_cigar_vector[next_cigar_pos] != 'D'
                {
                    let mut next_base: u8 = 0; // origin is -1
                    if read_pos + 1 < read_bases.len() {
                        next_base = read_bases[read_pos + 1];
                    }
                    self.out_of_regions_read_stats_collector
                        .collect_deleted_base(next_base);
                }
                alignment_vector[alignment_pos as usize] = '-';
                alignment_pos += 1;
            } else if cigar_char == 'N' {
                alignment_vector[alignment_pos as usize] = 'N';
                alignment_pos += 1;
            } else if cigar_char == 'S' {
                self.out_of_regions_read_stats_collector
                    .inc_clipping_content(read_pos);
                read_pos += 1;
            } else if cigar_char == 'H' {
                self.out_of_regions_read_stats_collector
                    .inc_clipping_content(read_pos);
            } else if cigar_char == 'P' {
                alignment_vector[alignment_pos as usize] = '-';
                alignment_pos += 1;
            }
        }

        if read_is_clipped {
            self.out_of_regions_read_stats_collector
                .inc_num_clipped_reads();
        }

        let num_mismatches = self.compute_read_mismatches(read);
        self.out_of_regions_read_stats_collector
            .inc_num_mismatches(num_mismatches as i32);

        match read.aux(b"MD") {
            Ok(value) => {
                if let Aux::U8(v) = value {
                    self.out_of_regions_read_stats_collector
                        .inc_edit_distance(v as i32);
                }
            }
            Err(e) => {}
        }

        alignment_vector
    }

    fn compute_read_mismatches(&mut self, read: &Record) -> usize {
        let mut num_mismatches = 0;
        let mut md_attr_vec: Vec<u8> = vec![];
        match read.aux(b"MD") {
            Ok(value) => {
                if let Aux::String(v) = value {
                    md_attr_vec = v.as_bytes().to_vec();
                }
            }
            Err(e) => return num_mismatches,
        }

        for i in 0..md_attr_vec.len() {
            let c = md_attr_vec[i] as char;
            if c == 'A' || c == 'T' || c == 'C' || c == 'G' {
                num_mismatches += 1;
            }
        }

        num_mismatches
    }

    fn compute_read_alignment(&self, read: &Record) -> Vec<char> {
        let alignment_length = read.cigar().end_pos() - read.pos();
        let read_length = read.seq_len();
        if alignment_length <= 0 || read_length == 0 {
            return vec![];
        }
        // precompute total size of alignment
        let mut total_size = 0;
        let cigar = read.cigar();
        let mut read_is_clipped = false;
        let mut read_has_deletion = false;
        let mut read_has_insertion = false;
        for (_i, c) in cigar.iter().enumerate() {
            total_size += c.len();
        }

        // compute extended cigar
        let mut extended_cigar_vector: Vec<char> = vec![];
        for (_, c) in cigar.iter().enumerate() {
            for j in 0..c.len() as usize {
                extended_cigar_vector.push(c.char());
            }
        }

        let mut alignment_vector: Vec<char> = vec![];
        let mut read_pos = 0;
        let read_bases = read.seq();

        for pos in 0..extended_cigar_vector.len() {
            let cigar_char = extended_cigar_vector[pos];
            let base = read_bases[read_pos];
            if cigar_char == 'M' || cigar_char == '=' {
                read_pos += 1;
                alignment_vector.push(base as char);
            } else if cigar_char == 'I' {
                read_pos += 1;
            } else if cigar_char == 'D' {
                alignment_vector.push('-');
            } else if cigar_char == 'N' {
                alignment_vector.push('N');
            } else if cigar_char == 'S' {
                read_pos += 1;
            } else if cigar_char == 'H' {
            } else if cigar_char == 'P' {
                alignment_vector.push('-');
            }
        }

        alignment_vector
    }

    fn process_read_alignment(
        &mut self,
        window_size: i64,
        window_start: i64,
        window_end: i64,
        selected_region: &BitVec,
        alignment: &Vec<char>,
        read_start: i64,
        read_end: i64,
        mapping_quality: u8,
    ) -> bool {
        // working variables
        let mut out_of_bounds_flag = false;

        let mut read_data = Self::get_window_data(window_start, &mut self.analysis_results);

        if read_start > read_end {
            self.ctx.inc_number_of_reads_with_start_greater_then_end();
        }
        if read_end > window_end {
            out_of_bounds_flag = true;
        }

        let mut relative = -1;
        let mut pos = -1;
        // run Read
        for j in read_start..read_end + 1 {
            relative = j - window_start;
            pos = j - read_start;
            if relative < 0 {
                continue;
            } else if relative >= window_size {
                break;
            } else if self.analyze_regions_flag {
                let inside_of_region_flag = selected_region.get(relative as usize).unwrap();

                if inside_of_region_flag {
                    if self.compute_outside_stats_flag {
                        read_data = Self::get_window_data(window_start, &mut self.analysis_results);
                    }
                } else if self.compute_outside_stats_flag {
                    read_data =
                        Self::get_window_data(window_start, &mut self.out_of_regions_results);
                } else {
                    continue;
                }
            }

            let nucleotide = alignment[pos as usize];

            // aligned bases
            read_data.number_of_mapped_bases += 1;
            if nucleotide != '-' && nucleotide != 'N' {
                // mapping quality
                read_data.acum_mapping_quality(relative as u64, mapping_quality as i32);
                // base stats
                read_data.acum_base(relative as u64, nucleotide);
            } else if nucleotide == 'N' {
                read_data.number_of_ns += 1;
            }
        }

        out_of_bounds_flag
    }

    fn propagate_read(
        &mut self,
        alignment: &Vec<char>,
        read_start: i64,
        read_end: i64,
        mapping_quality: u8,
    ) {
        let num_of_processed_windows = self.get_number_of_processed_windows();
        let num_of_windows = self.get_total_number_of_windows();

        let mut index = num_of_processed_windows + 1;
        let mut out_of_bounds = true;
        let mut adjacent_window: &BamGenomeWindow;
        while out_of_bounds && index < num_of_windows {
            let window_start = self
                .ctx
                .bam_stats
                .as_ref()
                .unwrap()
                .get_window_start(index as usize);
            // synchronized (lock) {
            adjacent_window = self.ctx.get_open_window(window_start);
            let window_size = adjacent_window.get_window_size();
            let window_start = adjacent_window.get_window_start();
            let window_end = adjacent_window.get_window_end();
            let selected_region = adjacent_window.get_region().clone();

            if self.compute_outside_stats_flag {
                self.ctx.get_outside_open_window(window_start);
            }
            // }

            // acum read
            out_of_bounds = self.process_read_alignment(
                window_size,
                window_start,
                window_end,
                &selected_region,
                alignment,
                read_start,
                read_end,
                mapping_quality,
            );
            index += 1;
        }
    }

    fn get_window_data(
        window_start: i64,
        results_map: &mut HashMap<i64, SingleReadData>,
    ) -> &mut SingleReadData {
        if results_map.contains_key(&window_start) {
            results_map.get_mut(&window_start).unwrap()
        } else {
            let data = SingleReadData::new(window_start);
            results_map.insert(window_start, data);
            results_map.get_mut(&window_start).unwrap()
        }
    }

    fn get_current_window_size(&self) -> i64 {
        self.ctx.get_current_window().unwrap().get_window_size()
    }

    fn get_current_window_start(&self) -> i64 {
        self.ctx.get_current_window().unwrap().get_window_start()
    }

    fn get_current_window_end(&self) -> i64 {
        self.ctx.get_current_window().unwrap().get_window_end()
    }

    fn get_current_window_selected_region(&self) -> &BitVec {
        self.ctx.get_current_window().unwrap().get_region()
    }

    fn get_number_of_processed_windows(&self) -> i32 {
        self.ctx
            .bam_stats
            .as_ref()
            .unwrap()
            .get_number_of_processed_windows()
    }

    fn get_total_number_of_windows(&self) -> i32 {
        self.ctx.bam_stats.as_ref().unwrap().get_number_of_windows()
    }

    fn get_read_stats_collector(&mut self, read: &Record) -> Option<&mut ReadStatsCollector> {
        Some(&mut self.read_stats_collector)
    }

    fn compute_num_mismatches(read: &Record) -> i32 {
        let mut num_mismatches = 0;

        // if let Some(md_attr) = read.aux(b"MD").and_then(|a| a.string()) {
        //     for c in md_attr.chars() {
        //         if c == 'A' || c == 'C' || c == 'G' || c == 'T' {
        //             num_mismatches += 1;
        //         }
        //     }
        // }

        num_mismatches
    }
}

struct FinalizeWindowTask<'a> {
    bam_stats: &'a mut BamStats,
    window: &'a mut BamGenomeWindow,
}

impl<'a> FinalizeWindowTask<'a> {
    pub fn new(_bam_stats: &'a mut BamStats, _window: &'a mut BamGenomeWindow) -> Self {
        Self {
            bam_stats: _bam_stats,
            window: _window,
        }
    }

    pub fn run(&mut self) -> i32 {
        self.window.compute_descriptors();
        self.bam_stats.add_window_information(self.window);

        // todo!() error handling
        return 0;
    }
}

#[derive(Debug, Clone)]
pub struct BamStatsAnalysis {
    bam_file: String,

    // reference
    reference_file: String,
    reference_available_flag: bool,
    reference: Vec<u8>,
    reference_size: i64,
    number_of_reference_config: usize,

    // currentWindow management
    number_of_windows: i32,
    effective_number_of_window: i32,
    /// the number of bases in a window
    window_size: i32,

    // coordinates transformer
    locator: GenomeLocator,

    // globals
    number_of_reads: usize,
    number_of_valid_reads: usize,
    number_of_secondary_alignments: usize,
    number_of_duplicates_skip: usize,
    number_of_correct_strand_reads: usize,
    number_of_problematic_reads: usize,
    number_of_reads_with_start_greater_than_end: usize,

    // statistics
    bam_stats: Option<BamStats>,
    // private Logger logger;
    // private LoggerThread loggerThread;

    // working variables
    current_window: Option<BamGenomeWindow>,
    open_windows: HashMap<i64, BamGenomeWindow>,
    thread_number: usize,
    num_reads_in_bunch: i32,
    progress: usize,
    min_homopolymer_size: i32,

    // nucleotide reporting
    out_dir: String,

    // gff support
    selected_regions_available_flag: bool,
    feature_file: String,
    num_of_selected_regions: usize,

    // inside region analysis
    inside_reference_sizes: usize,
    skip_marked_duplicates_flag: bool,
    skip_detected_duplicates_flag: bool,
    collect_intersecting_paired_end_reads_flag: bool,

    // outside region analysis
    compute_outside_stats_flag: bool,
    current_outside_window: Option<BamGenomeWindow>,
    open_outside_windows: HashMap<i64, BamGenomeWindow>,
    outside_bam_stats: Option<BamStats>,

    // read size
    acum_read_size: usize,
    max_read_size: usize,
    min_read_size: usize,

    // region
    selected_region_starts: Vec<usize>,
    selected_region_ends: Vec<usize>,
    region_overlap_lookup_table: RegionOverlapLookupTable,
    protocol: LibraryProtocol,

    bam_stats_collector: BamStatsCollector,
    outside_bam_stats_collector: BamStatsCollector,

    // chromosome
    chromosome_window_indexes: Vec<usize>,

    max_size_of_task_queue: usize,

    // reporting
    active_reporting_flag: bool,
    save_coverage_flag: bool,
    non_zero_coverage_only_flag: bool,
    path_to_coverage_report: String,
    window_block_size_to_report: usize,

    // analysis
    is_paired_data_flag: bool,
    results: Vec<TaskResult>,
    // Future<Integer> finalizeWindowResult;
    // #[serde(skip_serializing)]
    // #[serde(skip_deserializing)]
    finalize_window_result: i32,
    time_to_calc_overlappers: usize,
    pg_program: String,
    pg_command_string: String,
    reads_bunch: Vec<Record>,
    // // multi thread
    // worker_thread_pool: Option<ExecutorService<_,_>>,
}

impl BamStatsAnalysis {
    // lazy_static! {
    //     static ref GENOME_GC_CONTENT_MAP: HashMap<&'static str, &'static str> = {
    //         let mut map = HashMap::new();
    //         map.insert("HUMAN (hg19)", "species/human.hg19.gc_histogram.txt");
    //         map.insert("MOUSE (mm9)", "species/mouse.mm9.gc_histogram.txt");
    //         map.insert("MOUSE (mm10)", "species/mouse.mm10.gc_histogram.txt");
    //         map
    //     };
    // }

    pub const WARNING_ID_CHROMOSOME_NOT_FOUND: &str = "Some regions are not loaded";
    pub const WARNING_ID_NO_MAPPED_READS: &str = "NO MAPPED READS";
    pub const WARNING_NO_MARKED_DUPLICATES: &str = "No flagged duplicates are detected";
    pub const HUMAN_GENOME_ID: &str = "HUMAN (hg19)";
    pub const MOUSE_GENOME_ID: &str = "MOUSE (mm9)";
    pub const MM10_GENOME_ID: &str = "MOUSE (mm10)";
    pub const HUMAN_GENOME_NAME: &str = "human";
    pub const MOUSE_GENOME_NAME: &str = "mouse";

    pub fn new(_bam_file: String) -> Self {
        Self {
            bam_file: _bam_file,
            // reference
            reference_file: "".to_string(),
            reference_available_flag: false,
            reference: vec![],
            reference_size: 0,
            number_of_reference_config: 0,

            // currentWindow management
            number_of_windows: Constants::DEFAULT_NUMBER_OF_WINDOWS,
            effective_number_of_window: 0,
            window_size: 0,

            // coordinates transformer
            locator: GenomeLocator::new(),

            // globals
            number_of_reads: 0,
            number_of_valid_reads: 0,
            number_of_secondary_alignments: 0,
            number_of_duplicates_skip: 0,
            number_of_correct_strand_reads: 0,
            number_of_problematic_reads: 0,
            number_of_reads_with_start_greater_than_end: 0,

            // statistics
            bam_stats: None,
            // private Logger logger;
            // private LoggerThread loggerThread;

            // working variables
            current_window: None,
            open_windows: HashMap::new(),
            thread_number: 4,
            num_reads_in_bunch: Constants::DEFAULT_CHUNK_SIZE,
            progress: 0,
            min_homopolymer_size: Constants::DEFAULT_HOMOPOLYMER_SIZE,

            // nucleotide reporting
            out_dir: ".".to_string(),

            // gff support
            selected_regions_available_flag: false,
            feature_file: "".to_string(),
            num_of_selected_regions: 0,

            // inside region analysis
            inside_reference_sizes: 0,
            skip_marked_duplicates_flag: false,
            skip_detected_duplicates_flag: false,
            collect_intersecting_paired_end_reads_flag: false,

            // outside region analysis
            compute_outside_stats_flag: false,
            current_outside_window: None,
            open_outside_windows: HashMap::new(),
            outside_bam_stats: None,

            // read size
            acum_read_size: 0,
            max_read_size: 0,
            min_read_size: i32::MAX as usize,

            // region
            selected_region_starts: vec![],
            selected_region_ends: vec![],
            region_overlap_lookup_table: RegionOverlapLookupTable::new(),
            protocol: LibraryProtocol::NonStrandSpecific,

            bam_stats_collector: BamStatsCollector::new(),
            outside_bam_stats_collector: BamStatsCollector::new(),

            // chromosome
            chromosome_window_indexes: vec![],

            max_size_of_task_queue: 10, // don't know if it works

            // reporting
            active_reporting_flag: false,
            save_coverage_flag: false,
            non_zero_coverage_only_flag: true,
            path_to_coverage_report: "".to_string(),
            window_block_size_to_report: 50,

            // analysis
            is_paired_data_flag: false,
            results: vec![],
            // Future<Integer> finalizeWindowResult;
            finalize_window_result: -1,
            time_to_calc_overlappers: 0,
            pg_program: "".to_string(),
            pg_command_string: "".to_string(),
            reads_bunch: vec![],
        }

        // Set the region file
        // reporting
        // setup logging
    }

    pub fn to_hashmap(&self, header: &Header) -> HashMap<String, Vec<LinearMap<String, String>>> {
        let mut header_map = HashMap::default();

        lazy_static! {
            static ref REC_TYPE_RE: Regex = Regex::new(r"@([A-Z][A-Z])").unwrap();
            static ref TAG_RE: Regex = Regex::new(r"([A-Za-z][A-Za-z0-9]):([ -~]+)").unwrap();
        }

        let header_string = String::from_utf8(header.to_bytes()).unwrap();

        for line in header_string.split('\n').filter(|x| !x.is_empty()) {
            let parts: Vec<_> = line.split('\t').filter(|x| !x.is_empty()).collect();
            // assert!(rec_type_re.is_match(parts[0]));
            let record_type = REC_TYPE_RE
                .captures(parts[0])
                .unwrap()
                .get(1)
                .unwrap()
                .as_str()
                .to_owned();
            if record_type.eq("CO") {
                continue;
            }
            let mut field = LinearMap::default();
            for part in parts.iter().skip(1) {
                let cap_option = TAG_RE.captures(part);
                if cap_option.is_none() {
                    continue;
                }
                let cap = cap_option.unwrap();
                let tag = cap.get(1).unwrap().as_str().to_owned();
                let value = cap.get(2).unwrap().as_str().to_owned();
                field.insert(tag, value);
            }
            header_map
                .entry(record_type)
                .or_insert_with(Vec::new)
                .push(field);
        }
        header_map
    }

    fn load_and_init(&mut self) {
        let mut last_action_done = "Loading sam header...";
        println!("{}", last_action_done);
        let bam_reader = Reader::from_path(&self.bam_file).unwrap();

        let header = Header::from_template(bam_reader.header());

        let hash_header = self.to_hashmap(&header);
        if !hash_header.is_empty() {
            if hash_header.contains_key("HD") {
                let hd_records = &hash_header["HD"][0];
                if let Some(val) = hd_records.get("SO") {
                    if val != "coordinate" {
                        println!("According to header the BAM file is not sorted by coordinate!");
                    }
                } else {
                    println!("Non-standard header SortOrder value!");
                }
            } else {
                println!("@HD line is not presented in the BAM file header.");
            }
        }

        // load locator
        last_action_done = "Loading locator...";
        println!("{}", last_action_done);
        self.load_locator(&header);
        self.load_program_records(&header);

        // load reference
        last_action_done = "Loading reference...";
        println!("{}", last_action_done);
        self.load_reference();

        // init window set
        self.window_size = self.compute_window_size(self.reference_size, self.number_of_windows);
        let window_positions: Vec<i64> = self.compute_window_positions(self.window_size);
        self.effective_number_of_window = window_positions.len() as i32;

        if self.effective_number_of_window > Constants::DEFAULT_STABLIZED_WINDOW_PROPORTION {
            self.window_block_size_to_report = self.effective_number_of_window as usize / 10;
        }

        if self.skip_detected_duplicates_flag || self.skip_marked_duplicates_flag {
            let mut skip_duplicates_report = "";
            if self.skip_detected_duplicates_flag && self.skip_marked_duplicates_flag {
                skip_duplicates_report = "Both flagged and estimated by Qualimap algorithm duplicate alignments will be skipped...";
            } else if self.skip_marked_duplicates_flag {
                skip_duplicates_report = "Only flagged duplicate alignments will be skipped...";
            } else {
                skip_duplicates_report =
                    "Only duplicate alignments estimated by Qualimap algorithm will be skipped...";
            }
            println!("{}", skip_duplicates_report);
        }

        println!(
            "Number of windows: {}, effective number of windows: {}",
            self.number_of_windows, self.effective_number_of_window
        );
        println!("Chunk of reads size: {}", self.num_reads_in_bunch);
        println!("Number of threads: {}", self.thread_number);
        // initialize BamStats
        let mut bam_stats_instance = BamStats::new(
            "genome".to_string(),
            self.locator.clone(),
            self.reference_size,
            self.effective_number_of_window as usize,
        );

        bam_stats_instance.set_source_file(self.bam_file.clone());
        bam_stats_instance.set_window_references("w", window_positions);
        self.bam_stats = Some(bam_stats_instance);

        if self.collect_intersecting_paired_end_reads_flag {
            self.bam_stats_collector
                .enable_intersecting_reads_collection();
        }

        if self.save_coverage_flag {
            // todo!()
            // bamStats.activateCoverageReporting(pathToCoverageReport, nonZeroCoverageOnly);
        }

        // regions
        if self.selected_regions_available_flag {
            // todo!()
            // load selected regions

            // outside of regions stats
            if self.compute_outside_stats_flag {}
        }

        self.current_window = self.next_window(true);
    }
    pub fn run(&mut self) {
        let start_time = Instant::now();
        // let mut worker_thread_pool = Executors::new_fixed_thread_pool(self.thread_number as u32)
        //     .expect("Failed to create the thread pool");

        self.load_and_init();
        println!("Time to init: {:?}", Instant::now() - start_time);

        // run reads
        // todo!()        results = new ArrayList<Future<ProcessBunchOfReadsTask.Result>>(); timeToCalcOverlappers = 0;
        let mut bam_reader = Reader::from_path(&self.bam_file).unwrap();
        let mut record = Record::new();
        let mut process_time = Duration::seconds(0);
        let mut time_to_finalize_and_get_next_window = Duration::seconds(0);
        let mut time_to_analyze_reads_bunch = Duration::seconds(0);
        let mut time_to_finish_reads_bunch = Duration::seconds(0);
        let mut time_to_task_run = Duration::seconds(0);
        while let Some(result) = bam_reader.read(&mut record) {
            match result {
                Ok(_) => {
                    let start_time = Instant::now();
                    self.process_sequence(
                        &record,
                        &mut time_to_finalize_and_get_next_window,
                        &mut time_to_analyze_reads_bunch,
                        &mut time_to_finish_reads_bunch,
                        &mut time_to_task_run,
                    );
                    process_time += Instant::now() - start_time;
                }
                Err(_) => {
                    self.number_of_problematic_reads += 1;
                }
            }
        }

        let bam_stats_old = self.bam_stats.clone().unwrap();
        let last_process_start = Instant::now();
        if !self.reads_bunch.is_empty()
            || bam_stats_old.get_number_of_processed_windows()
                < bam_stats_old.get_number_of_windows()
        {
            let num_windows = bam_stats_old.get_number_of_windows();
            let last_position = bam_stats_old.get_window_end(num_windows as usize - 1) + 1;
            self.collect_analysis_results(
                &mut time_to_analyze_reads_bunch,
                &mut time_to_finish_reads_bunch,
                &mut time_to_task_run,
            );

            //finalize
            self.finalize_and_get_next_window(
                last_position,
                true,
                &mut time_to_finalize_and_get_next_window,
            );

            if self.selected_regions_available_flag && self.compute_outside_stats_flag {
                self.current_outside_window
                    .as_mut()
                    .unwrap()
                    .inverse_regions();
                self.finalize_and_get_next_outside_window(last_position, true);
            }
        }
        println!("Finish All");

        process_time += Instant::now() - last_process_start;
        println!("Time to process sequence: {:?}", process_time);
        println!(
            "Time to finalize current and get next: {:?}",
            time_to_finalize_and_get_next_window
        );
        println!("Time to analyze bunch: {:?}", time_to_analyze_reads_bunch);
        println!("Time to finish bunch: {:?}", time_to_finish_reads_bunch);
        println!("Time to task run: {:?}", time_to_task_run);

        // todo!() delete worker thread pool
        // workerThreadPool.shutdown();
        // workerThreadPool.awaitTermination(2, TimeUnit.MINUTES);

        let end_time = Instant::now();

        let bam_stats_final = self.bam_stats.as_mut().unwrap();
        println!(
            "Total processed windows:{}",
            bam_stats_final.get_number_of_processed_windows()
        );
        println!("Number of reads: {}", self.number_of_reads);
        println!("Number of valid reads: {}", self.number_of_valid_reads);
        println!(
            "Number of correct strand reads:{}",
            self.number_of_correct_strand_reads
        );

        if self.number_of_reads_with_start_greater_than_end > 0 {
            println!(
                "{}  read alignments have start greater than end",
                self.number_of_reads_with_start_greater_than_end
            );
        }

        if self.number_of_problematic_reads > 0 {
            println!(
                "SAMRecordParser marked {} problematic reads.",
                self.number_of_problematic_reads
            );
        }

        if self.collect_intersecting_paired_end_reads_flag {
            self.bam_stats_collector.finalize_alignment_info();
        }

        println!("\nInside of regions...");
        println!("{}", self.bam_stats_collector.report());

        if self.compute_outside_stats_flag {
            println!("\nOutside of regions...");
            println!("{}", self.outside_bam_stats_collector.report());
        }

        println!("Time taken to analyze reads: {:?}", end_time - start_time);

        if self.number_of_reads == 0 {
            panic!("The BAM file is empty or corrupt")
        }

        bam_stats_final.number_of_reads = self.number_of_reads as i64;
        bam_stats_final.num_duplicates_skipped = self.number_of_duplicates_skip as u64;
        bam_stats_final.set_skip_duplicates_mode(
            self.skip_marked_duplicates_flag,
            self.skip_detected_duplicates_flag,
        );
        bam_stats_final.number_of_secondary_alignments = self.number_of_secondary_alignments as i64;
        if self.skip_detected_duplicates_flag && self.number_of_duplicates_skip == 0 {
            // Todo!();
            // bamStats.addWarning(WARNING_NO_MARKED_DUPLICATES,
            //     "Make sure duplicate alignments are flagged in the BAM file or apply a different skip duplicates mode.");
        }

        let mut total_number_of_mapped_reads = self.bam_stats_collector.get_num_mapped_reads();
        let mut total_number_of_paired_reads = self.bam_stats_collector.get_num_paired_reads();
        let mut total_number_of_mapped_first_of_pair =
            self.bam_stats_collector.get_num_mapped_first_in_pair();
        let mut total_number_of_mapped_second_of_pair =
            self.bam_stats_collector.get_num_mapped_second_in_pair();
        let mut total_number_of_singletons = self.bam_stats_collector.get_num_singletons();
        let mut total_number_of_supp_alignments =
            self.bam_stats_collector.get_num_supplementary_alignments();
        let mapped_reads_in_region = total_number_of_mapped_reads > 0;

        if self.selected_regions_available_flag {
            bam_stats_final.num_selected_regions = self.num_of_selected_regions as i32;
            if self.compute_outside_stats_flag {
                // Size was calculated two times during the analysis
                self.inside_reference_sizes = self.inside_reference_sizes / 2;
            }
            bam_stats_final.in_region_reference_size = self.inside_reference_sizes as i64;

            // set bam_stats with metrics inside of region
            bam_stats_final.number_of_mapped_reads_in_regions = total_number_of_mapped_reads as i64;
            bam_stats_final.number_of_paired_reads_in_regions = total_number_of_paired_reads as i64;
            bam_stats_final.number_of_mapped_first_of_pair_in_regions =
                total_number_of_mapped_first_of_pair as i64;
            bam_stats_final.number_of_mapped_second_of_pair_in_regions =
                total_number_of_mapped_second_of_pair as i64;
            bam_stats_final.number_of_singletons_in_regions = total_number_of_singletons as i64;
            bam_stats_final.num_correct_strand_reads = self.number_of_correct_strand_reads as i64;

            // update totals from outside collectors
            total_number_of_mapped_reads = total_number_of_mapped_reads
                + self.outside_bam_stats_collector.get_num_mapped_reads();
            total_number_of_paired_reads =
                total_number_of_paired_reads + self.bam_stats_collector.get_num_paired_reads();
            total_number_of_mapped_first_of_pair = total_number_of_mapped_first_of_pair
                + self
                    .outside_bam_stats_collector
                    .get_num_mapped_first_in_pair();
            total_number_of_mapped_second_of_pair = total_number_of_mapped_second_of_pair
                + self
                    .outside_bam_stats_collector
                    .get_num_mapped_second_in_pair();
            total_number_of_singletons =
                total_number_of_singletons + self.outside_bam_stats_collector.get_num_singletons();
        }

        // set bam_stats with total information
        bam_stats_final.number_of_mapped_reads = total_number_of_mapped_reads as i64;
        bam_stats_final.number_of_paired_reads = total_number_of_paired_reads as i64;
        bam_stats_final.number_of_mapped_first_of_pair =
            total_number_of_mapped_first_of_pair as i64;
        bam_stats_final.number_of_mapped_second_of_pair =
            total_number_of_mapped_second_of_pair as i64;
        bam_stats_final.number_of_singletons = total_number_of_singletons as i64;
        bam_stats_final.number_of_supp_alignments = total_number_of_supp_alignments as i64;
        if self.collect_intersecting_paired_end_reads_flag {
            bam_stats_final.set_number_of_intersecting_read_pairs(
                self.bam_stats_collector.get_num_overlapping_read_pairs(),
                self.bam_stats_collector.get_num_overlapping_bases(),
            )
        }

        bam_stats_final.reference_size = self.reference_size;
        bam_stats_final.number_of_reference_contigs = self.locator.get_contigs().len() as i64;
        bam_stats_final.read_max_size = self.max_read_size as i32;
        bam_stats_final.read_min_size = self.min_read_size as i32;
        bam_stats_final.read_mean_size = self.acum_read_size as f64 / self.number_of_reads as f64;
        bam_stats_final.num_detected_duplicate_reads =
            self.bam_stats_collector.num_marked_duplicates;

        self.is_paired_data_flag = bam_stats_final.number_of_paired_reads > 0
            && (bam_stats_final.number_of_mapped_reads > bam_stats_final.number_of_singletons);

        if mapped_reads_in_region {
            println!("Computing descriptors...");
            bam_stats_final.compute_descriptors();
            println!("Computing per chromosome statistics...");
            bam_stats_final
                .compute_chromosome_stats(&self.locator, &self.chromosome_window_indexes);
            println!("Computing histograms...");
            bam_stats_final.compute_histograms();
        } else {
            println!("\nWARNING: number of mapped reads equals zero");
            // todo()!
            // bamStats.addWarning(WARNING_ID_NO_MAPPED_READS, "Total number of mapped reads or mapped reads in region equals zero.\n" +
            // "For more details, check the number of Unmapped reads.");
        }

        if self.selected_regions_available_flag && self.compute_outside_stats_flag {
            // Todo()!
            // outside_bam_stats
        }

        let overall_time = Instant::now();
        println!("Overall analysis time: {:?}", overall_time - start_time);
    }

    pub fn process_sequence(
        &mut self,
        record: &Record,
        time_to_finalize_and_get_next_window: &mut Duration,
        time_to_analyze_reads_bunch: &mut Duration,
        time_to_finish_reads_bunch: &mut Duration,
        time_to_task_run: &mut Duration,
    ) {
        let bam_stats = self.bam_stats.as_mut().unwrap();
        let current_window = self.current_window.as_ref().unwrap();

        let contig_id = record.tid();
        let contig_option = self.locator.get_contig(contig_id as usize);
        // compute absolute position

        let mut absolute_position = -1;
        if let Some(contig) = contig_option {
            absolute_position = self
                .locator
                .get_absolute_coordinates(contig.name(), record.pos() as i32 + 1);
        }
        // println!("position: {}", absolute_position);

        //compute read size
        let read_size = record.seq_len();
        self.acum_read_size += read_size;
        if read_size > self.max_read_size {
            self.max_read_size = read_size;
        }
        if read_size < self.min_read_size {
            self.min_read_size = read_size;
        }

        if record.is_secondary() {
            self.number_of_secondary_alignments += 1;
            return;
        }

        let novel_read_flag = (record.flags() & Constants::SAM_FLAG_SUPP_ALIGNMENT as u16) == 0;
        if novel_read_flag {
            self.number_of_reads += 1;
        }

        // record.is_quality_check_failed()
        // filter invalid reads
        let read_is_valid: bool = !record.is_quality_check_failed();
        if read_is_valid {
            // accumulate only mapped reads
            if record.is_unmapped() {
                return;
            }

            let mut insert_size = 0;
            if record.is_paired() {
                insert_size = record.insert_size()
            }

            if absolute_position < current_window.get_window_start() {
                panic!("The alignment file is unsorted.\nPlease sort the BAM file by coordinate.");
            }

            // let find_overlappers_start = Instant::now();
            if self.selected_regions_available_flag {
                // todo!()
            } else {
                // read.setAttribute(Constants.READ_IN_REGION, 0);
                if self
                    .bam_stats_collector
                    .update_stats_and_judge_is_dup(record)
                    && self.skip_marked_duplicates_flag
                {
                    self.number_of_duplicates_skip += 1;
                    return;
                }

                if bam_stats
                    .update_read_start_histogram_and_judge_is_detected_dup(absolute_position)
                    && self.skip_detected_duplicates_flag
                {
                    self.number_of_duplicates_skip += 1;
                    return;
                }
                if self.collect_intersecting_paired_end_reads_flag {
                    // todo()! finished
                    self.bam_stats_collector.collect_paired_read_info(record);
                }
                bam_stats.update_insert_size_histogram(insert_size as i32);
            }

            // self.time_to_calc_overlappers = std::time::Instant::now() - find_overlappers_start;

            // finalize current and get next window
            if absolute_position > current_window.get_window_end() {
                // println!("branch 1");
                //analyzeReads(readsBunch);
                self.collect_analysis_results(
                    time_to_analyze_reads_bunch,
                    time_to_finish_reads_bunch,
                    time_to_task_run,
                );

                //finalize
                self.finalize_and_get_next_window(
                    absolute_position,
                    true,
                    time_to_finalize_and_get_next_window,
                );

                if self.selected_regions_available_flag && self.compute_outside_stats_flag {
                    self.current_outside_window
                        .as_mut()
                        .unwrap()
                        .inverse_regions();
                    self.finalize_and_get_next_outside_window(absolute_position, true);
                }
            }

            if self.current_window.is_none() {
                //Some reads are out of reference bounds?
                return;
            }
            // todo!() that clone costs a lot of time
            self.reads_bunch.push(record.clone());

            if self.reads_bunch.len() >= self.num_reads_in_bunch as usize {
                // println!("branch 2");
                // if self.results.len() >= self.max_size_of_task_queue {
                //     self.collect_analysis_results(
                //         time_to_analyze_reads_bunch,
                //         time_to_finish_reads_bunch,
                //         time_to_task_run,
                //     );
                // } else {
                //     self.analyze_reads_bunch(time_to_task_run, time_to_analyze_reads_bunch);
                // }
                self.collect_analysis_results(
                    time_to_analyze_reads_bunch,
                    time_to_finish_reads_bunch,
                    time_to_task_run,
                );
            }

            self.number_of_valid_reads += 1;
        } else {
            self.number_of_problematic_reads += 1;
        }
    }

    fn finalize_and_get_next_window(
        &mut self,
        position: i64,
        detailed_flag: bool,
        time_to_finalize_and_get_next_window: &mut Duration,
    ) {
        // let mut last_window = self.current_window.as_ref().unwrap();
        while position > self.current_window.as_ref().unwrap().get_window_end() {
            let finalize_and_get_next_start = Instant::now();
            self.finalize_window(); // bamstats and current window has been updated in this function
            *time_to_finalize_and_get_next_window += Instant::now() - finalize_and_get_next_start;
            self.current_window = self.next_window(detailed_flag); // update current window
            if self.current_window.is_none() {
                break;
            }
        }
    }

    fn finalize_and_get_next_outside_window(&mut self, position: i64, detailed_flag: bool) {
        // let mut last_window = self.current_window.as_ref().unwrap();
        while position
            > self
                .current_outside_window
                .as_ref()
                .unwrap()
                .get_window_end()
        {
            self.finalize_outside_window(); // bamstats and current window has been updated in this function
            self.current_outside_window = self.next_outside_window(detailed_flag);
            if self.current_outside_window.is_none() {
                break;
            }
        }
    }

    fn finalize_window(&mut self) {
        // todo!() multi thread
        // if (finalizeWindowResult != null) {
        //     // We only run finalization of one window in parallel to prevent to many open windows
        //     finalizeWindowResult.get();
        // }
        let bam_stats = self.bam_stats.as_mut().unwrap();
        let window_start = bam_stats.get_current_window_start();

        // update self.current window
        let window = self.open_windows.get_mut(&window_start).unwrap().clone();
        self.current_window = Some(window);

        // remove the finished window
        self.open_windows.remove(&window_start);
        bam_stats.inc_processed_windows();

        // report progress
        let num_processed_windows = bam_stats.get_number_of_processed_windows();
        if num_processed_windows as usize % self.window_block_size_to_report == 0 {
            println!(
                "Processed {} out of {} windows...",
                { num_processed_windows },
                { self.effective_number_of_window }
            );
        }

        //System.out.println("Time taken to count overlappers: " + timeToCalcOverlappers);
        self.time_to_calc_overlappers = 0;

        // todo!() create a thread to handle Final window task
        let mut finalize_window_task =
            FinalizeWindowTask::new(bam_stats, self.current_window.as_mut().unwrap());
        self.finalize_window_result = finalize_window_task.run();
    }

    fn finalize_outside_window(&mut self) {
        // todo!() multi thread
        // if (finalizeWindowResult != null) {
        //     // We only run finalization of one window in parallel to prevent to many open windows
        //     finalizeWindowResult.get();
        // }
        let bam_stats = self.outside_bam_stats.as_mut().unwrap();
        let window_start = bam_stats.get_current_window_start();

        // update self.current window
        let window = self
            .open_outside_windows
            .get_mut(&window_start)
            .unwrap()
            .clone();
        self.current_outside_window = Some(window);

        // remove the finished window
        self.open_outside_windows.remove(&window_start);
        bam_stats.inc_processed_windows();

        // report progress
        let num_processed_windows = bam_stats.get_number_of_processed_windows();
        if num_processed_windows as usize % self.window_block_size_to_report == 0 {
            println!(
                "Processed {} out of {} windows...",
                { num_processed_windows },
                { self.effective_number_of_window }
            );
        }

        //System.out.println("Time taken to count overlappers: " + timeToCalcOverlappers);
        self.time_to_calc_overlappers = 0;

        // todo!() create a thread to handle Final window task
        let mut finalize_window_task =
            FinalizeWindowTask::new(bam_stats, self.current_outside_window.as_mut().unwrap());
        self.finalize_window_result = finalize_window_task.run();
    }

    fn collect_analysis_results(
        &mut self,
        time_to_analyze_reads_bunch: &mut Duration,
        time_to_finish_reads_bunch: &mut Duration,
        time_to_task_run: &mut Duration,
    ) {
        // start last bunch
        self.analyze_reads_bunch(time_to_task_run, time_to_analyze_reads_bunch);

        let finsih_task_start = Instant::now();
        // wait till all tasks are finished
        for task_result in &self.results {
            let data_sets = task_result.get_read_alignment_data();

            for single_read_data in data_sets {
                let window_start = single_read_data.get_window_start();
                let window = self.open_windows.get_mut(&window_start).unwrap();

                window.add_read_alignment_data(single_read_data);
            }

            self.bam_stats
                .as_mut()
                .unwrap()
                .add_read_stats_data(task_result.get_read_stats_collector());
            if self.selected_regions_available_flag && self.compute_outside_stats_flag {
                let outside_data_sets = task_result.get_out_of_region_reads_data();
                for single_read_data in outside_data_sets {
                    let window_start = single_read_data.get_window_start();
                    let window = self.open_outside_windows.get_mut(&window_start).unwrap();
                    window.add_read_alignment_data(single_read_data);
                }
            }
        }

        self.results.clear();
        *time_to_finish_reads_bunch += Instant::now() - finsih_task_start;
    }

    fn analyze_reads_bunch(
        &mut self,
        time_to_task_run: &mut Duration,
        time_to_analyze_reads_bunch: &mut Duration,
    ) {
        let last_bunch_start = Instant::now();

        let copy_bunch = self.reads_bunch.clone();

        let locator = self.locator.clone();
        let selected_region_flag = self.get_selected_regions_available_flag();
        let compute_outside_stats_flag = self.get_compute_outside_stats_flag();
        let min_homopolymer_size = self.get_min_homopolymer_size();
        // let arc_mutex_self = Arc::new(Mutex::new(self.borrow_mut()));
        let mut task = ProcessBunchOfReadsTask::new(
            copy_bunch,
            selected_region_flag,
            compute_outside_stats_flag,
            min_homopolymer_size,
            locator,
            self,
        );
        // let task_run_start = Instant::now();
        let task_result = task.run(time_to_task_run);
        // *time_to_task_run += Instant::now() - task_run_start;
        // let the_future = executor_service
        //     .submit_async(Box::new(move || {
        //         let mut task_result = task.run();
        //         println!("Long lasting computation finished");
        //         task_result
        //     }))
        //     .expect("Failed to submit function");
        // the_future.get();
        // let mut taks_result = task.run();
        self.results.push(task_result);
        self.reads_bunch.clear();
        *time_to_analyze_reads_bunch += Instant::now() - last_bunch_start;
    }

    fn load_locator(&mut self, header: &Header) {
        let hash_header = self.to_hashmap(&header);
        if hash_header.contains_key("SQ") {
            let sq_records = &hash_header["SQ"];
            for record in sq_records {
                let contig_name = record.get("SN").unwrap();
                let contig_length = record.get("LN").unwrap().parse::<i32>().unwrap();
                self.locator.add_contig(contig_name.clone(), contig_length);
            }
        }
    }

    fn load_program_records(&mut self, header: &Header) {
        let hash_header = self.to_hashmap(&header);
        if hash_header.contains_key("PG") {
            let pg_records = &hash_header["PG"];
            if !pg_records.is_empty() {
                let record = &pg_records[0];
                let program_id = record.get("ID").unwrap();
                if let Some(program_version) = record.get("VN") {
                    self.pg_program = format!("{} ({})", program_id, program_version);
                }
                if let Some(cmd) = record.get("CL") {
                    self.pg_command_string = cmd.clone();
                }
            }
        }
    }

    fn load_reference(&mut self) {
        if self.reference_available_flag {
            // todo!()
            // read custom reference file
            // update reference_size
            // update number_of_reference_config
        } else {
            self.reference_size = self.locator.get_total_size();
            self.number_of_reference_config = self.locator.get_contigs().len();
        }
    }

    fn compute_window_size(&self, reference_size: i64, number_of_windows: i32) -> i32 {
        let window_size = (reference_size as f64 / number_of_windows as f64).floor() as i32;
        if (reference_size as f64 / number_of_windows as f64) > window_size as f64 {
            return window_size + 1;
        }
        window_size
    }

    fn compute_window_positions(&mut self, window_size: i32) -> Vec<i64> {
        let contigs = self.locator.get_contigs();
        let mut window_starts: Vec<i64> = Vec::new();

        let mut start_pos = 1;
        let mut i = 0;
        let num_configs = contigs.len();
        self.chromosome_window_indexes.push(0);
        while start_pos < self.reference_size {
            window_starts.push(start_pos);
            start_pos += window_size as i64;
            while i < num_configs {
                let next_contig_start = contigs[i].end() + 1;
                if start_pos >= next_contig_start {
                    if next_contig_start < self.reference_size {
                        self.chromosome_window_indexes.push(window_starts.len());
                        if start_pos > next_contig_start {
                            window_starts.push(next_contig_start);
                        }
                    }
                    i += 1;
                } else {
                    break;
                }
            }
        }

        window_starts
    }

    fn next_window(&mut self, detailed_flag: bool) -> Option<BamGenomeWindow> {
        // init new current
        let mut new_current = None;

        let bam_stats = self.bam_stats.as_mut().unwrap();
        let num_processed_windows = bam_stats.get_number_of_processed_windows();
        let num_total_windows = bam_stats.get_number_of_windows();

        if num_processed_windows < num_total_windows {
            let window_start = bam_stats.get_window_start(num_processed_windows as usize);
            let window_end = bam_stats.get_current_window_end();
            let reference_size = bam_stats.get_reference_size();
            let current_window_name = bam_stats.get_current_window_name();

            if num_processed_windows < num_total_windows && window_start <= reference_size {
                if self.open_windows.contains_key(&window_start) {
                    let current_window = self
                        .open_windows
                        .get(&window_start)
                        .as_deref()
                        .unwrap()
                        .clone();
                    new_current = Some(current_window);
                } else {
                    bam_stats.increment_initialized_windows();
                    let current_window = self.init_window(
                        current_window_name,
                        window_start,
                        window_end.min(reference_size),
                        detailed_flag,
                    );
                    self.open_windows
                        .insert(window_start, current_window.clone());
                    new_current = Some(current_window);
                }
            }
        }

        new_current
    }

    fn next_outside_window(&mut self, detailed_flag: bool) -> Option<BamGenomeWindow> {
        // init new current
        let mut new_outside_window = None;

        let mut bam_stats = self.outside_bam_stats.as_mut().unwrap();
        let num_processed_windows = bam_stats.get_number_of_processed_windows();
        let num_total_windows = bam_stats.get_number_of_windows();

        if num_processed_windows < num_total_windows {
            let window_start = bam_stats.get_window_start(num_processed_windows as usize);
            let window_end = bam_stats.get_current_window_end();
            let reference_size = bam_stats.get_reference_size();
            let current_window_name = bam_stats.get_current_window_name();

            if num_processed_windows < num_total_windows && window_start <= reference_size {
                if self.open_outside_windows.contains_key(&window_start) {
                    let current_window = self
                        .open_outside_windows
                        .get(&window_start)
                        .unwrap()
                        .clone();
                    new_outside_window = Some(current_window);
                } else {
                    bam_stats.increment_initialized_windows();
                    let current_window = self.init_window(
                        current_window_name,
                        window_start,
                        window_end.min(reference_size),
                        detailed_flag,
                    );
                    self.open_outside_windows
                        .insert(window_start, current_window.clone());
                    new_outside_window = Some(current_window);
                }
            }
        }

        new_outside_window
    }

    fn init_window(
        &mut self,
        name: String,
        window_start: i64,
        window_end: i64,
        detailed_flag: bool,
    ) -> BamGenomeWindow {
        let mut mini_reference: Vec<u8> = vec![];
        if !self.reference.is_empty() {
            mini_reference =
                self.reference[(window_start - 1) as usize..(window_end - 1) as usize].to_vec()
        }

        let window = BamGenomeWindow::new(
            name,
            window_start,
            window_end,
            &mini_reference,
            detailed_flag,
        );

        if self.selected_regions_available_flag {
            // todo!()  calculateRegionsLookUpTableForWindow(w);
        }
        window
    }

    fn get_selected_regions_available_flag(&self) -> bool {
        self.selected_regions_available_flag
    }

    fn get_compute_outside_stats_flag(&self) -> bool {
        self.compute_outside_stats_flag
    }

    fn get_min_homopolymer_size(&self) -> i32 {
        self.min_homopolymer_size
    }

    fn get_locator(&self) -> &GenomeLocator {
        &self.locator
    }

    fn get_bam_stats(&self) -> Option<&BamStats> {
        if let Some(bam_stats) = &self.bam_stats {
            return Some(bam_stats);
        }
        None
    }

    fn get_mut_bam_stats(&mut self) -> Option<&mut BamStats> {
        if let Some(bam_stats) = &mut self.bam_stats {
            return Some(bam_stats);
        }
        None
    }

    fn get_outside_bam_stats(&self) -> Option<&BamStats> {
        if let Some(bam_stats) = &self.outside_bam_stats {
            return Some(bam_stats);
        }
        None
    }

    fn get_mut_outside_bam_stats(&mut self) -> Option<&mut BamStats> {
        if let Some(bam_stats) = &mut self.outside_bam_stats {
            return Some(bam_stats);
        }
        None
    }

    fn get_open_window(&mut self, window_start: i64) -> &BamGenomeWindow {
        if self.open_windows.contains_key(&window_start) {
            self.open_windows.get(&window_start).unwrap()
        } else {
            let bam_stats = self.get_mut_bam_stats().unwrap();
            let num_init_windows = bam_stats.get_number_of_initialized_windows();
            let window_name = bam_stats.get_window_name(num_init_windows as usize);
            let window_end = bam_stats.get_window_end(num_init_windows as usize);
            let reference_size = bam_stats.get_reference_size();
            bam_stats.increment_initialized_windows();
            let new_window = self.init_window(
                window_name,
                window_start,
                window_end.min(reference_size),
                true,
            );

            self.open_windows.insert(window_start, new_window);

            self.open_windows.get(&window_start).unwrap()
        }
    }

    fn get_outside_open_window(&mut self, window_start: i64) -> &BamGenomeWindow {
        if self.open_outside_windows.contains_key(&window_start) {
            self.open_outside_windows.get(&window_start).unwrap()
        } else {
            let bam_stats = self.get_mut_outside_bam_stats().unwrap();
            let num_init_windows = bam_stats.get_number_of_initialized_windows();
            let window_name = bam_stats.get_window_name(num_init_windows as usize);
            let window_end = bam_stats.get_window_end(num_init_windows as usize);
            let reference_size = bam_stats.get_reference_size();
            bam_stats.increment_initialized_windows();
            let new_window = self.init_window(
                window_name,
                window_start,
                window_end.min(reference_size),
                true,
            );

            self.open_outside_windows.insert(window_start, new_window);

            self.open_outside_windows.get(&window_start).unwrap()
        }
    }

    // fn get_current_window(&self) -> &Option<BamGenomeWindow> {
    //     &self.current_window
    // }

    fn inc_number_of_reads_with_start_greater_then_end(&mut self) {
        self.number_of_reads_with_start_greater_than_end += 1;
    }

    fn get_current_window(&self) -> Option<&BamGenomeWindow> {
        if let Some(window) = &self.current_window {
            return Some(window);
        }
        None
    }
}

pub struct InputParameters {}

impl InputParameters {}

pub struct CoverageAcrossReference {}

impl CoverageAcrossReference {}

pub struct CovergaeHistogram {}

impl CovergaeHistogram {}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct InputInfo {
    bam_file: String,
    analysis_date: String,
    number_of_windows: i32,
    size_of_homopolymer: i32,
    program: String,
    command_line: String,
    skip_duplicate_alignments: String,
    analyze_overlapping_paired_end_reads: String,
}

impl InputInfo {
    pub fn new() -> Self {
        Self {
            bam_file: "".to_string(),
            analysis_date: "".to_string(),
            number_of_windows: 0,
            size_of_homopolymer: 0,
            program: "".to_string(),
            command_line: "".to_string(),
            skip_duplicate_alignments: "".to_string(),
            analyze_overlapping_paired_end_reads: "".to_string(),
        }
    }

    fn extract(&mut self, bam_stats_analysis: &BamStatsAnalysis) {
        self.bam_file = bam_stats_analysis.bam_file.clone();
        // todo!() date
        // let format = format_description!("[year]-[month]-[day]");
        // let current = Instant::now();
        // let date = Date.now
        self.number_of_windows = bam_stats_analysis.number_of_windows;
        self.size_of_homopolymer = bam_stats_analysis.min_homopolymer_size;

        self.program = bam_stats_analysis.pg_program.clone();
        self.command_line = bam_stats_analysis.pg_command_string.clone();

        let mut skip_duplicate_aln_status = "false";
        if bam_stats_analysis.skip_detected_duplicates_flag {
            skip_duplicate_aln_status = "yes (only estimated)";
        } else if bam_stats_analysis.skip_marked_duplicates_flag {
            skip_duplicate_aln_status = "yes (only flagged)";
        } else if bam_stats_analysis.skip_detected_duplicates_flag
            && bam_stats_analysis.skip_marked_duplicates_flag
        {
            skip_duplicate_aln_status = "yes (both flagged and estimated)";
        }
        self.skip_duplicate_alignments = skip_duplicate_aln_status.to_string();
        self.analyze_overlapping_paired_end_reads = bam_stats_analysis
            .collect_intersecting_paired_end_reads_flag
            .to_string();
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct BasicStats {
    // gloals
    reference_size: usize,
    number_of_reads: usize,
    mapped_reads: String,
    unmapped_reads: String,
    mapped_paired_reads: String,
    first_mapped_read_in_pair: String,
    second_mapped_read_in_pair: String,
    both_mapped_reads_in_pair: String,
    singletons_mapped_reads_in_pair: String,
    number_of_secondary_alignments: usize,
    number_of_supplementary_alignments: String,
    read_min_max_mean_length: String,
    flagged_duplicated_reads: String,
    estimated_duplicated_reads: String,
    duplication_rate: String,
    clipped_reads: String,
    duplicated_reads_skipped: String,

    // ATCG content
    number_and_percent_of_a: String,
    number_and_percent_of_t: String,
    number_and_percent_of_c: String,
    number_and_percent_of_g: String,
    number_and_percent_of_n: String,
    gc_percent: String,

    // coverage
    mean_coverage: String,
    standard_deviation_of_coverage: String,
    mean_coverage_with_ignored_overlap: String,

    //mapping quality
    mean_mapping_quality: String,
    mean_insert_size: String,

    //insert size
    standard_deviation_of_insert_size: String,
    p25_median_p75_insert_size: String,
    general_error_rate: String,

    // mismatches and indels
    number_of_mismatches: usize,
    number_of_insertions: usize,
    percent_of_mapped_reads_with_at_least_one_insertion: String,
    number_of_deletions: usize,
    percent_of_mapped_reads_with_at_least_one_deletion: String,
    percent_of_homopolymer_indels: String,

    // Chromosome
    chromosome_stats: Vec<ChromosomeInfo>,
}

impl BasicStats {
    pub fn new() -> Self {
        Self {
            reference_size: 0,
            number_of_reads: 0,
            mapped_reads: "".to_string(),
            unmapped_reads: "".to_string(),
            mapped_paired_reads: "".to_string(),
            first_mapped_read_in_pair: "".to_string(),
            second_mapped_read_in_pair: "".to_string(),
            both_mapped_reads_in_pair: "".to_string(),
            singletons_mapped_reads_in_pair: "".to_string(),
            number_of_secondary_alignments: 0,
            number_of_supplementary_alignments: "".to_string(),
            read_min_max_mean_length: "".to_string(),
            flagged_duplicated_reads: "".to_string(),
            estimated_duplicated_reads: "".to_string(),
            duplication_rate: "".to_string(),
            clipped_reads: "".to_string(),
            duplicated_reads_skipped: "".to_string(),

            number_and_percent_of_a: "".to_string(),
            number_and_percent_of_t: "".to_string(),
            number_and_percent_of_c: "".to_string(),
            number_and_percent_of_g: "".to_string(),
            number_and_percent_of_n: "".to_string(),
            gc_percent: "".to_string(),
            mean_coverage: "".to_string(),
            standard_deviation_of_coverage: "".to_string(),
            mean_coverage_with_ignored_overlap: "".to_string(),
            mean_mapping_quality: "".to_string(),
            mean_insert_size: "".to_string(),
            standard_deviation_of_insert_size: "".to_string(),
            p25_median_p75_insert_size: "".to_string(),
            general_error_rate: "".to_string(),
            number_of_mismatches: 0,
            number_of_insertions: 0,
            percent_of_mapped_reads_with_at_least_one_insertion: "".to_string(),
            number_of_deletions: 0,
            percent_of_mapped_reads_with_at_least_one_deletion: "".to_string(),
            percent_of_homopolymer_indels: "".to_string(),

            chromosome_stats: vec![],
        }
    }

    fn extract(&mut self, bam_stats: &BamStats) {
        self.update_globals(bam_stats);
        // todo!() region analyze
        // if bam_stats.num_selected_regions > 0 {

        // }

        self.update_base_content(bam_stats);
        self.update_coverage(bam_stats);
        self.update_mapping_quality(bam_stats);
        self.update_insert_size(bam_stats);
        self.update_mismatches_and_indels(bam_stats);
        self.update_chromosome_info(bam_stats);
    }

    fn update_globals(&mut self, bam_stats: &BamStats) {
        self.reference_size = bam_stats.reference_size as usize;
        self.number_of_reads = bam_stats.number_of_reads as usize;

        let mapped_reads_num = bam_stats.number_of_mapped_reads;
        let mapped_reads_percent = mapped_reads_num as f64 / self.number_of_reads as f64;
        self.mapped_reads = format!(
            "{} / {:.2}%",
            mapped_reads_num,
            mapped_reads_percent * 100.0
        );

        self.unmapped_reads = format!(
            "{} / {:.2}%",
            self.number_of_reads - mapped_reads_num as usize,
            100.0 - mapped_reads_percent * 100.0
        );

        let mapped_paired_reads_num = bam_stats.number_of_paired_reads;
        let mapped_paired_reads_percent =
            mapped_paired_reads_num as f64 / self.number_of_reads as f64;
        self.mapped_paired_reads = format!(
            "{} / {:.2}%",
            mapped_paired_reads_num,
            mapped_paired_reads_percent * 100.0
        );

        if mapped_paired_reads_num > 0 {
            let number_of_mapped_first_of_pair = bam_stats.number_of_mapped_first_of_pair;
            let number_of_mapped_second_of_pair = bam_stats.number_of_mapped_second_of_pair;
            let number_of_singletons = bam_stats.number_of_singletons;
            let number_of_mapped_both_of_pair = mapped_paired_reads_num - number_of_singletons;
            self.first_mapped_read_in_pair = format!(
                "{} / {:.2}%",
                number_of_mapped_first_of_pair,
                number_of_mapped_first_of_pair as f64 / self.number_of_reads as f64 * 100.0
            );
            self.second_mapped_read_in_pair = format!(
                "{} / {:.2}%",
                number_of_mapped_second_of_pair,
                number_of_mapped_second_of_pair as f64 / self.number_of_reads as f64 * 100.0
            );
            self.both_mapped_reads_in_pair = format!(
                "{} / {:.2}%",
                number_of_mapped_both_of_pair,
                number_of_mapped_both_of_pair as f64 / self.number_of_reads as f64 * 100.0
            );

            self.singletons_mapped_reads_in_pair = format!(
                "{} / {:.2}%",
                number_of_singletons,
                number_of_singletons as f64 / self.number_of_reads as f64 * 100.0
            );
        }

        self.number_of_secondary_alignments = bam_stats.number_of_secondary_alignments as usize;
        if bam_stats.number_of_supp_alignments > 0 {
            self.number_of_supplementary_alignments = format!(
                "{} / {:.2}",
                bam_stats.number_of_supp_alignments,
                bam_stats.number_of_supp_alignments as f64 / self.number_of_reads as f64 * 100.0
            );
        }

        self.read_min_max_mean_length = format!(
            "{} / {} / {:.2}",
            bam_stats.read_min_size, bam_stats.read_max_size, bam_stats.read_mean_size
        );

        if bam_stats.num_selected_regions == 0 {
            if bam_stats.report_overlapping_read_pairs {
                // todo!()
                // add overlapping
            }
            let show_marked_duplicates = bam_stats.num_detected_duplicate_reads > 0
                || bam_stats.skip_duplicates_mode == SkipDuplicatesMode::ONLY_MARKED_DUPLICATES
                || bam_stats.skip_duplicates_mode == SkipDuplicatesMode::BOTH;
            if show_marked_duplicates {
                self.flagged_duplicated_reads = format!(
                    "{} / {:.2}%",
                    bam_stats.num_detected_duplicate_reads,
                    bam_stats.num_detected_duplicate_reads as f64 / self.number_of_reads as f64
                        * 100.0
                )
            }
            let show_estimate_duplicates = bam_stats.num_detected_duplicate_reads == 0
                || bam_stats.skip_duplicates_mode == SkipDuplicatesMode::ONLY_DETECTED_DUPLICATES
                || bam_stats.skip_duplicates_mode == SkipDuplicatesMode::BOTH;
            if show_estimate_duplicates {
                self.estimated_duplicated_reads = format!(
                    "{} / {:.2}%",
                    bam_stats.num_estimated_duplicate_reads,
                    bam_stats.num_estimated_duplicate_reads as f64 / self.number_of_reads as f64
                        * 100.0
                );
                self.duplication_rate = format!("{:.2}%", bam_stats.duplication_rate);
            }
        }

        self.clipped_reads = format!(
            "{} / {:.2}",
            bam_stats.num_clipped_reads,
            bam_stats.num_clipped_reads as f64 / self.number_of_reads as f64 * 100.0
        );

        if bam_stats.num_duplicates_skipped > 0 {
            self.duplicated_reads_skipped = format!(
                "{} / {:.2}%",
                bam_stats.num_duplicates_skipped,
                bam_stats.num_duplicates_skipped as f64 / self.number_of_reads as f64 * 100.0
            );
        }
    }

    fn update_base_content(&mut self, bam_stats: &BamStats) {
        self.number_and_percent_of_a = format!(
            "{} / {:.2}%",
            bam_stats.number_of_as, bam_stats.mean_a_relative_content
        );
        self.number_and_percent_of_t = format!(
            "{} / {:.2}%",
            bam_stats.number_of_ts, bam_stats.mean_t_relative_content
        );
        self.number_and_percent_of_c = format!(
            "{} / {:.2}%",
            bam_stats.number_of_cs, bam_stats.mean_c_relative_content
        );
        self.number_and_percent_of_g = format!(
            "{} / {:.2}%",
            bam_stats.number_of_gs, bam_stats.mean_g_relative_content
        );
        self.number_and_percent_of_n = format!(
            "{} / {:.2}%",
            bam_stats.number_of_ns, bam_stats.mean_n_relative_content
        );

        self.gc_percent = format!("{:.2}%", bam_stats.mean_gc_relative_content);
    }

    fn update_coverage(&mut self, bam_stats: &BamStats) {
        self.mean_coverage = format!("{:.4}", bam_stats.mean_coverage);
        self.standard_deviation_of_coverage = format!("{:.4}", bam_stats.std_coverage);

        if bam_stats.num_overlapping_read_pairs > 0 {
            self.mean_coverage_with_ignored_overlap =
                format!("{:.4}", bam_stats.adapted_mean_coverage);
        }
    }

    fn update_mapping_quality(&mut self, bam_stats: &BamStats) {
        self.mean_mapping_quality = format!("{:.2}", bam_stats.mean_mapping_quality_per_window);
    }

    fn update_insert_size(&mut self, bam_stats: &BamStats) {
        if bam_stats.mean_insert_size != 0.0 {
            self.mean_insert_size = format!("{:.2}", bam_stats.mean_insert_size);
            self.standard_deviation_of_insert_size = format!("{:.2}", bam_stats.std_insert_size);
            self.p25_median_p75_insert_size = format!(
                "{} / {} / {}",
                bam_stats.p25_insert_size, bam_stats.median_insert_size, bam_stats.p75_insert_size
            );
        }
    }

    fn update_mismatches_and_indels(&mut self, bam_stats: &BamStats) {
        let num_indels = bam_stats.num_insertions + bam_stats.num_deletions;
        let error_rate = bam_stats.get_error_rate();
        if num_indels > 0 || bam_stats.num_mismatches > 0 || error_rate > 0.0 {
            if error_rate > 0.0 {
                self.general_error_rate = format!("{:.2}%", error_rate * 100.0);
            }
            if bam_stats.num_mismatches > 0 {
                self.number_of_mismatches = bam_stats.num_mismatches as usize;
            }
            if num_indels > 0 {
                self.number_of_insertions = bam_stats.num_insertions as usize;
                self.percent_of_mapped_reads_with_at_least_one_insertion = format!(
                    "{:.2}%",
                    bam_stats.num_reads_with_insertion as f64
                        / bam_stats.number_of_mapped_reads as f64
                        * 100.0
                );
                self.number_of_deletions = bam_stats.num_deletions as usize;
                self.percent_of_mapped_reads_with_at_least_one_deletion = format!(
                    "{:.2}%",
                    bam_stats.num_reads_with_deletion as f64
                        / bam_stats.number_of_mapped_reads as f64
                        * 100.0
                );
                self.percent_of_homopolymer_indels = format!(
                    "{:.2}%",
                    bam_stats.get_homopolymer_indels_fraction() * 100.0
                );
            }
        }
    }

    fn update_chromosome_info(&mut self, bam_stats: &BamStats) {
        self.chromosome_stats = bam_stats.chromosome_stats.clone();
    }
}
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct CoverageStats {
    // coverage across reference
    coverage_across_reference_x_label: String,
    coverage_across_reference_x_data: Vec<f64>,
    coverage_across_reference_y_label: String,
    coverage_across_reference_y_data: Vec<f64>,
    coverage_across_reference_deviation_std_y_data: Vec<f64>,
    // balanced coverage histogram
    balanced_coverage_bar_names_x_label: String,
    balanced_coverage_bar_names_x_data: Vec<String>,
    balance_coverage_histogram_y_label: String,
    balance_coverage_histogram_y_data: XYVector,
    // Coverage Histogram 0-50x
    ranged_coverage_histogram_x_label: String,
    ranged_coverage_histogram_x_data: Vec<String>,
    ranged_coverage_histogram_y_label: String,
    ranged_coverage_histogram_y_data: XYVector,
    // Genome Fraction Coverage
    genome_fraction_coverage_x_label: String,
    genome_fraction_coverage_x_data: Vec<String>,
    genome_fraction_coverage_y_label: String,
    genome_fraction_coverage_y_data: XYVector,
}

impl CoverageStats {
    pub fn new() -> Self {
        Self {
            // coverage across reference
            coverage_across_reference_x_label: "".to_string(),
            coverage_across_reference_x_data: vec![],
            coverage_across_reference_y_label: "".to_string(),
            coverage_across_reference_y_data: vec![],
            coverage_across_reference_deviation_std_y_data: vec![],
            // balanced coverage histogram
            balanced_coverage_bar_names_x_label: "".to_string(),
            balanced_coverage_bar_names_x_data: vec![],
            balance_coverage_histogram_y_label: "".to_string(),
            balance_coverage_histogram_y_data: XYVector::new(),
            // Coverage Histogram 0-50x
            ranged_coverage_histogram_x_label: "".to_string(),
            ranged_coverage_histogram_x_data: vec![],
            ranged_coverage_histogram_y_label: "".to_string(),
            ranged_coverage_histogram_y_data: XYVector::new(),
            // Genome Fraction Coverage
            genome_fraction_coverage_x_label: "".to_string(),
            genome_fraction_coverage_x_data: vec![],
            genome_fraction_coverage_y_label: "".to_string(),
            genome_fraction_coverage_y_data: XYVector::new(),
        }
    }

    fn extract(&mut self, bam_stats: &BamStats) {
        // coverage across reference
        self.coverage_across_reference_x_label = "Position (bp)".to_string();
        self.coverage_across_reference_x_data =
            vec![0.0; bam_stats.get_number_of_windows() as usize];
        for i in 0..bam_stats.get_number_of_windows() as usize {
            self.coverage_across_reference_x_data[i] =
                (bam_stats.get_window_start(i) + bam_stats.get_window_end(i)) as f64 / 2.0;
        }

        self.coverage_across_reference_y_label = "Coverage (X)".to_string();
        self.coverage_across_reference_y_data = bam_stats.coverage_across_reference.clone();
        self.coverage_across_reference_deviation_std_y_data =
            bam_stats.std_coverage_across_reference.clone();

        // balanced coverage histogram
        self.balance_coverage_histogram_y_label = "Number of genomic locations".to_string();
        self.balance_coverage_histogram_y_data = bam_stats.balanced_coverage_histogram.clone();

        self.balanced_coverage_bar_names_x_label = "Coverage (X)".to_string();
        for i in 0..bam_stats.balanced_coverage_bar_names.len() {
            self.balanced_coverage_bar_names_x_data
                .push(bam_stats.balanced_coverage_bar_names[&(i as i64)].clone());
        }

        // Coverage Histogram 0-50x
        let min_coverage = bam_stats.coverage_histogram.get(0).unwrap().get_x();
        if min_coverage < 40.0 {
            // This histogram only makes sense for low coverage samples
            self.ranged_coverage_histogram_x_label = "Coverage (X)".to_string();
            self.ranged_coverage_histogram_y_label = "Number of genomic locations".to_string();

            for index in 0..bam_stats.coverage_histogram.items.len() {
                if index > 50 {
                    break;
                }
                let item = &bam_stats.coverage_histogram.items[index];
                self.ranged_coverage_histogram_y_data
                    .add_item(XYItem::new(item.xy_item.get_x(), item.xy_item.get_y()));
                self.ranged_coverage_histogram_x_data
                    .push(item.xy_item.get_x().to_string());
            }
        }

        // Genome Fraction Coverage
        self.genome_fraction_coverage_y_data = bam_stats.coverage_quotes.clone();
        self.genome_fraction_coverage_x_label = "Coverage (X)".to_string();
        self.genome_fraction_coverage_y_label = "Fraction of reference (%)".to_string();

        for i in 0..self.genome_fraction_coverage_y_data.items.len() {
            self.genome_fraction_coverage_x_data
                .push((i + 1).to_string());
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct DuplicationStats {
    duplication_rate_histogram_x_label: String,
    duplication_rate_histogram_x_data: Vec<String>,
    duplication_rate_histogram_y_label: String,
    duplication_rate_histogram_y_data: XYVector,
}

impl DuplicationStats {
    pub fn new() -> Self {
        Self {
            duplication_rate_histogram_x_label: "".to_string(),
            duplication_rate_histogram_x_data: vec![],
            duplication_rate_histogram_y_label: "".to_string(),
            duplication_rate_histogram_y_data: XYVector::new(),
        }
    }

    fn extract(&mut self, bam_stats: &BamStats) {
        self.duplication_rate_histogram_x_label = "Duplication rate".to_string();
        self.duplication_rate_histogram_y_label = "Number of loci".to_string();

        self.duplication_rate_histogram_y_data = bam_stats.unique_read_starts_histogram.clone();

        for i in 0..self.duplication_rate_histogram_y_data.items.len() {
            self.duplication_rate_histogram_x_data
                .push((i + 1).to_string());
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct MappedReadsNucleotideContent {
    reads_nuleotide_a_content_xy_data: XYVector,
    reads_nuleotide_t_content_xy_data: XYVector,
    reads_nuleotide_c_content_xy_data: XYVector,
    reads_nuleotide_g_content_xy_data: XYVector,
    reads_nuleotide_n_content_xy_data: XYVector,
}

impl MappedReadsNucleotideContent {
    pub fn new() -> Self {
        Self {
            reads_nuleotide_a_content_xy_data: XYVector::new(),
            reads_nuleotide_t_content_xy_data: XYVector::new(),
            reads_nuleotide_c_content_xy_data: XYVector::new(),
            reads_nuleotide_g_content_xy_data: XYVector::new(),
            reads_nuleotide_n_content_xy_data: XYVector::new(),
        }
    }

    fn extract(&mut self, bam_stats: &BamStats) {
        self.reads_nuleotide_a_content_xy_data = bam_stats.reads_as_histogram.clone();
        self.reads_nuleotide_t_content_xy_data = bam_stats.reads_ts_histogram.clone();
        self.reads_nuleotide_c_content_xy_data = bam_stats.reads_cs_histogram.clone();
        self.reads_nuleotide_g_content_xy_data = bam_stats.reads_gs_histogram.clone();
        self.reads_nuleotide_n_content_xy_data = bam_stats.reads_ns_histogram.clone();
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct MappedReadsGCContent {
    reads_gc_content_histogram_x_label: String,
    reads_gc_content_histogram_y_label: String,
    reads_gc_content_histogram_xy_data: XYVector,
}

impl MappedReadsGCContent {
    pub fn new() -> Self {
        Self {
            reads_gc_content_histogram_x_label: "".to_string(),
            reads_gc_content_histogram_y_label: "".to_string(),
            reads_gc_content_histogram_xy_data: XYVector::new(),
        }
    }

    fn extract(&mut self, bam_stats: &BamStats) {
        // todo!() comparison of reference genome gc content
        self.reads_gc_content_histogram_x_label = "GC Content (%)".to_string();
        self.reads_gc_content_histogram_y_label = "Fraction of reads".to_string();

        self.reads_gc_content_histogram_xy_data = bam_stats.get_gc_content_histogram();
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct ClippingProfileStats {
    clipping_profile_x_label: String,
    clipping_profile_y_label: String,
    clipping_profile_xy_data: XYVector,
}
impl ClippingProfileStats {
    pub fn new() -> Self {
        Self {
            clipping_profile_x_label: "".to_string(),
            clipping_profile_y_label: "".to_string(),
            clipping_profile_xy_data: XYVector::new(),
        }
    }

    fn extract(&mut self, bam_stats: &BamStats) {
        self.clipping_profile_x_label = "Read position (bp)".to_string();
        self.clipping_profile_y_label = " Clipped bases (%)".to_string();
        self.clipping_profile_xy_data = bam_stats.reads_clipping_profile_histogram.clone();
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct HomopolymerIndelStats {
    poly_a: i32,
    poly_c: i32,
    poly_g: i32,
    poly_t: i32,
    poly_n: i32,
    non_poly: i32,
}
impl HomopolymerIndelStats {
    pub fn new() -> Self {
        Self {
            poly_a: 0,
            poly_c: 0,
            poly_g: 0,
            poly_t: 0,
            poly_n: 0,
            non_poly: 0,
        }
    }

    fn extract(&mut self, bam_stats: &BamStats) {
        self.poly_a = bam_stats.homopolymer_indels_data[0];
        self.poly_c = bam_stats.homopolymer_indels_data[1];
        self.poly_g = bam_stats.homopolymer_indels_data[2];
        self.poly_t = bam_stats.homopolymer_indels_data[3];
        self.poly_n = bam_stats.homopolymer_indels_data[4];
        self.non_poly = bam_stats.homopolymer_indels_data[5];
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct MappingQualityStats {
    // mapping quality across reference
    mapping_quality_across_reference_x_label: String,
    mapping_quality_across_reference_x_data: Vec<f64>,
    mapping_quality_across_reference_y_label: String,
    mapping_quality_across_reference_y_data: Vec<f64>,
    // mapping quality histogram
    mapping_quality_histogram_x_label: String,
    mapping_quality_histogram_y_label: String,
    mapping_quality_histogram_xy_data: XYVector,
}
impl MappingQualityStats {
    pub fn new() -> Self {
        Self {
            // mapping quality across reference
            mapping_quality_across_reference_x_label: "".to_string(),
            mapping_quality_across_reference_x_data: vec![],
            mapping_quality_across_reference_y_label: "".to_string(),
            mapping_quality_across_reference_y_data: vec![],
            // mapping quality histogram
            mapping_quality_histogram_x_label: "".to_string(),
            mapping_quality_histogram_y_label: "".to_string(),
            mapping_quality_histogram_xy_data: XYVector::new(),
        }
    }

    fn extract(&mut self, bam_stats: &BamStats) {
        // mapping quality across reference
        self.mapping_quality_across_reference_x_label = "Position (bp)".to_string();
        self.mapping_quality_across_reference_y_label = "Mapping quality".to_string();
        self.mapping_quality_across_reference_y_data =
            bam_stats.mapping_quality_across_reference.clone();
        self.mapping_quality_across_reference_x_data =
            vec![0.0; bam_stats.get_number_of_windows() as usize];
        for i in 0..bam_stats.get_number_of_windows() as usize {
            self.mapping_quality_across_reference_x_data[i] =
                (bam_stats.get_window_start(i) + bam_stats.get_window_end(i)) as f64 / 2.0;
        }

        // mapping quality histogram
        self.mapping_quality_histogram_x_label = "Mapping quality".to_string();
        self.mapping_quality_histogram_y_label = "Number of genomic locations".to_string();
        self.mapping_quality_histogram_xy_data = bam_stats.mapping_quality_histogram.clone();
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct InsertSizeStats {
    // insert size across reference
    insert_size_across_reference_x_label: String,
    insert_size_across_reference_x_data: Vec<f64>,
    insert_size_across_reference_y_label: String,
    insert_size_across_reference_y_data: Vec<f64>,
    // insert size histogram
    insert_size_histogram_x_label: String,
    insert_size_histogram_y_label: String,
    insert_size_histogram_xy_data: XYVector,
}
impl InsertSizeStats {
    pub fn new() -> Self {
        Self {
            // insert size across reference
            insert_size_across_reference_x_label: "".to_string(),
            insert_size_across_reference_x_data: vec![],
            insert_size_across_reference_y_label: "".to_string(),
            insert_size_across_reference_y_data: vec![],
            // insert size histogram
            insert_size_histogram_x_label: "".to_string(),
            insert_size_histogram_y_label: "".to_string(),
            insert_size_histogram_xy_data: XYVector::new(),
        }
    }

    fn extract(&mut self, bam_stats: &BamStats) {
        // insert size across reference
        self.insert_size_across_reference_x_label = "Position (bp)".to_string();
        self.insert_size_across_reference_y_label = "Insert size (bp)".to_string();
        self.insert_size_across_reference_y_data = bam_stats.insert_size_across_reference.clone();
        self.insert_size_across_reference_x_data =
            vec![0.0; bam_stats.get_number_of_windows() as usize];
        for i in 0..bam_stats.get_number_of_windows() as usize {
            self.insert_size_across_reference_x_data[i] =
                (bam_stats.get_window_start(i) + bam_stats.get_window_end(i)) as f64 / 2.0;
        }

        // insert size histogram
        self.insert_size_histogram_x_label = "Insert size (bp)".to_string();
        self.insert_size_histogram_y_label = "Number of reads".to_string();
        self.insert_size_histogram_xy_data = bam_stats.insert_size_histogram.clone();
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Qualimap {
    pub input_info: InputInfo,
    pub basic_stats: BasicStats,
    pub coverage_stats: CoverageStats,
    pub duplication_stats: DuplicationStats,
    pub mapping_read_nucleotide_content: MappedReadsNucleotideContent,
    pub mapping_read_gc_content: MappedReadsGCContent,
    pub clipping_profile_stats: ClippingProfileStats,
    pub homopolymer_indel_stats: HomopolymerIndelStats,
    pub mapping_quality_stats: MappingQualityStats,
    pub insert_size_stats: InsertSizeStats,
}

impl Qualimap {
    pub fn new() -> Self {
        Self {
            input_info: InputInfo::new(),
            basic_stats: BasicStats::new(),
            coverage_stats: CoverageStats::new(),
            duplication_stats: DuplicationStats::new(),
            mapping_read_nucleotide_content: MappedReadsNucleotideContent::new(),
            mapping_read_gc_content: MappedReadsGCContent::new(),
            clipping_profile_stats: ClippingProfileStats::new(),
            homopolymer_indel_stats: HomopolymerIndelStats::new(),
            mapping_quality_stats: MappingQualityStats::new(),
            insert_size_stats: InsertSizeStats::new(),
        }
    }

    pub fn run(&mut self, bam_file: String) {
        let mut bam_stats_analysis = BamStatsAnalysis::new(bam_file);
        println!("Running BAM file analysis...");
        bam_stats_analysis.run();
        println!("End of bam qc");
        println!("Computing report...");

        let bam_stats = bam_stats_analysis.get_bam_stats().unwrap();
        self.input_info.extract(&bam_stats_analysis);
        self.basic_stats.extract(bam_stats);
        self.coverage_stats.extract(bam_stats);
        self.duplication_stats.extract(bam_stats);
        self.mapping_read_nucleotide_content.extract(bam_stats);
        self.mapping_read_gc_content.extract(bam_stats);

        if bam_stats.clipping_is_present() {
            self.clipping_profile_stats.extract(bam_stats);
        }

        self.homopolymer_indel_stats.extract(bam_stats);
        self.mapping_quality_stats.extract(bam_stats);

        if bam_stats_analysis.is_paired_data_flag {
            self.insert_size_stats.extract(bam_stats);
        }
    }
}
