#![feature(mutex_unlock)]
use crate::qc::config::bam_config::{QualimapConfig, QualimapConstants};
use crate::qc::util::{BamRecord, RegionOverlapLookupTable};
use bit_vec::BitVec;
use executor_service::{ExecutorService, Executors, Future};
use noodles::{bed, gff, gtf};
extern crate num_cpus;
use csv::WriterBuilder;
use lazy_static::lazy_static;
use linear_map::LinearMap;
use math::stats;
use regex::Regex;
use rust_htslib::bam::record::{Aux, AuxArray, Cigar};
use rust_htslib::bam::{Header, Read, Reader, Record};
use rust_htslib::tpool;
use serde::{Deserialize, Serialize};
use std::error::Error;
use std::fmt::Debug;
use std::fs::File;
use std::io::prelude::*;
use std::io::{self, BufReader};
use std::path::Path;
use std::sync::Mutex;
use std::{
    collections::HashMap,
    f64::consts::{E, PI},
    vec,
};
use time::{ Duration, Instant};

lazy_static! {
    static ref OPEN_WINDOWS: Mutex<HashMap<i64, BamGenomeWindow>> = {
        let mut map: HashMap<i64, BamGenomeWindow> = HashMap::new();
        Mutex::new(map)
    };
    static ref OPEN_OUTSIDE_WINDOWS: Mutex<HashMap<i64, BamGenomeWindow>> = {
        let mut map: HashMap<i64, BamGenomeWindow> = HashMap::new();
        Mutex::new(map)
    };
    static ref GLOBAL_EXECUTOR: Mutex<ExecutorService<Box<dyn FnOnce() -> TaskResult + Send + 'static>, TaskResult>> = Mutex::new(Executors::new_fixed_thread_pool(4).expect("Failed to create the thread pool")); // 这里的数字表示线程池的大小
}
static mut NUMBER_OF_INITIALIZED_WINDOWS: i32 = 0;
static mut NUMBER_OF_INITIALIZED_OUTSIDE_WINDOWS: i32 = 0;
static mut NUMBER_OF_PROCESSED_WINDOWS: i32 = 0;
static mut NUMBER_OF_PROCESSED_OUTSIDE_WINDOWS: i32 = 0;
static mut NUMBER_OF_TOTAL_WINDOWS: i32 = 0;
static mut WINDOW_NANMES: Vec<String> = vec![];
static mut OUTSIDE_WINDOW_NANMES: Vec<String> = vec![];
static mut WINDOW_STARTS: Vec<i64> = vec![];
static mut WINDOW_ENDS: Vec<i64> = vec![];
static mut CUSTOM_REFERENCES: Vec<u8> = vec![];
static mut SELECTED_REGION_STARTS: Vec<i64> = vec![];
static mut SELECTED_REGION_ENDS: Vec<i64> = vec![];

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
    max_read_starts_per_position: i32,
    current_read_start_position: i64,
    read_start_counter: i32,
    read_starts_histogram: Vec<i64>,
}

impl ReadStartsHistogram {
    fn new() -> Self {
        Self {
            max_read_starts_per_position: QualimapConstants::DEFAULT_DUPL_RATE_HIST_MAX,
            current_read_start_position: -1,
            read_start_counter: 1,
            read_starts_histogram: vec![
                0;
                QualimapConstants::DEFAULT_DUPL_RATE_HIST_MAX as usize + 1
            ],
        }
    }

    pub fn set_max_read_starts_per_position(&mut self, num: i32) {
        self.max_read_starts_per_position = num;
        self.read_starts_histogram = vec![0; num as usize + 1];
    }

    pub fn get_max_read_starts_per_position(&self) -> i32 {
        self.max_read_starts_per_position
    }
    pub fn update(&mut self, position: i64) -> bool {
        if position == self.current_read_start_position {
            self.read_start_counter += 1;
        } else {
            let hist_pos = if self.read_start_counter < self.max_read_starts_per_position {
                self.read_start_counter
            } else {
                self.max_read_starts_per_position
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
            // array.reserve(expected_size - size+1);
            // for _ in 0..new_size - old_size {
            //     array.push(0);
            // }
            for _ in 0..expected_size - size + 1 {
                // tothink!() may be it's expected_size-size ,original is wrong
                array.push(0);
            }
        }
    }

    fn set_window_references(&mut self, prefix: &str, window_positions: Vec<i64>) {
        unsafe {
            for i in 0..NUMBER_OF_TOTAL_WINDOWS as usize {
                WINDOW_NANMES.push(format!("{}_{}", prefix, i + 1));
                WINDOW_STARTS.push(window_positions[i]);

                if i + 1 == NUMBER_OF_TOTAL_WINDOWS as usize {
                    WINDOW_ENDS.push(self.reference_size);
                } else {
                    WINDOW_ENDS.push(window_positions[i + 1] - 1);
                }
            }
        }
    }

    fn set_outside_window_references(&mut self, prefix: &str, window_positions: Vec<i64>) {
        unsafe {
            for i in 0..NUMBER_OF_TOTAL_WINDOWS as usize {
                OUTSIDE_WINDOW_NANMES.push(format!("{}_{}", prefix, i + 1));
                if WINDOW_STARTS.is_empty() && WINDOW_ENDS.is_empty() {
                    WINDOW_STARTS.push(window_positions[i]);

                    if i + 1 == NUMBER_OF_TOTAL_WINDOWS as usize {
                        WINDOW_ENDS.push(self.reference_size);
                    } else {
                        WINDOW_ENDS.push(window_positions[i + 1] - 1);
                    }
                }
            }
        }
    }

    fn set_dup_rate_max_read_starts(&mut self, num: i32) {
        self.read_starts_histogram
            .set_max_read_starts_per_position(num);
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
        (1.0 - self.homopolymer_indels_data[5] as f64
            / (self.num_insertions + self.num_deletions) as f64)
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

    fn get_reads_with_insertion_percentage(&self) -> f64 {
        self.num_reads_with_insertion as f64 / self.number_of_mapped_reads as f64 * 100.0
    }

    fn get_num_deletions(&self) -> usize {
        return self.num_deletions as usize;
    }

    fn get_reads_with_deletions_percentage(&self) -> f64 {
        self.num_reads_with_deletion as f64 / self.number_of_mapped_reads as f64 * 100.0
    }

    fn get_percentage_of_mapped_reads(&self) -> f64 {
        self.number_of_mapped_reads as f64 / self.number_of_reads as f64 * 100.0
    }

    fn get_percentage_of_supp_alignments(&self) -> f64 {
        self.number_of_supp_alignments as f64 / self.number_of_reads as f64 * 100.0
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
        } else if !self.mapping_quality_histogram_map.contains_key(&key) {
            self.mapping_quality_histogram_map.insert(key, 1);
        } else {
            let mut value = self.mapping_quality_histogram_map.get(&key).unwrap();
            let update_value = *value + 1;
            self.mapping_quality_histogram_map.insert(key, update_value);
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
            // update cache
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
        unsafe {
            let chromosome_count = chromosome_window_index.len();
            let contig_records = locator.get_contigs();

            self.chromosome_stats = Vec::with_capacity(chromosome_count);

            for k in 0..chromosome_count {
                let first_window_index = chromosome_window_index[k];
                let last_window_index = if k + 1 < chromosome_count {
                    chromosome_window_index[k + 1] - 1
                } else {
                    NUMBER_OF_TOTAL_WINDOWS as usize - 1
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

        for i in 1..self
            .read_starts_histogram
            .get_max_read_starts_per_position() as usize
            + 1
        {
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
        // let total_size = self.reads_as_data.len()
        //     + self.reads_ts_data.len()
        //     + self.reads_cs_data.len()
        //     + self.reads_gs_data.len()
        //     + self.reads_ns_data.len();

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
    window_inside_reference_size: i64,

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
            window_inside_reference_size: 0,
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

    fn init_window(
        name: String,
        window_start: i64,
        window_end: i64,
        selected_region_available: bool,
        detailed_flag: bool,
    ) -> Self {
        let mut mini_reference: Vec<u8> = vec![];
        unsafe {
            if !CUSTOM_REFERENCES.is_empty() {
                mini_reference = CUSTOM_REFERENCES
                    [(window_start - 1) as usize..(window_end - 1) as usize]
                    .to_vec()
            }
        }
        let mut window = BamGenomeWindow::new(
            name,
            window_start,
            window_end,
            &mini_reference,
            detailed_flag,
        );
        // region analysis
        if selected_region_available {
            Self::calculate_regions_lookup_table_for_window(&mut window);
        }
        window
    }

    fn calculate_regions_lookup_table_for_window(w: &mut BamGenomeWindow) {
        let window_start = w.start;
        let window_end = w.end;
        let window_size = w.window_size as usize;

        let mut bit_set = BitVec::from_elem(window_size, false);
        unsafe {
            let num_regions = SELECTED_REGION_STARTS.len();

            for i in 0..num_regions {
                let region_start = SELECTED_REGION_STARTS[i];
                let region_end = SELECTED_REGION_ENDS[i];

                if region_start == -1 {
                    continue;
                }

                if region_start >= window_start && region_start <= window_end {
                    let end = std::cmp::min(window_end, region_end);
                    Self::set_range_flag(
                        &mut bit_set,
                        (region_start - window_start) as usize,
                        (end - window_start + 1) as usize,
                        true,
                    );
                } else if region_end >= window_start && region_end <= window_end {
                    Self::set_range_flag(
                        &mut bit_set,
                        0,
                        (region_end - window_start + 1) as usize,
                        true,
                    );
                } else if region_start <= window_start && region_end >= window_end {
                    Self::set_range_flag(
                        &mut bit_set,
                        0,
                        (window_end - window_start + 1) as usize,
                        true,
                    );
                }
            }
        }
        w.selected_regions = bit_set;
        w.selected_regions_available_flag = true;
    }

    fn set_range_flag(bit_set: &mut BitVec, start_index: usize, end_index: usize, flag: bool) {
        for index in start_index..end_index {
            bit_set.set(index, flag);
        }
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

    pub fn merge_insert_size(&mut self, correct_insert_size: i32, acum_insert_size: f64) {
        self.correct_insert_sizes += correct_insert_size;
        self.acum_insert_size += acum_insert_size;
    }
    /// Calculate the average metrics
    pub fn compute_descriptors(&mut self) {
        // todo!() add effective_window_length to inside reference size
        self.effective_window_length = self.window_size;
        if self.selected_regions_available_flag {
            self.effective_window_length =
                self.selected_regions.iter().filter(|x| *x).count() as i64
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
        // let mut inverse_selected_regions = BitVec::new();
        // for x in self.selected_regions.iter() {
        //     inverse_selected_regions.push(!x);
        // }
        // self.selected_regions = inverse_selected_regions;
    }

    fn get_region(&self) -> &BitVec {
        &self.selected_regions
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
                first_read_end_pos: read.cigar().end_pos() as i32,
                second_read_start_pos: 0,
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
    pub fn new(_window_start: i64, window_size: i64) -> Self {
        Self {
            number_of_sequenced_bases: 0,
            number_of_mapped_bases: 0,
            number_of_as: 0,
            number_of_ts: 0,
            number_of_cs: 0,
            number_of_gs: 0,
            number_of_ns: 0,
            coverage_data: Vec::with_capacity(window_size as usize),
            mapping_quality_data: Vec::with_capacity(window_size as usize),

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
            'N' => {
                self.number_of_ns += 1;
            }
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
            array.reserve(new_size - old_size);
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

    fn reset_counters(&mut self) {
        self.prev_base = 0;
        self.homopolymer_size = 1;
        self.prev_base_inside_indel_region_flag = false;
        self.homopolymer_starts_inside_indel_region_flag = false;
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
    correct_insert_sizes: i32,
    acum_insert_size: f64,
    number_of_reads_with_start_greater_than_end: usize,
}
impl TaskResult {
    pub fn new() -> Self {
        Self {
            reads_data: None,
            out_region_reads_data: None,
            read_stats_collector: None,
            out_region_read_stats_collector: None,
            correct_insert_sizes: 0,
            acum_insert_size: 0.0,
            number_of_reads_with_start_greater_than_end: 0,
        }
    }

    pub fn acum_insert_size(&mut self, insert_size: i64) {
        if insert_size > 0 {
            self.correct_insert_sizes += 1;
            self.acum_insert_size += insert_size.abs() as f64;
        }
    }

    pub fn get_read_stats_collector(&self) -> &ReadStatsCollector {
        self.read_stats_collector.as_ref().unwrap()
    }

    pub fn get_out_region_read_stats_collector(&self) -> &ReadStatsCollector {
        self.out_region_read_stats_collector.as_ref().unwrap()
    }

    pub fn set_read_stats_collector(&mut self, read_stats_collector: ReadStatsCollector) {
        self.read_stats_collector = Some(read_stats_collector);
    }

    pub fn set_out_of_region_read_stats_collector(
        &mut self,
        read_stats_collector: ReadStatsCollector,
    ) {
        self.out_region_read_stats_collector = Some(read_stats_collector);
    }

    pub fn set_global_reads_data(&mut self, reads_data: Vec<SingleReadData>) {
        self.reads_data = Some(reads_data);
    }

    pub fn set_number_of_reads_with_start_greater_than_end(&mut self, num_reads: usize) {
        self.number_of_reads_with_start_greater_than_end = num_reads;
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

pub struct ProcessBunchOfReadsTask {
    // current_window: BamGenomeWindow,
    compute_insert_size_flag: bool,
    analyze_regions_flag: bool,
    compute_outside_stats_flag: bool,
    analysis_results: HashMap<i64, SingleReadData>,
    out_of_regions_results: HashMap<i64, SingleReadData>,
    read_stats_collector: ReadStatsCollector,
    out_of_regions_read_stats_collector: ReadStatsCollector,
    locator: GenomeLocator,
    number_of_reads_with_start_greater_than_end: usize,

    // window infomation
    current_window_start: i64,
    current_window_end: i64,
    current_window_region: BitVec,
}

impl ProcessBunchOfReadsTask {
    const CIGAR_M: char = 'M';
    const CIGAR_EQ: char = '=';

    pub fn new(
        selected_region_flag: bool,
        compute_outside_stats_flag: bool,
        min_homopolymer_size: i32,
        _locator: GenomeLocator,
        _current_window_start: i64,
        _current_window_end: i64,
        _current_window_region: BitVec,
    ) -> Self {
        Self {
            compute_insert_size_flag: true,
            analyze_regions_flag: selected_region_flag,
            compute_outside_stats_flag: compute_outside_stats_flag,
            analysis_results: HashMap::new(),
            out_of_regions_results: HashMap::new(),
            read_stats_collector: ReadStatsCollector::new(min_homopolymer_size),
            out_of_regions_read_stats_collector: ReadStatsCollector::new(min_homopolymer_size),
            locator: _locator,
            number_of_reads_with_start_greater_than_end: 0,
            current_window_start: _current_window_start,
            current_window_end: _current_window_end,
            current_window_region: _current_window_region,
        }
    }

    fn complement(base: u8) -> u8 {
        match base {
            b'A' => b'T',
            b'T' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            _ => base,
        }
    }
    fn reverse_complete(bases: &mut Vec<u8>) {
        let last_index = bases.len() - 1;

        let mut i = 0;
        let mut j = last_index;

        while i < j {
            let tmp = Self::complement(bases[i]);
            bases[i] = Self::complement(bases[j]);
            bases[j] = tmp;
            i += 1;
            j -= 1;
        }

        if bases.len() % 2 == 1 {
            bases[i] = Self::complement(bases[i]);
        }
    }
    pub fn run(&mut self, reads_bunch: Vec<BamRecord>) -> TaskResult {
        let mut task_result = TaskResult::new();

        for bam_record in &reads_bunch {
            let record = bam_record.record();

            let contig_id = record.tid();
            let contig_option = self.locator.get_contig(contig_id as usize);
            // compute absolute position
            let mut read_absolute_start = -1;
            if let Some(contig) = contig_option {
                read_absolute_start = self
                    .locator
                    .get_absolute_coordinates(contig.name(), record.pos() as i32 + 1);
            }

            let mut alignment: Vec<char> = vec![];
            // compute alignment
            let read_in_region_flag = bam_record.in_region_flag();
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
                // update insert size
                task_result.acum_insert_size(record.insert_size());
            }

            let mapping_quality = record.mapq();
            let read_absolute_end = read_absolute_start + alignment.len() as i64 - 1;

            // acum read
            let region = self.current_window_region.clone();
            let out_of_bounds_flag = self.process_read_alignment(
                self.current_window_end - self.current_window_start + 1,
                self.current_window_start,
                self.current_window_end,
                &region,
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

        let mut reads_data = vec![];

        for val in self.analysis_results.values() {
            reads_data.push(val.clone());
        }

        task_result.set_number_of_reads_with_start_greater_than_end(
            self.number_of_reads_with_start_greater_than_end,
        );
        task_result.set_global_reads_data(reads_data);
        task_result.set_read_stats_collector(self.read_stats_collector.clone());

        if self.analyze_regions_flag && self.compute_outside_stats_flag {
            self.out_of_regions_read_stats_collector.save_gc();

            let mut out_of_region_reads_data = vec![];
            for val in self.out_of_regions_results.values() {
                out_of_region_reads_data.push(val.clone());
            }

            task_result.set_out_of_region_reads_data(out_of_region_reads_data);
            task_result.set_out_of_region_read_stats_collector(
                self.out_of_regions_read_stats_collector.clone(),
            );
        }
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
        let mut read_bases = read.seq().as_bytes();

        // read.getReadNegativeStrandFlag()
        if read.is_reverse() {
            Self::reverse_complete(&mut read_bases);
        }

        // function collect base takes 0.15~0.17 seconds, this cycle costs about 0.4 seconds
        self.read_stats_collector.reset_counters();
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

        for (__, c) in cigar.iter().enumerate() {
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

        let mut read_bases = read.seq().as_bytes();
        self.out_of_regions_read_stats_collector.reset_counters();
        if read.is_reverse() {
            Self::reverse_complete(&mut read_bases);
        }

        // function collect base takes 0.15~0.17 seconds, this cycle costs about 0.4 seconds
        let mut read_pos = 0;
        let mut alignment_pos = 0;
        let mut alignment_vector: Vec<char> = vec![0 as char; alignment_length as usize];

        for pos in 0..extended_cigar_vector.len() {
            let cigar_char = extended_cigar_vector[pos];

            if cigar_char == 'M' || cigar_char == '=' {
                let base = read_bases[read_pos];
                self.out_of_regions_read_stats_collector
                    .collect_base(read_pos, base, false);
                read_pos += 1;
                alignment_vector[alignment_pos as usize] = base as char;
                alignment_pos += 1;
            } else if cigar_char == 'I' {
                let base = read_bases[read_pos];
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

        match read.aux(b"NM") {
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

        if alignment_length < 0 {
            return vec![];
        }
        // precompute total size of alignment
        let mut total_size = 0;
        let cigar = read.cigar();

        for (_i, c) in cigar.iter().enumerate() {
            total_size += c.len();
        }

        // compute extended cigar
        let mut extended_cigar_vector: Vec<char> = vec![' '; total_size as usize];
        let mut start_pos = 0;
        let mut end_pos = 0;
        for (_, c) in cigar.iter().enumerate() {
            end_pos = start_pos + c.len();
            for pos in start_pos..end_pos {
                extended_cigar_vector[pos as usize] = c.char();
            }
            start_pos = end_pos;
        }

        let mut read_pos = 0;
        let read_bases = read.seq();
        let mut alignment_pos = 0;
        let mut alignment_vector: Vec<char> = vec![0 as char; alignment_length as usize];

        for pos in 0..extended_cigar_vector.len() {
            let cigar_char = extended_cigar_vector[pos];

            if cigar_char == 'M' || cigar_char == '=' {
                // get base
                let base = read_bases[read_pos];
                read_pos += 1;
                alignment_vector[alignment_pos as usize] = base as char;
                // set base
                alignment_pos += 1;
            } else if cigar_char == 'I' {
                read_pos += 1;
            } else if cigar_char == 'D' {
                alignment_vector[alignment_pos as usize] = '-';
                alignment_pos += 1;
            } else if cigar_char == 'N' {
                alignment_vector[alignment_pos as usize] = 'N';
                alignment_pos += 1;
            } else if cigar_char == 'S' {
                read_pos += 1;
            } else if cigar_char == 'H' {
            } else if cigar_char == 'P' {
                alignment_vector[alignment_pos as usize] = '-';
                alignment_pos += 1;
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

        let mut in_region_read_data =
            Self::get_window_data(window_start, window_size, &mut self.analysis_results);

        let mut out_region_read_data = None;
        if self.compute_outside_stats_flag {
            out_region_read_data = Some(Self::get_window_data(
                window_start,
                window_size,
                &mut self.out_of_regions_results,
            ));
        }

        let mut read_data = &mut in_region_read_data;

        if read_start > read_end {
            self.number_of_reads_with_start_greater_than_end += 1
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
                        read_data = &mut in_region_read_data;
                    }
                } else if self.compute_outside_stats_flag {
                    read_data = out_region_read_data.as_mut().unwrap();
                } else {
                    continue;
                }
            }

            let nucleotide = alignment[pos as usize];

            // aligned bases
            read_data.number_of_mapped_bases += 1;
            if nucleotide != '-' {
                // mapping quality
                read_data.acum_mapping_quality(relative as u64, mapping_quality as i32);
                // base stats
                read_data.acum_base(relative as u64, nucleotide);
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
        let num_of_processed_windows = unsafe { NUMBER_OF_PROCESSED_WINDOWS };
        let num_of_windows = unsafe { NUMBER_OF_TOTAL_WINDOWS };

        let mut index = num_of_processed_windows + 1;
        let mut out_of_bounds = true;
        let mut adjacent_window: &BamGenomeWindow;
        unsafe {
            while out_of_bounds && index < num_of_windows {
                let window_start = WINDOW_STARTS[index as usize];

                // check if we should insert a new open window
                Self::check_open_windows(window_start, &self.locator, self.analyze_regions_flag);
                let open_windows = OPEN_WINDOWS.lock().unwrap();
                adjacent_window = open_windows.get(&window_start).unwrap();

                if self.compute_outside_stats_flag {
                    // check if we should insert a new open outside window
                    Self::check_open_outside_windows(
                        window_start,
                        &self.locator,
                        self.analyze_regions_flag,
                    );
                    let open_outside_windows = OPEN_OUTSIDE_WINDOWS.lock().unwrap();
                    open_outside_windows.get(&window_start).unwrap();
                }

                // acum read
                let window_size = adjacent_window.get_window_size();
                let window_start = adjacent_window.get_window_start();
                let window_end = adjacent_window.get_window_end();
                let selected_region = adjacent_window.get_region().clone();
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
    }

    fn check_open_windows(window_start: i64, locator: &GenomeLocator, analyze_regions_flag: bool) {
        unsafe {
            let mut open_windows = OPEN_WINDOWS.lock().unwrap();

            if open_windows.contains_key(&window_start) {
                return;
            } else {
                let num_init_windows = NUMBER_OF_INITIALIZED_WINDOWS;
                let window_name = WINDOW_NANMES[num_init_windows as usize].clone();
                let window_end = WINDOW_ENDS[num_init_windows as usize];
                // todo!() it didn't consider custom refernce file
                let reference_size = locator.get_total_size();
                // update number of initialized windows
                NUMBER_OF_INITIALIZED_WINDOWS += 1;
                let new_window = BamGenomeWindow::init_window(
                    window_name,
                    window_start,
                    window_end.min(reference_size),
                    analyze_regions_flag,
                    true,
                );

                open_windows.insert(window_start, new_window);
            }
        }
    }

    fn check_open_outside_windows(
        window_start: i64,
        locator: &GenomeLocator,
        analyze_regions_flag: bool,
    ) {
        unsafe {
            let mut open_outside_windows = OPEN_OUTSIDE_WINDOWS.lock().unwrap();

            if open_outside_windows.contains_key(&window_start) {
                return;
            } else {
                let num_init_windows = NUMBER_OF_INITIALIZED_OUTSIDE_WINDOWS;
                let window_name = OUTSIDE_WINDOW_NANMES[num_init_windows as usize].clone();
                let window_end = WINDOW_ENDS[num_init_windows as usize];
                // todo!() it didn't consider custom refernce file
                let reference_size = locator.get_total_size();
                // update number of initialized windows
                NUMBER_OF_INITIALIZED_OUTSIDE_WINDOWS += 1;
                let new_window = BamGenomeWindow::init_window(
                    window_name,
                    window_start,
                    window_end.min(reference_size),
                    analyze_regions_flag,
                    true,
                );

                open_outside_windows.insert(window_start, new_window);
            }
        }
    }

    fn get_window_data(
        window_start: i64,
        window_size: i64,
        results_map: &mut HashMap<i64, SingleReadData>,
    ) -> &mut SingleReadData {
        if results_map.contains_key(&window_start) {
            results_map.get_mut(&window_start).unwrap()
        } else {
            let data = SingleReadData::new(window_start, window_size);
            results_map.insert(window_start, data);
            results_map.get_mut(&window_start).unwrap()
        }
    }

    fn get_read_stats_collector(&mut self, read: &Record) -> Option<&mut ReadStatsCollector> {
        Some(&mut self.read_stats_collector)
    }
}

struct FinalizeWindowTask<'a> {
    bam_stats: &'a mut BamStats,
    current_window_start: i64,
}

impl<'a> FinalizeWindowTask<'a> {
    pub fn new(_bam_stats: &'a mut BamStats, _current_window_start: i64) -> Self {
        Self {
            bam_stats: _bam_stats,
            current_window_start: _current_window_start,
        }
    }

    pub fn run_inside(&mut self) -> i32 {
        let mut open_windows = OPEN_WINDOWS.lock().unwrap();
        let current_window = open_windows.get_mut(&self.current_window_start).unwrap();
        current_window.compute_descriptors();
        self.bam_stats.add_window_information(current_window);

        // todo!() error handling
        return 0;
    }

    pub fn run_outside(&mut self) -> i32 {
        let mut open_outside_windows = OPEN_OUTSIDE_WINDOWS.lock().unwrap();
        let current_outside_window = open_outside_windows
            .get_mut(&self.current_window_start)
            .unwrap();
        current_outside_window.compute_descriptors();
        self.bam_stats
            .add_window_information(current_outside_window);

        // todo!() error handling
        return 0;
    }
}

pub struct BamStatsAnalysis {
    bam_file: String,

    // reference
    reference_file: String,
    reference_available_flag: bool,
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
    current_window_start: i64,
    current_window_end: i64,
    current_window_region: Option<BitVec>,
    thread_number: usize,
    num_reads_in_bunch: i32,
    progress: usize,
    min_homopolymer_size: i32,
    dup_rate_max_read_starts: i32,

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
    current_outside_window_start: i64,
    current_outside_window_end: i64,
    current_outside_window_region: Option<BitVec>,
    open_outside_windows: HashMap<i64, BamGenomeWindow>,
    outside_bam_stats: Option<BamStats>,

    // read size
    acum_read_size: usize,
    max_read_size: usize,
    min_read_size: usize,

    // region
    selected_region_starts: Vec<i64>,
    selected_region_ends: Vec<i64>,
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
    results: Vec<Future<TaskResult>>,
    // Future<Integer> finalizeWindowResult;
    // #[serde(skip_serializing)]
    // #[serde(skip_deserializing)]
    finalize_window_result: i32,
    time_to_calc_overlappers: usize,
    pg_program: String,
    pg_command_string: String,
    // // multi thread
    thread_pool: ExecutorService<Box<dyn FnOnce() -> TaskResult + Send + 'static>, TaskResult>,
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

    pub fn new(_bam_file: String, qualimap_config: &QualimapConfig) -> Self {
        let thread_num = qualimap_config.get_thread_num();
        let window_num = qualimap_config.get_window_num();
        let bun_size = qualimap_config.get_bunch_size();
        let min_homopolymer_size = qualimap_config.get_min_homopolymer_size();
        let dup_rate_max_read_starts = qualimap_config.get_dup_rate_max_read_starts();
        let _feature_file = qualimap_config.get_feature_file();

        let outside_analyze_flag = qualimap_config.get_outside_region_analyze_flag();
        let lib_protocal = qualimap_config.get_lib_protocol();

        let overlap_analyze_flag = qualimap_config.get_collect_overlap_flag();

        let skip_dup_flag = qualimap_config.get_skip_duplicate_flag();
        let skip_dup_mode = qualimap_config.get_skip_duplicate_mode();
        let mut _skip_marked_duplicates_flag = false;
        let mut _skip_detected_duplicates_flag = false;

        if skip_dup_flag {
            if skip_dup_mode == "flagged" {
                _skip_detected_duplicates_flag = true;
            } else if skip_dup_mode == "estimated" {
                _skip_marked_duplicates_flag = true;
            } else if skip_dup_mode == "both" {
                _skip_detected_duplicates_flag = true;
                _skip_marked_duplicates_flag = true;
            }
        }

        let mut _selected_regions_available_flag = false;
        if !_feature_file.is_empty() {
            _selected_regions_available_flag = true;
        }
        Self {
            bam_file: _bam_file,
            // reference
            reference_file: "".to_string(),
            reference_available_flag: false,
            reference_size: 0,
            number_of_reference_config: 0,

            // currentWindow management
            number_of_windows: window_num as i32,
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
            current_window_start: -1,
            current_window_end: -1,
            current_window_region: None,
            thread_number: thread_num,
            num_reads_in_bunch: bun_size as i32,
            progress: 0,
            min_homopolymer_size: min_homopolymer_size as i32,
            dup_rate_max_read_starts: dup_rate_max_read_starts as i32,

            // nucleotide reporting
            out_dir: ".".to_string(),

            // gff support
            selected_regions_available_flag: _selected_regions_available_flag,
            feature_file: _feature_file,
            num_of_selected_regions: 0,

            // inside region analysis
            inside_reference_sizes: 0,
            skip_marked_duplicates_flag: _skip_marked_duplicates_flag,
            skip_detected_duplicates_flag: _skip_detected_duplicates_flag,
            collect_intersecting_paired_end_reads_flag: overlap_analyze_flag,

            // outside region analysis
            compute_outside_stats_flag: outside_analyze_flag,
            current_outside_window_start: -1,
            current_outside_window_end: -1,
            current_outside_window_region: None,
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
            protocol: LibraryProtocol::get_protocol_by_name(lib_protocal),

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
            thread_pool: Executors::new_fixed_thread_pool(thread_num as u32)
                .expect("Failed to create the thread pool"),
        }
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
        self.load_locator(&hash_header);
        self.load_program_records(&hash_header);

        // load reference
        last_action_done = "Loading reference...";
        println!("{}", last_action_done);
        self.load_reference();

        // init window set
        self.window_size = self.compute_window_size(self.reference_size, self.number_of_windows);
        let window_positions: Vec<i64> = self.compute_window_positions(self.window_size);
        self.effective_number_of_window = window_positions.len() as i32;
        unsafe {
            NUMBER_OF_TOTAL_WINDOWS = self.effective_number_of_window;
        }

        if self.effective_number_of_window > QualimapConstants::DEFAULT_STABLIZED_WINDOW_PROPORTION
        {
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
        println!(
            "Current number of threads: {}, and max number support in your device is: {}",
            self.thread_number,
            num_cpus::get()
        );
        // initialize BamStats
        let mut bam_stats_instance = BamStats::new(
            "genome".to_string(),
            self.locator.clone(),
            self.reference_size,
            self.effective_number_of_window as usize,
        );

        bam_stats_instance.set_source_file(self.bam_file.clone());
        bam_stats_instance.set_window_references("w", window_positions.clone());
        bam_stats_instance.set_dup_rate_max_read_starts(self.dup_rate_max_read_starts);
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
            // load selected regions
            self.load_selected_regions();

            // outside of regions stats
            if self.compute_outside_stats_flag {
                // initialize outside bamStats
                let mut outside_bam_stats_instance = BamStats::new(
                    "outside".to_string(),
                    self.locator.clone(),
                    self.reference_size,
                    self.effective_number_of_window as usize,
                );

                outside_bam_stats_instance.set_source_file(self.bam_file.clone());
                outside_bam_stats_instance.set_outside_window_references("out_w", window_positions);
                outside_bam_stats_instance
                    .set_dup_rate_max_read_starts(self.dup_rate_max_read_starts);
                self.outside_bam_stats = Some(outside_bam_stats_instance);

                // get next outside window
                self.next_outside_window(true);

                // todo()!
                // if(activeReporting) {
                //     outsideBamStats.activateWindowReporting(outdir + "/outside_window.txt");
                // }

                // if(saveCoverage){
                //     outsideBamStats.activateCoverageReporting(outdir + "/outside_coverage.txt", nonZeroCoverageOnly);
                // }

                // we have twice more data from the bunch, so the queue is limited now
                self.max_size_of_task_queue /= 2;
            }
        }

        // get next window
        self.next_window(true);
    }

    fn load_bed_regions(&mut self) -> io::Result<()> {
        let mut regions_with_missing_chromosomes_count = 0;

        let mut reader = File::open(&self.feature_file)
            .map(BufReader::new)
            .map(bed::Reader::new)?;

        println!("Initializing regions from {}.....", &self.feature_file);
        for result in reader.records::<3>() {
            let _ = result?;
            self.num_of_selected_regions += 1;
        }
        if self.num_of_selected_regions == 0 {
            panic!("Failed to load selected regions.");
        }
        println!("Found {} regions", self.num_of_selected_regions);

        self.selected_region_starts = vec![0; self.num_of_selected_regions];
        self.selected_region_ends = vec![0; self.num_of_selected_regions];

        println!("Filling region references... ");
        reader = File::open(&self.feature_file)
            .map(BufReader::new)
            .map(bed::Reader::new)?;
        let mut pos = -1;
        let mut index = 0;

        for result in reader.records::<3>() {
            let record = result?;
            pos = self.locator.get_absolute_coordinates(
                record.reference_sequence_name(),
                record.start_position().get() as i32,
            );
            if pos == -1 {
                self.selected_region_starts[index] = -1;
                self.selected_region_ends[index] = -1;
                regions_with_missing_chromosomes_count += 1;
                continue;
            }

            let region_length = record.end_position().get() - record.start_position().get() + 1;
            self.selected_region_starts[index] = pos;
            self.selected_region_ends[index] = pos + region_length as i64 - 1;

            let record_string = record.to_string();
            let fields: Vec<&str> = record_string.split("\t").collect();
            let mut is_negative_strand = false;

            // NOTE: strand-specificity is analyzed in BAM file, but not for 3-column BED
            if fields.len() >= 6 && fields[5] == "-" {
                is_negative_strand = true;
            }
            self.region_overlap_lookup_table.put_region(
                record.start_position().get(),
                record.end_position().get(),
                record.reference_sequence_name(),
                !is_negative_strand,
            );
            index += 1;
        }

        if regions_with_missing_chromosomes_count == self.num_of_selected_regions {
            panic!("The feature file with regions can not be associated with the BAM file.\n
            Please check, if the chromosome names match in the annotation file and the alignment file.");
        } else if regions_with_missing_chromosomes_count > 0 {
            println!(
                "{} regions were skipped because chromosome name was not found in the BAM file.",
                regions_with_missing_chromosomes_count
            );
        }
        Ok(())
    }

    fn load_gff_regions(&mut self) -> io::Result<()> {
        let mut regions_with_missing_chromosomes_count = 0;

        let mut reader = File::open(&self.feature_file)
            .map(BufReader::new)
            .map(gff::Reader::new)?;
        println!("Initializing regions from {}.....", &self.feature_file);
        for result in reader.records() {
            let _ = result?;
            self.num_of_selected_regions += 1;
        }
        if self.num_of_selected_regions == 0 {
            panic!("Failed to load selected regions.");
        }
        println!("Found {} regions", self.num_of_selected_regions);

        self.selected_region_starts = vec![0; self.num_of_selected_regions];
        self.selected_region_ends = vec![0; self.num_of_selected_regions];

        println!("Filling region references... ");
        reader = File::open(&self.feature_file)
            .map(BufReader::new)
            .map(gff::Reader::new)?;
        let mut pos = -1;
        let mut index = 0;

        for result in reader.records() {
            let record = result?;
            pos = self.locator.get_absolute_coordinates(
                record.reference_sequence_name(),
                record.start().get() as i32,
            );
            if pos == -1 {
                self.selected_region_starts[index] = -1;
                self.selected_region_ends[index] = -1;
                regions_with_missing_chromosomes_count += 1;
                continue;
            }

            let region_length = record.end().get() - record.start().get() + 1;
            self.selected_region_starts[index] = pos;
            self.selected_region_ends[index] = pos + region_length as i64 - 1;

            let mut is_negative_strand = false;

            if record.strand().as_ref() == "-" {
                is_negative_strand = true;
            }
            self.region_overlap_lookup_table.put_region(
                record.start().get(),
                record.end().get(),
                record.reference_sequence_name(),
                !is_negative_strand,
            );
            index += 1;
        }
        if regions_with_missing_chromosomes_count == self.num_of_selected_regions {
            panic!("The feature file with regions can not be associated with the BAM file.\n
            Please check, if the chromosome names match in the annotation file and the alignment file.");
        } else if regions_with_missing_chromosomes_count > 0 {
            println!(
                "{} regions were skipped because chromosome name was not found in the BAM file.",
                regions_with_missing_chromosomes_count
            );
        }
        Ok(())
    }

    fn load_gtf_regions(&mut self) -> io::Result<()> {
        let mut regions_with_missing_chromosomes_count = 0;

        let mut reader = File::open(&self.feature_file)
            .map(BufReader::new)
            .map(gtf::Reader::new)?;
        println!("Initializing regions from {}.....", &self.feature_file);
        for result in reader.records() {
            let _ = result?;
            self.num_of_selected_regions += 1;
        }
        if self.num_of_selected_regions == 0 {
            panic!("Failed to load selected regions.");
        }
        println!("Found {} regions", self.num_of_selected_regions);

        self.selected_region_starts = vec![0; self.num_of_selected_regions];
        self.selected_region_ends = vec![0; self.num_of_selected_regions];

        println!("Filling region references... ");
        reader = File::open(&self.feature_file)
            .map(BufReader::new)
            .map(gtf::Reader::new)?;
        let mut pos = -1;
        let mut index = 0;

        for result in reader.records() {
            let record = result?;
            pos = self.locator.get_absolute_coordinates(
                record.reference_sequence_name(),
                record.start().get() as i32,
            );
            if pos == -1 {
                self.selected_region_starts[index] = -1;
                self.selected_region_ends[index] = -1;
                regions_with_missing_chromosomes_count += 1;
                continue;
            }

            let region_length = record.end().get() - record.start().get() + 1;
            self.selected_region_starts[index] = pos;
            self.selected_region_ends[index] = pos + region_length as i64 - 1;

            let mut is_negative_strand = false;

            if record.strand().is_some() && record.strand().unwrap().as_ref() == "-" {
                is_negative_strand = true;
            }
            self.region_overlap_lookup_table.put_region(
                record.start().get(),
                record.end().get(),
                record.reference_sequence_name(),
                !is_negative_strand,
            );
            index += 1;
        }
        if regions_with_missing_chromosomes_count == self.num_of_selected_regions {
            panic!("The feature file with regions can not be associated with the BAM file.\n
                Please check, if the chromosome names match in the annotation file and the alignment file.");
        } else if regions_with_missing_chromosomes_count > 0 {
            println!(
                "{} regions were skipped because chromosome name was not found in the BAM file.",
                regions_with_missing_chromosomes_count
            );
        }
        Ok(())
    }

    fn load_selected_regions(&mut self) -> io::Result<()> {
        let gff_re = Regex::new(r"\.gff$").unwrap();
        let gtf_re = Regex::new(r"\.gtf$").unwrap();
        let bed_re = Regex::new(r"\.bed$").unwrap();

        if bed_re.is_match(&self.feature_file) {
            self.load_bed_regions();
        } else if gff_re.is_match(&self.feature_file) {
            self.load_gff_regions();
        } else if gtf_re.is_match(&self.feature_file) {
            self.load_gtf_regions();
        } else {
            panic!("Unknown feature file format. Please provide file in GFF/GTF or BED format.")
        }

        unsafe {
            SELECTED_REGION_STARTS = self.selected_region_starts.clone();
            SELECTED_REGION_ENDS = self.selected_region_ends.clone();
        }

        Ok(())
    }

    pub fn run(&mut self) {
        let start_time = Instant::now();

        self.load_and_init();
        println!("Time to init: {:?}", Instant::now() - start_time);

        let mut bam_reader = Reader::from_path(&self.bam_file).unwrap();
        let thread_pool = tpool::ThreadPool::new(num_cpus::get() as u32).unwrap();
        bam_reader.set_thread_pool(&thread_pool);

        let mut record = Record::new();
        let mut process_time = Duration::seconds(0);
        let mut time_to_finalize_and_get_next_window = Duration::seconds(0);
        let mut time_to_analyze_reads_bunch = Duration::seconds(0);
        let mut time_to_finish_reads_bunch = Duration::seconds(0);
        let mut time_to_task_run = Duration::seconds(0);
        let mut inner_circle = Duration::seconds(0);
        let mut reads_bunch: Vec<BamRecord> = Vec::with_capacity(self.num_reads_in_bunch as usize);
        while let Some(result) = bam_reader.read(&mut record) {
            let inner_circle_start = Instant::now();
            match result {
                Ok(_) => {
                    let start_time = Instant::now();
                    self.process_sequence(
                        &record,
                        &mut reads_bunch,
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
            inner_circle += Instant::now() - inner_circle_start;
        }
        println!(
            "Time to read from bam files: {:?}",
            Instant::now() - start_time - inner_circle
        );

        let last_process_start = Instant::now();
        let number_of_processed_windows = unsafe { NUMBER_OF_PROCESSED_WINDOWS };
        if reads_bunch.is_empty() || number_of_processed_windows < self.effective_number_of_window {
            let num_windows = self.effective_number_of_window;
            let last_position = unsafe { WINDOW_ENDS[num_windows as usize - 1] + 1 };
            self.collect_analysis_results(
                &mut reads_bunch,
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
                // add lock
                let mut open_outside_windows = OPEN_OUTSIDE_WINDOWS.lock().unwrap();
                let current_outside_window = open_outside_windows
                    .get_mut(&self.current_outside_window_start)
                    .unwrap();
                current_outside_window.inverse_regions();
                // release lock
                drop(open_outside_windows);

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

        let bam_stats_final = self.bam_stats.as_mut().unwrap();
        println!("Total processed windows:{}", unsafe {
            NUMBER_OF_PROCESSED_WINDOWS
        });
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
            println!("WARNING: make sure duplicate alignments are flagged in the BAM file or apply a different skip duplicates mode.")
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
            total_number_of_paired_reads = total_number_of_paired_reads
                + self.outside_bam_stats_collector.get_num_paired_reads();
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
            println!(
                "Total number of mapped reads or mapped reads in region equals zero.\n
                        For more details, check the number of Unmapped reads."
            );
        }

        if self.selected_regions_available_flag && self.compute_outside_stats_flag {
            let outside_bam_stats_final = self.outside_bam_stats.as_mut().unwrap();
            outside_bam_stats_final.reference_size = self.reference_size;
            outside_bam_stats_final.number_of_reference_contigs =
                self.locator.get_contigs().len() as i64;
            outside_bam_stats_final.num_selected_regions = self.num_of_selected_regions as i32;
            outside_bam_stats_final.in_region_reference_size =
                self.reference_size as i64 - self.inside_reference_sizes as i64;

            outside_bam_stats_final.number_of_reads = self.number_of_reads as i64;
            outside_bam_stats_final.number_of_secondary_alignments =
                self.number_of_secondary_alignments as i64;

            outside_bam_stats_final.number_of_mapped_reads = total_number_of_mapped_reads as i64;
            outside_bam_stats_final.number_of_paired_reads = total_number_of_paired_reads as i64;
            outside_bam_stats_final.number_of_mapped_first_of_pair =
                total_number_of_mapped_first_of_pair as i64;
            outside_bam_stats_final.number_of_mapped_second_of_pair =
                total_number_of_mapped_second_of_pair as i64;
            outside_bam_stats_final.number_of_singletons = total_number_of_singletons as i64;

            // set bam_stats with metrics inside of region
            outside_bam_stats_final.number_of_mapped_reads_in_regions =
                self.outside_bam_stats_collector.get_num_mapped_reads() as i64;
            outside_bam_stats_final.number_of_paired_reads_in_regions =
                self.outside_bam_stats_collector.get_num_paired_reads() as i64;
            outside_bam_stats_final.number_of_mapped_first_of_pair_in_regions =
                self.outside_bam_stats_collector
                    .get_num_mapped_first_in_pair() as i64;
            outside_bam_stats_final.number_of_mapped_second_of_pair_in_regions =
                self.outside_bam_stats_collector
                    .get_num_mapped_second_in_pair() as i64;
            outside_bam_stats_final.number_of_singletons_in_regions =
                self.outside_bam_stats_collector.get_num_singletons() as i64;

            outside_bam_stats_final.read_max_size = self.max_read_size as i32;
            outside_bam_stats_final.read_min_size = self.min_read_size as i32;
            outside_bam_stats_final.read_mean_size =
                self.acum_read_size as f64 / self.number_of_reads as f64;

            outside_bam_stats_final.num_detected_duplicate_reads =
                self.outside_bam_stats_collector.num_marked_duplicates;

            if self.outside_bam_stats_collector.get_num_mapped_reads() > 0 {
                println!("Computing descriptors for outside regions...");
                outside_bam_stats_final.compute_descriptors();
                println!("Computing per chromosome statistics for outside regions...");
                outside_bam_stats_final
                    .compute_chromosome_stats(&self.locator, &self.chromosome_window_indexes);
                println!("Computing histograms for outside regions...");
                outside_bam_stats_final.compute_histograms();
            } else {
                println!("\nWARNING: number of mapped reads outside of regions equals zero");
                println!(
                    "Total number of mapped reads or mapped reads outside of  region equals zero.\n
                                For more details, check the number of Unmapped reads."
                );
            }
        }

        let overall_time = Instant::now();
        println!("Overall analysis time: {:?}", overall_time - start_time);
    }

    pub fn process_sequence(
        &mut self,
        record: &Record,
        reads_bunch: &mut Vec<BamRecord>,
        time_to_finalize_and_get_next_window: &mut Duration,
        time_to_analyze_reads_bunch: &mut Duration,
        time_to_finish_reads_bunch: &mut Duration,
        time_to_task_run: &mut Duration,
    ) {
        let contig_id = record.tid();
        let contig_option = self.locator.get_contig(contig_id as usize);
        // compute absolute position

        let mut absolute_position = -1;
        if let Some(contig) = contig_option {
            absolute_position = self
                .locator
                .get_absolute_coordinates(contig.name(), record.pos() as i32 + 1);
        }

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

        let novel_read_flag =
            (record.flags() & QualimapConstants::SAM_FLAG_SUPP_ALIGNMENT as u16) == 0;
        if novel_read_flag {
            self.number_of_reads += 1;
        }

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

            if absolute_position < self.current_window_start {
                panic!("The alignment file is unsorted.\nPlease sort the BAM file by coordinate.");
            }

            let mut bam_record = BamRecord::new(record.to_owned());

            if self.selected_regions_available_flag {
                let read_overlaps_regions_flag = self.judge_read_overlap_regions(record);
                if read_overlaps_regions_flag {
                    // set in region flag
                    bam_record.set_in_region_flag(true);

                    let bam_stats = self.bam_stats.as_mut().unwrap();
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
                        self.bam_stats_collector.collect_paired_read_info(record);
                    }

                    bam_stats.update_insert_size_histogram(insert_size as i32);
                } else {
                    // set out of region flag
                    bam_record.set_in_region_flag(false);

                    let outside_bam_stats = self.outside_bam_stats.as_mut().unwrap();
                    if self
                        .outside_bam_stats_collector
                        .update_stats_and_judge_is_dup(record)
                        && self.skip_marked_duplicates_flag
                    {
                        self.number_of_duplicates_skip += 1;
                        return;
                    }
                    if self.compute_outside_stats_flag {
                        if outside_bam_stats.update_read_start_histogram_and_judge_is_detected_dup(
                            absolute_position,
                        ) && self.skip_detected_duplicates_flag
                        {
                            self.number_of_duplicates_skip += 1;
                            return;
                        }
                        outside_bam_stats.update_insert_size_histogram(insert_size as i32);
                    }
                }
            } else {
                // set in region flag
                bam_record.set_in_region_flag(true);

                let bam_stats = self.bam_stats.as_mut().unwrap();
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
                    self.bam_stats_collector.collect_paired_read_info(record);
                }
                bam_stats.update_insert_size_histogram(insert_size as i32);
            }

            // finalize current and get next window
            if absolute_position > self.current_window_end {
                //analyzeReads(readsBunch);
                self.collect_analysis_results(
                    reads_bunch,
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
                    // add lock
                    let mut open_outside_windows = OPEN_OUTSIDE_WINDOWS.lock().unwrap();
                    let current_outside_window = open_outside_windows
                        .get_mut(&self.current_outside_window_start)
                        .unwrap();
                    current_outside_window.inverse_regions();
                    // release lock
                    drop(open_outside_windows);

                    self.finalize_and_get_next_outside_window(absolute_position, true);
                }
            }

            if self.current_window_start == -1 {
                //Some reads are out of reference bounds?
                return;
            }
            reads_bunch.push(bam_record);

            if reads_bunch.len() >= self.num_reads_in_bunch as usize {
                if self.results.len() >= self.max_size_of_task_queue {
                    self.collect_analysis_results(
                        reads_bunch,
                        time_to_analyze_reads_bunch,
                        time_to_finish_reads_bunch,
                        time_to_task_run,
                    );
                } else {
                    self.analyze_reads_bunch(
                        reads_bunch,
                        time_to_task_run,
                        time_to_analyze_reads_bunch,
                    );
                }
            }

            self.number_of_valid_reads += 1;
        } else {
            self.number_of_problematic_reads += 1;
        }
    }

    fn judge_read_overlap_regions(&mut self, read: &Record) -> bool {
        let contig_id = read.tid();
        let contig = self.locator.get_contig(contig_id as usize).unwrap();

        if self.protocol == LibraryProtocol::NonStrandSpecific {
            return self.region_overlap_lookup_table.overlaps(
                read.pos() as usize + 1,
                read.cigar().end_pos() as usize,
                contig.name(),
            );
        } else {
            let read_has_forward_strand = !read.is_reverse();
            let forward_transcript_strand_is_expected = match (read.is_paired(), &self.protocol) {
                (true, LibraryProtocol::StrandSpecificForward) => {
                    (read.is_first_in_template() && read_has_forward_strand)
                        || (read.is_last_in_template() && !read_has_forward_strand)
                }
                (true, LibraryProtocol::StrandSpecificReverse) => {
                    (read.is_first_in_template() && !read_has_forward_strand)
                        || (read.is_last_in_template() && read_has_forward_strand)
                }
                _ => {
                    self.protocol == LibraryProtocol::StrandSpecificForward
                        && read_has_forward_strand
                }
            };

            let overlap_result = self.region_overlap_lookup_table.overlaps_with_strand(
                read.pos() as usize + 1,
                read.cigar().end_pos() as usize,
                contig.name(),
                forward_transcript_strand_is_expected,
            );

            if overlap_result.is_strand_matches() {
                self.number_of_correct_strand_reads += 1;
            }

            return overlap_result.is_interval_overlap();
        }
    }

    fn finalize_and_get_next_window(
        &mut self,
        position: i64,
        detailed_flag: bool,
        time_to_finalize_and_get_next_window: &mut Duration,
    ) {
        let mut current_window_end = self.current_window_end;
        while position > current_window_end {
            let finalize_and_get_next_start = Instant::now();
            // bamstats and current window has been updated in this function
            self.finalize_window();
            *time_to_finalize_and_get_next_window += Instant::now() - finalize_and_get_next_start;
            // update current window
            self.next_window(detailed_flag);

            if self.current_window_start == -1 {
                break;
            }
            current_window_end = self.current_window_end;
        }
    }

    fn finalize_and_get_next_outside_window(&mut self, position: i64, detailed_flag: bool) {
        let mut current_outside_window_end = self.current_outside_window_end;
        while position > current_outside_window_end {
            // bamstats and current window has been updated in this function
            self.finalize_outside_window();
            // update current outside window
            self.next_outside_window(detailed_flag);

            if self.current_outside_window_start == -1 {
                break;
            }
            current_outside_window_end = self.current_outside_window_end;
        }
    }

    fn finalize_window(&mut self) {
        let bam_stats = self.bam_stats.as_mut().unwrap();

        let mut finalize_window_task =
            FinalizeWindowTask::new(bam_stats, self.current_window_start);
        self.finalize_window_result = finalize_window_task.run_inside();

        unsafe {
            NUMBER_OF_PROCESSED_WINDOWS += 1;

            // report progress
            let num_processed_windows = NUMBER_OF_PROCESSED_WINDOWS;
            if num_processed_windows as usize % self.window_block_size_to_report == 0 {
                println!(
                    "Processed {} out of {} windows...",
                    { num_processed_windows },
                    { self.effective_number_of_window }
                );
            }
        }
        self.time_to_calc_overlappers = 0;

        // remove the finished window
        OPEN_WINDOWS
            .lock()
            .unwrap()
            .remove(&self.current_window_start);
    }

    fn finalize_outside_window(&mut self) {
        let outside_bam_stats = self.outside_bam_stats.as_mut().unwrap();

        let mut finalize_window_task =
            FinalizeWindowTask::new(outside_bam_stats, self.current_outside_window_start);
        self.finalize_window_result = finalize_window_task.run_outside();

        unsafe {
            NUMBER_OF_PROCESSED_OUTSIDE_WINDOWS += 1;
        }

        // remove the finished window
        OPEN_OUTSIDE_WINDOWS
            .lock()
            .unwrap()
            .remove(&self.current_outside_window_start);
    }

    fn collect_analysis_results(
        &mut self,
        reads_bunch: &mut Vec<BamRecord>,
        time_to_analyze_reads_bunch: &mut Duration,
        time_to_finish_reads_bunch: &mut Duration,
        time_to_task_run: &mut Duration,
    ) {
        // start last bunch
        self.analyze_reads_bunch(reads_bunch, time_to_task_run, time_to_analyze_reads_bunch);

        let finsih_task_start = Instant::now();

        let mut number_of_reads_with_start_greater_than_end = 0;
        // wait till all tasks are finished
        for future_result in &self.results {
            let task_result = future_result.get().expect("Couldn't get a result");
            // todo!() merge insert size may have a bug caused by boundary
            let mut open_windows = OPEN_WINDOWS.lock().unwrap();
            let mut open_outside_windows = OPEN_OUTSIDE_WINDOWS.lock().unwrap();
            let current_window = open_windows.get_mut(&self.current_window_start).unwrap();
            current_window.merge_insert_size(
                task_result.correct_insert_sizes,
                task_result.acum_insert_size,
            );
            // acum reads with start greater than end in a task
            number_of_reads_with_start_greater_than_end +=
                task_result.number_of_reads_with_start_greater_than_end;

            let data_sets = task_result.get_read_alignment_data();
            // self.merge_number_of_reads_with_start_greater_then_end()
            for single_read_data in data_sets {
                let window_start = single_read_data.get_window_start();
                let window = open_windows.get_mut(&window_start).unwrap();
                window.add_read_alignment_data(single_read_data);
            }

            // update bam stats infomation from single read data
            self.bam_stats
                .as_mut()
                .unwrap()
                .add_read_stats_data(task_result.get_read_stats_collector());

            if self.selected_regions_available_flag && self.compute_outside_stats_flag {
                let outside_data_sets = task_result.get_out_of_region_reads_data();
                // update window infomation from single read data
                for single_read_data in outside_data_sets {
                    let window_start = single_read_data.get_window_start();
                    let window = open_outside_windows.get_mut(&window_start).unwrap();
                    window.add_read_alignment_data(single_read_data);
                }
                // update bam stats infomation from single read data
                self.outside_bam_stats
                    .as_mut()
                    .unwrap()
                    .add_read_stats_data(task_result.get_out_region_read_stats_collector());
            }
        }
        // update number_of_reads_with_start_greater_than_end
        self.merge_number_of_reads_with_start_greater_then_end(
            number_of_reads_with_start_greater_than_end,
        );

        self.results.clear();
        *time_to_finish_reads_bunch += Instant::now() - finsih_task_start;
    }

    fn analyze_reads_bunch(
        &mut self,
        reads_bunch: &mut Vec<BamRecord>,
        time_to_task_run: &mut Duration,
        time_to_analyze_reads_bunch: &mut Duration,
    ) {
        let last_bunch_start = Instant::now();
        let locator = self.locator.clone();
        let selected_region_flag = self.get_selected_regions_available_flag();
        let compute_outside_stats_flag = self.get_compute_outside_stats_flag();
        let min_homopolymer_size = self.get_min_homopolymer_size();

        let current_window_start = self.current_window_start;
        let current_window_end = self.current_window_end;
        let current_window_region = self.current_window_region.as_ref().unwrap().clone();

        let reads_bunch_copy = reads_bunch.clone();
        let mut task = ProcessBunchOfReadsTask::new(
            selected_region_flag,
            compute_outside_stats_flag,
            min_homopolymer_size,
            locator,
            current_window_start,
            current_window_end,
            current_window_region,
        );
        let future_result = self
            .thread_pool
            .submit_async(Box::new(move || task.run(reads_bunch_copy)))
            .expect("Failed to submit function");

        self.results.push(future_result);

        reads_bunch.clear();
        *time_to_analyze_reads_bunch += Instant::now() - last_bunch_start;
    }

    fn load_locator(&mut self, hash_header: &HashMap<String, Vec<LinearMap<String, String>>>) {
        if hash_header.contains_key("SQ") {
            let sq_records = &hash_header["SQ"];
            for record in sq_records {
                let contig_name = record.get("SN").unwrap();
                let contig_length = record.get("LN").unwrap().parse::<i32>().unwrap();
                self.locator.add_contig(contig_name.clone(), contig_length);
            }
        }
    }

    fn load_program_records(
        &mut self,
        hash_header: &HashMap<String, Vec<LinearMap<String, String>>>,
    ) {
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
        let mut window_starts: Vec<i64> = vec![];

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

    fn next_window(&mut self, detailed_flag: bool) {
        // reset current window
        self.current_window_start = -1;
        self.current_window_end = -1;

        let bam_stats = self.bam_stats.as_mut().unwrap();

        let mut open_windows = OPEN_WINDOWS.lock().unwrap();

        unsafe {
            let num_processed_windows = NUMBER_OF_PROCESSED_WINDOWS;
            let num_total_windows = NUMBER_OF_TOTAL_WINDOWS;
            if num_processed_windows < num_total_windows {
                let window_start = WINDOW_STARTS[num_processed_windows as usize];
                let window_end = WINDOW_ENDS[num_processed_windows as usize];
                let reference_size = bam_stats.get_reference_size();
                let current_window_name = WINDOW_NANMES[num_processed_windows as usize].clone();

                if window_start <= reference_size {
                    if open_windows.contains_key(&window_start) {
                        let window = open_windows.get(&window_start).unwrap();
                        self.current_window_start = window_start;
                        self.current_window_end = window_end;
                        self.current_window_region = Some(window.get_region().clone());
                    } else {
                        NUMBER_OF_INITIALIZED_WINDOWS += 1;

                        let current_window = BamGenomeWindow::init_window(
                            current_window_name,
                            window_start,
                            window_end.min(reference_size),
                            self.selected_regions_available_flag,
                            detailed_flag,
                        );
                        self.current_window_region = Some(current_window.get_region().clone());
                        open_windows.insert(window_start, current_window);
                        self.current_window_start = window_start;
                        self.current_window_end = window_end;
                    }
                }
            }
        }
    }

    fn next_outside_window(&mut self, detailed_flag: bool) {
        self.current_outside_window_start = -1;
        self.current_outside_window_end = -1;

        let bam_stats = self.outside_bam_stats.as_mut().unwrap();

        let mut open_outside_windows = OPEN_OUTSIDE_WINDOWS.lock().unwrap();

        unsafe {
            let num_processed_windows = NUMBER_OF_PROCESSED_OUTSIDE_WINDOWS;
            let num_total_windows = NUMBER_OF_TOTAL_WINDOWS;
            if num_processed_windows < num_total_windows {
                let window_start = WINDOW_STARTS[num_processed_windows as usize];
                let window_end = WINDOW_ENDS[num_processed_windows as usize];
                let reference_size = bam_stats.get_reference_size();
                let current_window_name =
                    OUTSIDE_WINDOW_NANMES[num_processed_windows as usize].clone();

                if window_start <= reference_size {
                    if open_outside_windows.contains_key(&window_start) {
                        let outside_window = open_outside_windows.get(&window_start).unwrap();
                        self.current_outside_window_start = window_start;
                        self.current_outside_window_end = window_end;
                        self.current_outside_window_region =
                            Some(outside_window.get_region().clone());
                    } else {
                        NUMBER_OF_INITIALIZED_OUTSIDE_WINDOWS += 1;

                        let current_outside_window = BamGenomeWindow::init_window(
                            current_window_name,
                            window_start,
                            window_end.min(reference_size),
                            self.selected_regions_available_flag,
                            detailed_flag,
                        );
                        self.current_outside_window_region =
                            Some(current_outside_window.get_region().clone());
                        open_outside_windows.insert(window_start, current_outside_window);
                        self.current_outside_window_start = window_start;
                        self.current_outside_window_end = window_end;
                    }
                }
            }
        }
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

    fn merge_number_of_reads_with_start_greater_then_end(&mut self, num_reads: usize) {
        self.number_of_reads_with_start_greater_than_end += num_reads;
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
pub struct DetailedStats {
    // input
    bam_file: String,

    // reference
    num_of_bases: usize,
    num_of_contigs: usize,

    // globals
    num_of_windows: usize,
    num_of_reads: usize,
    num_of_mapped_reads: usize,
    percent_of_mapped_reads: f64,
    num_of_supplementary_alignments: usize,
    percent_of_supplementary_alignments: f64,
    num_of_second_alignments: usize,
    num_of_paired_reads: usize,
    num_of_mapped_reads_first_in_pair: usize,
    num_of_mapped_reads_second_in_pair: usize,
    num_of_mapped_reads_both_in_pair: usize,
    num_of_mapped_reads_singletons: usize,
    num_of_overlapping_reads: usize,
    num_of_mapped_bases: usize,
    num_of_sequenced_bases: usize,
    num_of_aligned_bases: usize,

    num_of_flagged_duplicate_reads: usize,
    num_of_estimated_duplicate_reads: usize,
    duplication_rate: f64,

    // region
    num_of_selected_regions: usize,
    num_of_mapped_reads_in_regions: usize,

    // insert size
    mean_insert_size: f64,
    std_insert_size: f64,
    median_insert_size: usize,

    // mapping quality
    mean_mapping_quality: f64,

    // ACTG content
    num_of_as: usize,
    percent_of_as: f64,
    num_of_cs: usize,
    percent_of_cs: f64,
    num_of_ts: usize,
    percent_of_ts: f64,
    num_of_gs: usize,
    percent_of_gs: f64,
    num_of_ns: usize,
    percent_of_ns: f64,
    gc_content: f64,

    // mismatches and indels
    general_error_rate: f64,
    num_of_indels: usize,
    num_of_mismatches: usize,
    num_of_insertions: usize,
    num_of_deletions: usize,
    mapped_reads_with_insertions_fraction: f64,
    mapped_reads_with_deletion_fraction: f64,
    homopolymer_indels_fraction: f64,

    // Coverage
    mean_coverage_data: f64,
    std_coverage_data: f64,
    num_overlap_reads_pairs: usize,
    paired_end_adapted_mean_coverage: f64,
    coverage_quotes: XYVector,
    // Coverage per contig
    chromosome_stats: Vec<ChromosomeInfo>,
}

impl DetailedStats {
    pub fn new() -> Self {
        Self {
            // input
            bam_file: "".to_string(),

            // reference
            num_of_bases: 0,
            num_of_contigs: 0,

            // globals
            num_of_windows: 0,
            num_of_reads: 0,
            num_of_mapped_reads: 0,
            percent_of_mapped_reads: 0.0,
            num_of_supplementary_alignments: 0,
            percent_of_supplementary_alignments: 0.0,
            num_of_second_alignments: 0,
            num_of_paired_reads: 0,
            num_of_mapped_reads_first_in_pair: 0,
            num_of_mapped_reads_second_in_pair: 0,
            num_of_mapped_reads_both_in_pair: 0,
            num_of_mapped_reads_singletons: 0,
            num_of_overlapping_reads: 0,
            num_of_mapped_bases: 0,
            num_of_sequenced_bases: 0,
            num_of_aligned_bases: 0,

            num_of_flagged_duplicate_reads: 0,
            num_of_estimated_duplicate_reads: 0,
            duplication_rate: 0.0,

            // region
            num_of_selected_regions: 0,
            num_of_mapped_reads_in_regions: 0,

            // insert size
            mean_insert_size: 0.0,
            std_insert_size: 0.0,
            median_insert_size: 0,

            // mapping quality
            mean_mapping_quality: 0.0,

            // ACTG content
            num_of_as: 0,
            percent_of_as: 0.0,
            num_of_cs: 0,
            percent_of_cs: 0.0,
            num_of_ts: 0,
            percent_of_ts: 0.0,
            num_of_gs: 0,
            percent_of_gs: 0.0,
            num_of_ns: 0,
            percent_of_ns: 0.0,
            gc_content: 0.0,

            // mismatches and indels
            general_error_rate: 0.0,
            num_of_indels: 0,
            num_of_mismatches: 0,
            num_of_insertions: 0,
            num_of_deletions: 0,
            mapped_reads_with_insertions_fraction: 0.0,
            mapped_reads_with_deletion_fraction: 0.0,
            homopolymer_indels_fraction: 0.0,

            // Coverage
            mean_coverage_data: 0.0,
            std_coverage_data: 0.0,
            num_overlap_reads_pairs: 0,
            paired_end_adapted_mean_coverage: 0.0,
            coverage_quotes: XYVector::new(),
            // Coverage per contig
            chromosome_stats: vec![],
        }
    }

    fn extract(&mut self, bam_stats: &BamStats) {
        unsafe {
            // input
            self.bam_file = bam_stats.source_file.clone();

            // reference
            self.num_of_bases = bam_stats.reference_size as usize;
            self.num_of_contigs = bam_stats.number_of_reference_contigs as usize;

            // globals
            self.num_of_windows = NUMBER_OF_TOTAL_WINDOWS as usize;
            self.num_of_reads = bam_stats.number_of_reads as usize;
            self.num_of_mapped_reads = bam_stats.number_of_mapped_reads as usize;
            self.percent_of_mapped_reads = bam_stats.get_percentage_of_mapped_reads() as f64;
            self.num_of_supplementary_alignments = bam_stats.number_of_supp_alignments as usize;
            self.percent_of_supplementary_alignments =
                bam_stats.get_percentage_of_supp_alignments() as f64;
            self.num_of_second_alignments = bam_stats.number_of_secondary_alignments as usize;
            self.num_of_paired_reads = bam_stats.number_of_paired_reads as usize;
            self.num_of_mapped_reads_first_in_pair =
                bam_stats.number_of_mapped_first_of_pair as usize;
            self.num_of_mapped_reads_second_in_pair =
                bam_stats.number_of_mapped_second_of_pair as usize;

            self.num_of_mapped_reads_singletons = bam_stats.number_of_singletons as usize;
            self.num_of_mapped_reads_both_in_pair =
                self.num_of_paired_reads - self.num_of_mapped_reads_singletons;
            self.num_of_overlapping_reads = bam_stats.num_overlapping_read_pairs as usize;
            self.num_of_mapped_bases = bam_stats.number_of_mapped_bases as usize;
            self.num_of_sequenced_bases = bam_stats.number_of_sequenced_bases as usize;
            self.num_of_aligned_bases = bam_stats.number_of_aligned_bases as usize;

            self.num_of_flagged_duplicate_reads = bam_stats.num_detected_duplicate_reads as usize;
            self.num_of_estimated_duplicate_reads =
                bam_stats.num_estimated_duplicate_reads as usize;
            self.duplication_rate = bam_stats.duplication_rate;

            // region
            self.num_of_selected_regions = bam_stats.num_selected_regions as usize;
            self.num_of_mapped_reads_in_regions =
                bam_stats.number_of_mapped_reads_in_regions as usize;

            // insert size
            self.mean_insert_size = bam_stats.mean_insert_size;
            self.std_insert_size = bam_stats.std_insert_size;
            self.median_insert_size = bam_stats.median_insert_size as usize;

            // mapping quality
            self.mean_mapping_quality = bam_stats.mean_mapping_quality_per_window;

            // ACTG content
            self.num_of_as = bam_stats.number_of_as as usize;
            self.num_of_cs = bam_stats.number_of_cs as usize;
            self.num_of_ts = bam_stats.number_of_ts as usize;
            self.num_of_gs = bam_stats.number_of_gs as usize;
            self.num_of_ns = bam_stats.number_of_ns as usize;
            self.percent_of_as = bam_stats.mean_a_relative_content;
            self.percent_of_ts = bam_stats.mean_t_relative_content;
            self.percent_of_cs = bam_stats.mean_c_relative_content;
            self.percent_of_gs = bam_stats.mean_g_relative_content;
            self.percent_of_ns = bam_stats.mean_n_relative_content;
            self.gc_content = bam_stats.mean_gc_relative_content;

            // mismatches and indels
            self.general_error_rate = bam_stats.get_error_rate();
            self.num_of_indels = bam_stats.get_num_indels();
            self.num_of_mismatches = bam_stats.num_mismatches as usize;
            self.num_of_insertions = bam_stats.get_num_insertions();
            self.num_of_deletions = bam_stats.get_num_deletions();
            self.mapped_reads_with_insertions_fraction =
                bam_stats.get_reads_with_insertion_percentage();
            self.mapped_reads_with_deletion_fraction =
                bam_stats.get_reads_with_deletions_percentage();
            self.homopolymer_indels_fraction = bam_stats.get_homopolymer_indels_fraction() * 100.0;

            // Coverage
            self.mean_coverage_data = bam_stats.mean_coverage;
            self.std_coverage_data = bam_stats.std_coverage;
            self.num_overlap_reads_pairs = bam_stats.num_overlapping_read_pairs as usize;
            self.paired_end_adapted_mean_coverage = bam_stats.adapted_mean_coverage;
            self.coverage_quotes = bam_stats.coverage_quotes.clone();
            // Coverage per contig
            self.chromosome_stats = bam_stats.chromosome_stats.clone();
        }
    }

    pub fn export_to_txt(&self, output_dir: &str) -> Result<(), Box<dyn Error>> {
        let filepath = Path::new(output_dir).join(format!("{}.txt", "genome_results"));
        let mut report = File::create(filepath)?;
        writeln!(report, "BamQC report")?;
        writeln!(report, "-----------------------------------")?;
        writeln!(report)?;

        writeln!(report, ">>>>>>> Input")?;
        writeln!(report)?;
        writeln!(report, "     bam file = {}", self.bam_file.clone())?;
        writeln!(
            report,
            "     outfile = {}",
            format!("{}/{}.txt", output_dir, "genome_results")
        )?;
        writeln!(report)?;
        writeln!(report)?;

        writeln!(report, ">>>>>>> Reference")?;
        writeln!(report)?;
        writeln!(report, "     number of bases = {}", self.num_of_bases)?;
        writeln!(report, "     number of contigs = {}", self.num_of_contigs)?;
        writeln!(report)?;
        writeln!(report)?;

        writeln!(report, ">>>>>>> Globals")?;
        writeln!(report)?;
        writeln!(report, "     number of windows = {}", self.num_of_windows)?;
        writeln!(report)?;
        writeln!(report, "     number of reads = {}", self.num_of_reads)?;
        writeln!(
            report,
            "     number of mapped reads = {} ({})",
            self.num_of_mapped_reads,
            format!("{:.2}%", self.percent_of_mapped_reads)
        )?;

        if self.num_of_supplementary_alignments > 0 {
            writeln!(
                report,
                "     number of supplementary alignments = {} ({})",
                self.num_of_supplementary_alignments,
                format!("{:.2}%", self.percent_of_supplementary_alignments)
            )?;
        }

        writeln!(
            report,
            "     number of secondary alignments = {}",
            self.num_of_second_alignments
        )?;

        writeln!(report)?;

        if self.num_of_paired_reads > 0 {
            writeln!(
                report,
                "     number of mapped paired reads (first in pair) = {}",
                self.num_of_mapped_reads_first_in_pair
            )?;
            writeln!(
                report,
                "     number of mapped paired reads (second in pair) = {}",
                self.num_of_mapped_reads_second_in_pair
            )?;
            writeln!(
                report,
                "     number of mapped paired reads (both in pair) = {}",
                self.num_of_paired_reads - self.num_of_mapped_reads_singletons
            )?;
            writeln!(
                report,
                "     number of mapped paired reads (singletons) = {}",
                self.num_of_mapped_reads_singletons
            )?;

            if self.num_overlap_reads_pairs > 0 {
                writeln!(
                    report,
                    "     number of overlapping read pairs = {}",
                    self.num_overlap_reads_pairs
                )?;
            }

            writeln!(report)?;
        }

        writeln!(
            report,
            "     number of mapped bases = {} bp",
            self.num_of_mapped_bases
        )?;
        writeln!(
            report,
            "     number of sequenced bases = {} bp",
            self.num_of_sequenced_bases
        )?;
        writeln!(
            report,
            "     number of aligned bases = {} bp",
            self.num_of_aligned_bases
        )?;

        if self.num_of_flagged_duplicate_reads > 0 {
            writeln!(
                report,
                "     number of duplicated reads (flagged) = {}",
                self.num_of_flagged_duplicate_reads
            )?;
        } else {
            writeln!(
                report,
                "     number of duplicated reads (estimated) = {}",
                self.num_of_estimated_duplicate_reads
            )?;
            writeln!(
                report,
                "     duplication rate = {}",
                format!("{:.2}%", self.duplication_rate)
            )?;
        }

        if self.num_of_mapped_reads == 0
            || (self.num_of_selected_regions > 0 && self.num_of_mapped_reads_in_regions == 0)
        {
            return Ok(());
        }

        writeln!(report)?;
        writeln!(report)?;

        // todo!() region analyze

        // insert size
        writeln!(report, ">>>>>>> Insert size")?;
        writeln!(report)?;
        writeln!(
            report,
            "     mean insert size = {}",
            format!("{:.4}", self.mean_insert_size)
        )?;
        writeln!(
            report,
            "     std insert size = {}",
            format!("{:.4}", self.std_insert_size)
        )?;
        writeln!(
            report,
            "     median insert size = {}",
            self.median_insert_size
        )?;
        writeln!(report)?;
        writeln!(report)?;

        // mapping quality
        writeln!(report, ">>>>>>> Mapping quality")?;
        writeln!(report)?;
        writeln!(
            report,
            "     mean mapping quality = {}",
            format!("{:.4}", self.mean_mapping_quality)
        )?;
        writeln!(report)?;
        writeln!(report)?;

        // actg content
        writeln!(report, ">>>>>>> ACTG content")?;
        writeln!(report)?;
        writeln!(
            report,
            "     number of A's = {} bp ({})",
            self.num_of_as,
            format!("{:.2}%", self.percent_of_as)
        )?;
        writeln!(
            report,
            "     number of C's = {} bp ({})",
            self.num_of_cs,
            format!("{:.2}%", self.percent_of_cs)
        )?;
        writeln!(
            report,
            "     number of T's = {} bp ({})",
            self.num_of_ts,
            format!("{:.2}%", self.percent_of_ts)
        )?;
        writeln!(
            report,
            "     number of G's = {} bp ({})",
            self.num_of_gs,
            format!("{:.2}%", self.percent_of_gs)
        )?;
        writeln!(
            report,
            "     number of N's = {} bp ({})",
            self.num_of_ns,
            format!("{:.2}%", self.percent_of_ns)
        )?;
        writeln!(report)?;
        writeln!(
            report,
            "     GC percentage = {}",
            format!("{:.2}%", self.gc_content)
        )?;

        writeln!(report)?;
        writeln!(report)?;

        // Mismatches and indels
        writeln!(report, ">>>>>>> Mismatches and indels")?;
        writeln!(report)?;
        writeln!(
            report,
            "    general error rate = {}",
            format!("{:.4}", self.general_error_rate)
        )?;
        writeln!(
            report,
            "    number of mismatches = {}",
            self.num_of_mismatches
        )?;

        if self.num_of_indels > 0 {
            writeln!(
                report,
                "    number of insertions = {}",
                self.num_of_insertions
            )?;
            writeln!(
                report,
                "    mapped reads with insertion percentage = {}",
                format!("{:.2}%", self.mapped_reads_with_insertions_fraction)
            )?;
            writeln!(
                report,
                "    number of deletions = {}",
                self.num_of_deletions
            )?;
            writeln!(
                report,
                "    mapped reads with deletion percentage = {}",
                format!("{:.2}%", self.mapped_reads_with_deletion_fraction)
            )?;
            writeln!(
                report,
                "    homopolymer indels = {}",
                format!("{:.2}%", self.homopolymer_indels_fraction)
            )?;
        }

        writeln!(report)?;
        writeln!(report)?;

        // coverageData
        writeln!(report, ">>>>>>> Coverage")?;
        writeln!(report)?;
        writeln!(
            report,
            "     mean coverageData = {:.4}X",
            self.mean_coverage_data
        )?;
        writeln!(
            report,
            "     std coverageData = {:.4}X",
            self.std_coverage_data
        )?;

        if self.num_overlap_reads_pairs > 0 {
            writeln!(
                report,
                "     paired-end adapted mean coverageData = {}X",
                format!("{:.2}%", self.paired_end_adapted_mean_coverage)
            )?;
        }

        writeln!(report)?;
        for i in 0..self.coverage_quotes.get_size() {
            writeln!(
                report,
                "     There is a {:.2}% of reference with a coverage >= {}X",
                self.coverage_quotes.get(i).unwrap().get_y(),
                self.coverage_quotes.get(i).unwrap().get_x()
            )?;
        }
        writeln!(report)?;

        // Coverage per contig
        writeln!(report, ">>>>>>> Coverage per contig")?;
        writeln!(report)?;

        for chromosome_info in &self.chromosome_stats {
            writeln!(
                report,
                "\t{}\t{}\t{}\t{}\t{}",
                chromosome_info.get_name(),
                chromosome_info.get_length(),
                chromosome_info.get_num_bases(),
                chromosome_info.get_cov_mean(),
                chromosome_info.get_cov_std()
            )?;
        }
        writeln!(report)?;

        // report.close()?;
        Ok(())
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
    // all coverage histogram
    all_coverage_histogram_x_label: String,
    all_coverage_histogram_x_data: Vec<f64>,
    all_coverage_histogram_y_label: String,
    all_coverage_histogram_y_data: Vec<f64>,
    // balanced coverage histogram
    balanced_coverage_bar_names_x_label: String,
    balanced_coverage_bar_names_x_data: Vec<String>,
    balance_coverage_histogram_y_label: String,
    balance_coverage_histogram_y_data: Vec<f64>,
    // Coverage Histogram 0-50x
    ranged_coverage_histogram_x_label: String,
    ranged_coverage_histogram_x_data: Vec<f64>,
    ranged_coverage_histogram_y_label: String,
    ranged_coverage_histogram_y_data: Vec<f64>,
    // Genome Fraction Coverage
    genome_fraction_coverage_x_label: String,
    genome_fraction_coverage_x_data: Vec<f64>,
    genome_fraction_coverage_y_label: String,
    genome_fraction_coverage_y_data: Vec<f64>,
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
            // all coverage histogram
            all_coverage_histogram_x_label: "".to_string(),
            all_coverage_histogram_x_data: vec![],
            all_coverage_histogram_y_label: "".to_string(),
            all_coverage_histogram_y_data: vec![],
            // balanced coverage histogram
            balanced_coverage_bar_names_x_label: "".to_string(),
            balanced_coverage_bar_names_x_data: vec![],
            balance_coverage_histogram_y_label: "".to_string(),
            balance_coverage_histogram_y_data: vec![],
            // Coverage Histogram 0-50x
            ranged_coverage_histogram_x_label: "".to_string(),
            ranged_coverage_histogram_x_data: vec![],
            ranged_coverage_histogram_y_label: "".to_string(),
            ranged_coverage_histogram_y_data: vec![],
            // Genome Fraction Coverage
            genome_fraction_coverage_x_label: "".to_string(),
            genome_fraction_coverage_x_data: vec![],
            genome_fraction_coverage_y_label: "".to_string(),
            genome_fraction_coverage_y_data: vec![],
        }
    }

    pub fn export_to_txt(&self, output_dir: &str) -> Result<(), Box<dyn Error>> {
        self.export_coverage_across_reference(output_dir);
        self.export_coverage_histogram(output_dir);
        Ok(())
    }

    fn export_coverage_histogram(&self, output_dir: &str) -> Result<(), Box<dyn Error>> {
        let filepath = Path::new(output_dir).join(format!("{}.txt", "coverage_histogram"));
        let mut f = File::create(filepath).unwrap();
        let mut writer = WriterBuilder::new().delimiter(b'\t').from_writer(f);
        writer.write_record(&["#Coverage", "Number of genomic locations"])?;

        for (x, y) in self
            .all_coverage_histogram_x_data
            .iter()
            .zip(self.all_coverage_histogram_y_data.iter())
        {
            writer.write_record(&[x.to_string(), y.to_string()])?
        }
        Ok(())
    }

    fn export_coverage_across_reference(&self, output_dir: &str) -> Result<(), Box<dyn Error>> {
        let filepath = Path::new(output_dir).join(format!("{}.txt", "coverage_across_reference"));
        let mut f = File::create(filepath).unwrap();
        let mut writer = WriterBuilder::new().delimiter(b'\t').from_writer(f);
        writer.write_record(&["#Position (bp)", "Coverage", "Std"])?;

        for (x, (y, std)) in self.coverage_across_reference_x_data.iter().zip(
            self.coverage_across_reference_y_data
                .iter()
                .zip(self.coverage_across_reference_deviation_std_y_data.iter()),
        ) {
            writer.write_record(&[x.to_string(), y.to_string(), std.to_string()])?
        }
        Ok(())
    }

    fn extract(&mut self, bam_stats: &BamStats) {
        unsafe {
            // coverage across reference
            self.coverage_across_reference_x_label = "Position (bp)".to_string();
            self.coverage_across_reference_x_data = vec![0.0; NUMBER_OF_TOTAL_WINDOWS as usize];
            for i in 0..NUMBER_OF_TOTAL_WINDOWS as usize {
                unsafe {
                    self.coverage_across_reference_x_data[i] =
                        (WINDOW_STARTS[i] + WINDOW_ENDS[i]) as f64 / 2.0;
                }
            }

            self.coverage_across_reference_y_label = "Coverage (X)".to_string();
            self.coverage_across_reference_y_data = bam_stats.coverage_across_reference.clone();
            self.coverage_across_reference_deviation_std_y_data =
                bam_stats.std_coverage_across_reference.clone();

            // all coverage histogram
            self.all_coverage_histogram_y_label = "Number of genomic locations".to_string();
            self.all_coverage_histogram_x_label = "Coverage (X)".to_string();
            self.all_coverage_histogram_x_data = bam_stats.coverage_histogram.get_x_vector();
            self.all_coverage_histogram_y_data = bam_stats.coverage_histogram.get_y_vector();

            // balanced coverage histogram
            self.balance_coverage_histogram_y_label = "Number of genomic locations".to_string();
            self.balance_coverage_histogram_y_data =
                bam_stats.balanced_coverage_histogram.get_y_vector();

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
                        .push(item.xy_item.get_y());
                    self.ranged_coverage_histogram_x_data
                        .push(item.xy_item.get_x());
                }
            }

            // Genome Fraction Coverage
            self.genome_fraction_coverage_y_data = bam_stats.coverage_quotes.get_y_vector();
            self.genome_fraction_coverage_x_data = bam_stats.coverage_quotes.get_x_vector();
            self.genome_fraction_coverage_x_label = "Coverage (X)".to_string();
            self.genome_fraction_coverage_y_label = "Fraction of reference (%)".to_string();
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct DuplicationStats {
    duplication_rate_histogram_x_label: String,
    duplication_rate_histogram_y_label: String,
    duplication_rate_histogram_x_data: Vec<f64>,
    duplication_rate_histogram_y_data: Vec<f64>,
}

impl DuplicationStats {
    pub fn new() -> Self {
        Self {
            duplication_rate_histogram_x_label: "".to_string(),
            duplication_rate_histogram_y_label: "".to_string(),
            duplication_rate_histogram_x_data: vec![],
            duplication_rate_histogram_y_data: vec![],
        }
    }

    fn extract(&mut self, bam_stats: &BamStats) {
        self.duplication_rate_histogram_x_label = "Duplication rate".to_string();
        self.duplication_rate_histogram_y_label = "Number of loci".to_string();
        self.duplication_rate_histogram_x_data =
            bam_stats.unique_read_starts_histogram.get_x_vector();
        self.duplication_rate_histogram_y_data =
            bam_stats.unique_read_starts_histogram.get_y_vector();
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct MappedReadsNucleotideContent {
    reads_nuleotide_a_content_x_data: Vec<f64>,
    reads_nuleotide_a_content_y_data: Vec<f64>,
    reads_nuleotide_t_content_x_data: Vec<f64>,
    reads_nuleotide_t_content_y_data: Vec<f64>,
    reads_nuleotide_c_content_x_data: Vec<f64>,
    reads_nuleotide_c_content_y_data: Vec<f64>,
    reads_nuleotide_g_content_x_data: Vec<f64>,
    reads_nuleotide_g_content_y_data: Vec<f64>,
    reads_nuleotide_n_content_x_data: Vec<f64>,
    reads_nuleotide_n_content_y_data: Vec<f64>,
}

impl MappedReadsNucleotideContent {
    pub fn new() -> Self {
        Self {
            reads_nuleotide_a_content_x_data: vec![],
            reads_nuleotide_a_content_y_data: vec![],
            reads_nuleotide_t_content_x_data: vec![],
            reads_nuleotide_t_content_y_data: vec![],
            reads_nuleotide_c_content_x_data: vec![],
            reads_nuleotide_c_content_y_data: vec![],
            reads_nuleotide_g_content_x_data: vec![],
            reads_nuleotide_g_content_y_data: vec![],
            reads_nuleotide_n_content_x_data: vec![],
            reads_nuleotide_n_content_y_data: vec![],
        }
    }

    fn extract(&mut self, bam_stats: &BamStats) {
        self.reads_nuleotide_a_content_x_data = bam_stats.reads_as_histogram.get_x_vector();
        self.reads_nuleotide_a_content_y_data = bam_stats.reads_as_histogram.get_y_vector();
        self.reads_nuleotide_t_content_x_data = bam_stats.reads_ts_histogram.get_x_vector();
        self.reads_nuleotide_t_content_y_data = bam_stats.reads_ts_histogram.get_y_vector();
        self.reads_nuleotide_c_content_x_data = bam_stats.reads_cs_histogram.get_x_vector();
        self.reads_nuleotide_c_content_y_data = bam_stats.reads_cs_histogram.get_y_vector();
        self.reads_nuleotide_g_content_x_data = bam_stats.reads_gs_histogram.get_x_vector();
        self.reads_nuleotide_g_content_y_data = bam_stats.reads_gs_histogram.get_y_vector();
        self.reads_nuleotide_n_content_x_data = bam_stats.reads_ns_histogram.get_x_vector();
        self.reads_nuleotide_n_content_y_data = bam_stats.reads_ns_histogram.get_y_vector();
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct MappedReadsGCContent {
    reads_gc_content_histogram_x_label: String,
    reads_gc_content_histogram_y_label: String,
    reads_gc_content_histogram_x_data: Vec<f64>,
    reads_gc_content_histogram_y_data: Vec<f64>,
}

impl MappedReadsGCContent {
    pub fn new() -> Self {
        Self {
            reads_gc_content_histogram_x_label: "".to_string(),
            reads_gc_content_histogram_y_label: "".to_string(),
            reads_gc_content_histogram_x_data: vec![],
            reads_gc_content_histogram_y_data: vec![],
        }
    }

    pub fn export_to_txt(&self, output_dir: &str) -> Result<(), Box<dyn Error>> {
        self.export_gc_content_histogram(output_dir);
        Ok(())
    }

    fn export_gc_content_histogram(&self, output_dir: &str) -> Result<(), Box<dyn Error>> {
        let filepath =
            Path::new(output_dir).join(format!("{}.txt", "mapped_reads_gc-content_distribution"));
        let mut f = File::create(filepath).unwrap();
        let mut writer = WriterBuilder::new().delimiter(b'\t').from_writer(f);
        writer.write_record(&["#GC Content (%)", "Sample"])?;

        for (x, y) in self
            .reads_gc_content_histogram_x_data
            .iter()
            .zip(self.reads_gc_content_histogram_y_data.iter())
        {
            writer.write_record(&[x.to_string(), y.to_string()])?
        }
        Ok(())
    }

    fn extract(&mut self, bam_stats: &BamStats) {
        // todo!() comparison of reference genome gc content
        self.reads_gc_content_histogram_x_label = "GC Content (%)".to_string();
        self.reads_gc_content_histogram_y_label = "Fraction of reads".to_string();

        self.reads_gc_content_histogram_x_data =
            bam_stats.get_gc_content_histogram().get_x_vector();
        self.reads_gc_content_histogram_y_data =
            bam_stats.get_gc_content_histogram().get_y_vector();
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct ClippingProfileStats {
    clipping_profile_x_label: String,
    clipping_profile_y_label: String,
    clipping_profile_x_data: Vec<f64>,
    clipping_profile_y_data: Vec<f64>,
}
impl ClippingProfileStats {
    pub fn new() -> Self {
        Self {
            clipping_profile_x_label: "".to_string(),
            clipping_profile_y_label: "".to_string(),
            clipping_profile_x_data: vec![],
            clipping_profile_y_data: vec![],
        }
    }

    fn extract(&mut self, bam_stats: &BamStats) {
        self.clipping_profile_x_label = "Read position (bp)".to_string();
        self.clipping_profile_y_label = " Clipped bases (%)".to_string();
        self.clipping_profile_x_data = bam_stats.reads_clipping_profile_histogram.get_x_vector();
        self.clipping_profile_y_data = bam_stats.reads_clipping_profile_histogram.get_y_vector();
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
    mapping_quality_histogram_x_data: Vec<f64>,
    mapping_quality_histogram_y_data: Vec<f64>,
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
            mapping_quality_histogram_x_data: vec![],
            mapping_quality_histogram_y_data: vec![],
        }
    }

    fn extract(&mut self, bam_stats: &BamStats) {
        unsafe {
            // mapping quality across reference
            self.mapping_quality_across_reference_x_label = "Position (bp)".to_string();
            self.mapping_quality_across_reference_y_label = "Mapping quality".to_string();
            self.mapping_quality_across_reference_y_data =
                bam_stats.mapping_quality_across_reference.clone();
            self.mapping_quality_across_reference_x_data =
                vec![0.0; NUMBER_OF_TOTAL_WINDOWS as usize];
            for i in 0..NUMBER_OF_TOTAL_WINDOWS as usize {
                self.mapping_quality_across_reference_x_data[i] =
                    (WINDOW_STARTS[i] + WINDOW_ENDS[i]) as f64 / 2.0;
            }

            // mapping quality histogram
            self.mapping_quality_histogram_x_label = "Mapping quality".to_string();
            self.mapping_quality_histogram_y_label = "Number of genomic locations".to_string();
            self.mapping_quality_histogram_x_data =
                bam_stats.mapping_quality_histogram.get_x_vector();
            self.mapping_quality_histogram_y_data =
                bam_stats.mapping_quality_histogram.get_y_vector();
        }
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
    insert_size_histogram_x_data: Vec<f64>,
    insert_size_histogram_y_data: Vec<f64>,
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
            insert_size_histogram_x_data: vec![],
            insert_size_histogram_y_data: vec![],
        }
    }

    pub fn export_to_txt(&self, output_dir: &str) -> Result<(), Box<dyn Error>> {
        self.export_insert_size_across_reference(output_dir);
        self.export_insert_size_histogram(output_dir);
        Ok(())
    }

    fn export_insert_size_histogram(&self, output_dir: &str) -> Result<(), Box<dyn Error>> {
        let filepath = Path::new(output_dir).join(format!("{}.txt", "insert_size_histogram"));
        let mut f = File::create(filepath).unwrap();
        let mut writer = WriterBuilder::new().delimiter(b'\t').from_writer(f);
        writer.write_record(&["#Insert size (bp)", "insert size"])?;

        for (x, y) in self
            .insert_size_histogram_x_data
            .iter()
            .zip(self.insert_size_histogram_y_data.iter())
        {
            writer.write_record(&[x.to_string(), y.to_string()])?
        }
        Ok(())
    }

    fn export_insert_size_across_reference(&self, output_dir: &str) -> Result<(), Box<dyn Error>> {
        let filepath =
            Path::new(output_dir).join(format!("{}.txt", "insert_size_across_reference"));
        let mut f = File::create(filepath).unwrap();

        let mut writer = WriterBuilder::new().delimiter(b'\t').from_writer(f);
        writer.write_record(&["#Position (bp)", "insert size"])?;

        for (x, y) in self
            .insert_size_across_reference_x_data
            .iter()
            .zip(self.insert_size_across_reference_y_data.iter())
        {
            writer.write_record(&[x.to_string(), y.to_string()])?
        }
        Ok(())
    }

    fn extract(&mut self, bam_stats: &BamStats) {
        unsafe {
            // insert size across reference
            self.insert_size_across_reference_x_label = "Position (bp)".to_string();
            self.insert_size_across_reference_y_label = "Insert size (bp)".to_string();
            self.insert_size_across_reference_y_data =
                bam_stats.insert_size_across_reference.clone();
            self.insert_size_across_reference_x_data = vec![0.0; NUMBER_OF_TOTAL_WINDOWS as usize];
            for i in 0..NUMBER_OF_TOTAL_WINDOWS as usize {
                unsafe {
                    self.insert_size_across_reference_x_data[i] =
                        (WINDOW_STARTS[i] + WINDOW_ENDS[i]) as f64 / 2.0;
                }
            }

            // insert size histogram
            self.insert_size_histogram_x_label = "Insert size (bp)".to_string();
            self.insert_size_histogram_y_label = "Number of reads".to_string();
            self.insert_size_histogram_x_data = bam_stats.insert_size_histogram.get_x_vector();
            self.insert_size_histogram_y_data = bam_stats.insert_size_histogram.get_y_vector();
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Qualimap {
    pub input_info: InputInfo,
    pub basic_stats: BasicStats,
    #[serde(skip_serializing)]
    pub detailed_stats: DetailedStats,
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
            detailed_stats: DetailedStats::new(),
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

    pub fn run(&mut self, bam_file: &str, qualimap_config: &QualimapConfig) {
        let mut bam_stats_analysis = BamStatsAnalysis::new(bam_file.to_string(), qualimap_config);

        println!("Running BAM file analysis...");
        bam_stats_analysis.run();
        println!("End of bam qc");
        println!("Computing report...");

        let bam_stats = bam_stats_analysis.get_bam_stats().unwrap();
        self.input_info.extract(&bam_stats_analysis);
        self.basic_stats.extract(bam_stats);
        self.detailed_stats.extract(bam_stats);
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

    pub fn export(&self,output_dir:&str) -> Result<(), Box<dyn Error>>{
        let _ =self.detailed_stats.export_to_txt(output_dir);
        let _ =self.coverage_stats.export_to_txt(output_dir);
        let _ =self.insert_size_stats.export_to_txt(output_dir);
        let _ =self.mapping_read_gc_content.export_to_txt(output_dir);

        Ok(())
    }
}
