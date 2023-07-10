use crate::qc::config::bam_config::{QualimapConfig, QualimapConstants};
use linear_map::LinearMap;
use math::stats;
use regex::Regex;
use rust_htslib::bam::record::{Aux};
use rust_htslib::bam::{Header, Read, Reader, Record};
use rust_htslib::tpool;
use serde::{Deserialize, Serialize};
use std::sync::Arc;
use std::thread;
use std::{
    collections::HashMap,
    f64::consts::{E, PI},
    vec,
};
use time::{Date, Duration, Instant};

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
struct BamGenomeWindow {
    start: usize,
    end: usize,
    number_of_as: usize,
    number_of_ts: usize,
    number_of_cs: usize,
    number_of_gs: usize,
    number_of_ns: usize,
    number_of_mapped_bases: usize,
    correct_insert_sizes: usize,
    acum_insert_sizes: f64,
    acum_mapping_quality: f64,
}
impl BamGenomeWindow {
    pub fn new(_start: usize, _end: usize) -> Self {
        Self {
            start: _start,
            end: _end,
            number_of_as: 0,
            number_of_ts: 0,
            number_of_cs: 0,
            number_of_gs: 0,
            number_of_ns: 0,
            number_of_mapped_bases: 0,
            correct_insert_sizes: 0,
            acum_insert_sizes: 0.0,
            acum_mapping_quality: 0.0,
        }
    }

    fn merge(&mut self, other: &Self) {
        self.number_of_as += other.number_of_as;
        self.number_of_ts += other.number_of_ts;
        self.number_of_cs += other.number_of_cs;
        self.number_of_gs += other.number_of_gs;
        self.number_of_ns += other.number_of_ns;
        self.number_of_mapped_bases += other.number_of_mapped_bases;
        self.correct_insert_sizes += other.correct_insert_sizes;
        self.acum_insert_sizes += other.acum_insert_sizes;
        self.acum_mapping_quality += other.acum_mapping_quality;
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
struct FastBamAnalysis {
    bam_file: String,
    thread_num: usize,
    number_of_windows: usize,
    total_num_of_records: usize,

    window_size: usize,
    effective_number_of_window: usize,
    window_starts: Vec<usize>,
    window_ends: Vec<usize>,
    windows: Vec<BamGenomeWindow>,

    number_of_problematic_reads: usize,
    number_of_valid_reads: usize,
    number_of_correct_strand_reads: usize,

    sum_coverage_squared: usize,
    sum_coverage: usize,

    // reference
    reference_file: String,
    reference_available_flag: bool,
    reference_size: usize,
    number_of_reference_config: usize,
    acum_read_size: usize,
    max_read_size: usize,
    min_read_size: usize,
    number_of_secondary_alignments: usize,
    number_of_reads: usize,
    selected_regions_available_flag: bool,
    num_supplementary_alignments: usize,
    num_mapped_reads: usize,
    num_paired_reads: usize,
    num_mapped_first_in_pair: usize,
    num_mapped_second_in_pair: usize,
    num_singletons: usize,
    num_marked_duplicates: usize,

    skip_marked_duplicates_flag: bool,
    num_of_duplicates_skip: usize,
    skip_detected_duplicates_flag: bool,
    num_estimated_duplicate_reads: usize,
    collect_intersecting_paired_end_reads_flag: bool,

    // read starts (duplicate)
    read_starts_histogram: ReadStartsHistogram,
    unique_read_starts_histogram: XYVector,
    duplication_rate: f64,

    // coordinates transformer
    locator: GenomeLocator,
    num_insertions: usize,
    num_read_with_insertion: usize,
    num_deletions: usize,
    num_read_with_deletion: usize,

    edit_distance: usize,
    num_mismatches: usize,
    num_clipped_reads: usize,

    // read content
    reads_a_content: Vec<usize>,
    reads_c_content: Vec<usize>,
    reads_g_content: Vec<usize>,
    reads_t_content: Vec<usize>,
    reads_n_content: Vec<usize>,
    reads_as_histogram: XYVector,
    reads_cs_histogram: XYVector,
    reads_gs_histogram: XYVector,
    reads_ts_histogram: XYVector,
    reads_ns_histogram: XYVector,
    a_content_across_reference: Vec<f64>,
    t_content_across_reference: Vec<f64>,
    c_content_across_reference: Vec<f64>,
    g_content_across_reference: Vec<f64>,
    n_content_across_reference: Vec<f64>,
    gc_content_across_reference: Vec<f64>,
    a_relative_content_across_reference: Vec<f64>,
    c_relative_content_across_reference: Vec<f64>,
    t_relative_content_across_reference: Vec<f64>,
    g_relative_content_across_reference: Vec<f64>,
    n_relative_content_across_reference: Vec<f64>,
    gc_relative_content_across_reference: Vec<f64>,

    // clipping content
    reads_clipping_content: Vec<usize>,
    reads_clipping_profile_histogram: XYVector,

    reads_gc_content: Vec<f32>,
    sample_count: usize,
    homopolymer_indels: Vec<usize>,

    num_bases: usize,
    num_gc: usize,
    min_homopolymer_size: usize,
    prev_base_inside_indel_region_flag: bool,
    homopolymer_starts_inside_indel_region_flag: bool,
    prev_base: u8,
    homopolymer_size: usize,

    num_read_with_correct_insert_sizes: usize,
    acum_insert_size: f64,
    number_of_reads_with_start_greater_than_end: usize,

    number_of_sequenced_bases: usize,
    number_of_mapped_bases: usize,
    number_of_as: usize,
    number_of_ts: usize,
    number_of_cs: usize,
    number_of_gs: usize,
    number_of_ns: usize,

    pg_program: String,
    pg_command_string: String,

    // gc content histogram
    gc_content_histogram: Vec<f64>,

    coverage_per_position: Vec<Vec<i64>>,
    coverage_across_reference: Vec<f64>,
    std_coverage_across_reference: Vec<f64>,
    coverage_per_window: Vec<usize>,
    coverage_squared_per_window: Vec<usize>,
    coverage_histogram_map: HashMap<i64, i64>,
    coverage_histogram_cache: Vec<i64>,
    coverage_histogram: XYVector,
    acum_coverage_histogram: XYVector,
    balanced_coverage_histogram: XYVector,
    balanced_coverage_bar_names: HashMap<i64, String>,
    max_coverage_quota: i32,
    coverage_quotes: XYVector,

    mapping_quality_per_position: Vec<Vec<i64>>,
    mapping_quality_across_reference: Vec<f64>,
    mapping_quality_histogram_map: HashMap<i64, i64>,
    mapping_quality_histogram_cache: Vec<i64>,
    mapping_quality_histogram: XYVector,

    // insert size
    insert_size_across_reference: Vec<f64>,
    mean_insert_size: f64,
    p25_insert_size: i32,
    median_insert_size: i32,
    p75_insert_size: i32,
    std_insert_size: f64,
    insert_size_histogram: XYVector,
    insert_size_histogram_map: HashMap<i64, i64>,
    insert_size_histogram_cache: Vec<i64>,
    insert_size_array: Vec<usize>,
}

impl FastBamAnalysis {
    const INITIAL_SIZE: usize = 64;
    const NUM_BINS: usize = 1000;
    pub const CACHE_SIZE: usize = 2000;
    const SMOOTH_DISTANCE: i32 = 0;

    pub fn new(_bam_file: String, qualimap_config: &QualimapConfig) -> Self {
        let _thread_num = qualimap_config.get_thread_num();
        let _window_num = qualimap_config.get_window_num();

        Self {
            bam_file: _bam_file,
            thread_num: _thread_num,
            number_of_windows: _window_num,
            total_num_of_records: 0,

            window_size: 0,
            effective_number_of_window: 0,
            window_starts: vec![],
            window_ends: vec![],
            windows: vec![],

            number_of_problematic_reads: 0,
            number_of_valid_reads: 0,
            number_of_correct_strand_reads: 0,
            sum_coverage_squared: 0,
            sum_coverage: 0,

            reference_file: "".to_string(),
            reference_available_flag: false,
            reference_size: 0,
            number_of_reference_config: 0,
            locator: GenomeLocator::new(),
            acum_read_size: 0,
            max_read_size: 0,
            min_read_size: i32::MAX as usize,
            number_of_secondary_alignments: 0,
            number_of_reads: 0,
            selected_regions_available_flag: false,
            num_supplementary_alignments: 0,
            num_mapped_reads: 0,
            num_paired_reads: 0,
            num_mapped_first_in_pair: 0,
            num_mapped_second_in_pair: 0,
            num_singletons: 0,
            num_marked_duplicates: 0,

            skip_marked_duplicates_flag: false,
            num_of_duplicates_skip: 0,
            skip_detected_duplicates_flag: false,
            num_estimated_duplicate_reads: 0,
            collect_intersecting_paired_end_reads_flag: false,

            // read starts
            read_starts_histogram: ReadStartsHistogram::new(),
            unique_read_starts_histogram: XYVector::new(),
            duplication_rate: 0.0,

            insert_size_array: vec![],

            num_insertions: 0,
            num_read_with_insertion: 0,
            num_deletions: 0,
            num_read_with_deletion: 0,

            edit_distance: 0,
            num_mismatches: 0,
            num_clipped_reads: 0,

            // read content
            reads_a_content: vec![0; Self::INITIAL_SIZE],
            reads_c_content: vec![0; Self::INITIAL_SIZE],
            reads_g_content: vec![0; Self::INITIAL_SIZE],
            reads_t_content: vec![0; Self::INITIAL_SIZE],
            reads_n_content: vec![0; Self::INITIAL_SIZE],
            reads_as_histogram: XYVector::new(),
            reads_cs_histogram: XYVector::new(),
            reads_gs_histogram: XYVector::new(),
            reads_ts_histogram: XYVector::new(),
            reads_ns_histogram: XYVector::new(),
            a_content_across_reference: vec![],
            t_content_across_reference: vec![],
            c_content_across_reference: vec![],
            g_content_across_reference: vec![],
            n_content_across_reference: vec![],
            a_relative_content_across_reference: vec![],
            c_relative_content_across_reference: vec![],
            t_relative_content_across_reference: vec![],
            g_relative_content_across_reference: vec![],
            n_relative_content_across_reference: vec![],

            // gc content
            gc_content_across_reference: vec![],
            gc_relative_content_across_reference: vec![],
            gc_content_histogram: vec![0.0; Self::NUM_BINS + 1],

            // clipping content
            reads_clipping_profile_histogram: XYVector::new(),
            reads_clipping_content: vec![0; Self::INITIAL_SIZE],

            reads_gc_content: vec![],
            sample_count: 0,
            homopolymer_indels: vec![0; 6],

            num_bases: 0,
            num_gc: 0,
            min_homopolymer_size: QualimapConstants::DEFAULT_HOMOPOLYMER_SIZE as usize,
            prev_base_inside_indel_region_flag: false,
            homopolymer_starts_inside_indel_region_flag: false,
            prev_base: 0,
            homopolymer_size: 1,

            num_read_with_correct_insert_sizes: 0,
            acum_insert_size: 0.0,
            number_of_reads_with_start_greater_than_end: 0,

            number_of_sequenced_bases: 0,
            number_of_mapped_bases: 0,
            number_of_as: 0,
            number_of_ts: 0,
            number_of_cs: 0,
            number_of_gs: 0,
            number_of_ns: 0,

            pg_program: "".to_string(),
            pg_command_string: "".to_string(),

            // coverage
            coverage_per_position: vec![],
            coverage_across_reference: vec![],
            std_coverage_across_reference: vec![],
            coverage_per_window: vec![],
            coverage_squared_per_window: vec![],
            coverage_histogram_map: HashMap::with_capacity(
                QualimapConstants::DEFAULT_NUMBER_OF_WINDOWS as usize,
            ),
            coverage_histogram_cache: vec![0; Self::CACHE_SIZE],
            coverage_histogram: XYVector::new(),
            acum_coverage_histogram: XYVector::new(),
            balanced_coverage_histogram: XYVector::new(),
            balanced_coverage_bar_names: HashMap::new(),
            max_coverage_quota: 50,
            coverage_quotes: XYVector::new(),

            // mapping quality
            mapping_quality_per_position: vec![],
            mapping_quality_across_reference: vec![],
            mapping_quality_histogram_map: HashMap::with_capacity(
                QualimapConstants::DEFAULT_NUMBER_OF_WINDOWS as usize,
            ),
            mapping_quality_histogram_cache: vec![0; Self::CACHE_SIZE],
            mapping_quality_histogram: XYVector::new(),

            // insert size
            insert_size_across_reference: vec![],
            mean_insert_size: 0.0,
            p25_insert_size: 0,
            median_insert_size: 0,
            p75_insert_size: 0,
            std_insert_size: 0.0,
            insert_size_histogram: XYVector::new(),
            insert_size_histogram_map: HashMap::with_capacity(
                QualimapConstants::DEFAULT_NUMBER_OF_WINDOWS as usize,
            ),
            insert_size_histogram_cache: vec![0; Self::CACHE_SIZE as usize],
        }
    }

    fn find_index(arr: &Vec<usize>, m: usize) -> Option<usize> {
        let mut left = 0;
        let mut right = arr.len();

        while left < right {
            let mid = left + (right - left) / 2;

            if arr[mid] <= m {
                left = mid + 1;
            } else {
                right = mid;
            }
        }

        if left > 0 && arr[left - 1] <= m {
            Some(left - 1)
        } else {
            None
        }
    }

    pub fn process_sequence(&mut self, record: &Record) {
        let contig_id = record.tid();
        let contig_option = self.locator.get_contig(contig_id as usize);

        // compute absolute position
        let mut absolute_position = -1; // start from 1 instead of 0
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
            self.number_of_valid_reads += 1;

            // accumulate only mapped reads
            if record.is_unmapped() {
                return;
            }

            let mut insert_size = 0;
            if record.is_paired() {
                insert_size = record.insert_size()
            }

            if self.selected_regions_available_flag {
                println!("selected regions flag is true");
                // todo!()
            } else {
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

                // marked duplicates
                if record.is_duplicate() {
                    self.num_marked_duplicates += 1;
                    if self.skip_marked_duplicates_flag {
                        self.num_of_duplicates_skip += 1;
                        return;
                    }
                }

                // estimated duplicates
                if self.update_read_start_histogram_and_judge_is_detected_dup(absolute_position)
                    && self.skip_detected_duplicates_flag
                {
                    self.num_of_duplicates_skip += 1;
                    return;
                }

                if self.collect_intersecting_paired_end_reads_flag {
                    // todo!()
                    // self.bam_stats_collector.collect_paired_read_info(record);
                }

                // update insert size array
                if insert_size > 0 {
                    self.insert_size_array.push(insert_size as usize);
                }
            }

            let read_in_region_flag = true;
            let mut alignment: Vec<char> = vec![];
            if read_in_region_flag {
                // it costs a lot of time, about 0.5~0.6 seconds
                alignment = self.compute_read_alignment_and_collect(record);
            }
            // todo!()
            // else if self.analyze_regions_flag && self.compute_outside_stats_flag {
            //     alignment = self.compute_outside_read_alignment_and_collect(record);
            // } else {
            //     alignment = self.compute_read_alignment(record);
            // }
            if alignment.is_empty() {
                return;
            }

            let mut index =
                Self::find_index(&self.window_starts, absolute_position as usize).unwrap();
            let mut current_window_start = self.window_starts[index];
            let mut current_window_end = self.window_ends[index];
            let mut current_window = &mut self.windows[index];

            if record.is_proper_pair() && insert_size > 0 {
                // update insert size number
                self.num_read_with_correct_insert_sizes += 1;
                self.acum_insert_size += insert_size.abs() as f64;

                current_window.acum_insert_sizes += insert_size.abs() as f64;
                current_window.correct_insert_sizes += 1;
            }

            let mapping_quality = record.mapq();
            let read_absolute_end = absolute_position + alignment.len() as i64 - 1;

            if absolute_position > read_absolute_end {
                self.number_of_reads_with_start_greater_than_end += 1;
            }

            // todo!() consider in region and outsdie region
            let mut pos = -1;
            let mut relative = -1;

            for j in absolute_position..=read_absolute_end {
                // handle if out of current bounds
                if j > current_window_end as i64 {
                    index += 1;
                    current_window_start = self.window_starts[index];
                    current_window_end = self.window_ends[index];
                    current_window = &mut self.windows[index];
                }
                // continue to process
                if j <= current_window_end as i64 {
                    // acum mapped base
                    self.number_of_mapped_bases += 1;
                    current_window.number_of_mapped_bases += 1;

                    pos = j - absolute_position;
                    relative = j - current_window_start as i64;
                    let nucleotide = alignment[pos as usize];

                    if nucleotide != '-' {
                        // base stats
                        self.number_of_sequenced_bases += 1;

                        // mapping quality
                        if mapping_quality != 0 {
                            self.mapping_quality_per_position[index][relative as usize] +=
                                mapping_quality as i64;
                            current_window.acum_mapping_quality += mapping_quality as f64;
                        }
                        // coverage
                        self.coverage_per_position[index][relative as usize] += 1;

                        // ATCG content
                        match nucleotide {
                            'A' => {
                                self.number_of_as += 1;
                                current_window.number_of_as += 1;
                            }
                            'C' => {
                                self.number_of_cs += 1;
                                current_window.number_of_cs += 1;
                            }
                            'T' => {
                                self.number_of_ts += 1;
                                current_window.number_of_ts += 1;
                            }
                            'G' => {
                                self.number_of_gs += 1;
                                current_window.number_of_gs += 1;
                            }
                            'N' => {
                                self.number_of_ns += 1;
                                current_window.number_of_ns += 1;
                            }
                            _ => {}
                        }
                    }
                }
            }
        } else {
            self.number_of_problematic_reads += 1;
        }
    }

    fn run(&mut self) {
        let mut bam_reader = Reader::from_path(&self.bam_file).unwrap();
        let thread_pool = tpool::ThreadPool::new(num_cpus::get() as u32).unwrap();
        bam_reader.set_thread_pool(&thread_pool);

        let mut record = Record::new();
        while let Some(result) = bam_reader.read(&mut record) {
            match result {
                Ok(_) => {
                    self.process_sequence(&record);
                }
                Err(_) => {
                    self.number_of_problematic_reads += 1;
                }
            }
        }

        // save the last time;
        self.save_gc();

        self.finish();
    }

    fn analyze_intervals(&mut self, interval: (usize, usize)) {
        let mut bam_reader = Reader::from_path(&self.bam_file).unwrap();
        let thread_pool = tpool::ThreadPool::new(num_cpus::get() as u32).unwrap();
        bam_reader.set_thread_pool(&thread_pool);

        println!("Processing interval:{}:{}", interval.0, interval.1);

        for (index, result) in bam_reader.records().by_ref().skip(interval.0).enumerate() {
            if index == interval.1 - interval.0 + 1 {
                break;
            }
            match result {
                Ok(record) => {
                    self.process_sequence(&record);
                }
                Err(_) => {
                    self.number_of_problematic_reads += 1;
                }
            }
        }

        // for (index, result) in bam_reader
        //     .records()
        //     .enumerate()
        //     .into_iter()
        //     .take_while(|(x, y)| *x >= interval.0 && *x <= interval.1)
        // {
        //     num += 1;
        //     match result {
        //         Ok(_) => {
        //             self.process_sequence(&record);
        //         }
        //         Err(_) => {
        //             self.number_of_problematic_reads += 1;
        //         }
        //     }
        // }

        // while let Some(result) = bam_reader.read(&mut record) {
        //     index += 1;
        //     if index < interval.0 as i64 || index > interval.1 as i64 {
        //         continue;
        //     }
        //     match result {
        //         Ok(_) => {
        //             self.process_sequence(&record);
        //         }
        //         Err(_) => {
        //             self.number_of_problematic_reads += 1;
        //         }
        //     }
        // }

        self.save_gc();

        // gc may have a lot of problems
        // if interval.1 == self.total_num_of_records - 1 {
        //     self.save_gc();
        // }
    }

    fn load_and_init(&mut self) {
        let mut last_action_done = "Loading sam header...";
        println!("{}", last_action_done);
        let mut bam_reader = Reader::from_path(&self.bam_file).unwrap();
        bam_reader.set_threads(num_cpus::get() as usize);
        let start_time = Instant::now();
        // set total records
        self.total_num_of_records = bam_reader.records().by_ref().count();
        println!("time to calculate total{:?}", Instant::now() - start_time);

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
        self.compute_window_positions(self.window_size);
        self.effective_number_of_window = self.window_starts.len();
    }

    fn simple_load_and_init(&mut self, locater: GenomeLocator) {
        self.locator = locater;

        self.load_reference();

        // init window set
        self.window_size = self.compute_window_size(self.reference_size, self.number_of_windows);
        self.compute_window_positions(self.window_size);
        self.effective_number_of_window = self.window_starts.len();
    }

    fn report_inside_region_info(&self) -> String {
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

    pub fn merge(&mut self, other: &Self) {
        for i in 0..self.effective_number_of_window {
            self.windows[i].merge(&other.windows[i]);
            let window_size = self.windows[i].end - self.windows[i].start + 1;

            for (a, b) in self.mapping_quality_per_position[i]
                .iter_mut()
                .zip(other.mapping_quality_per_position[i].iter())
            {
                *a += b;
            }

            for (a, b) in self.coverage_per_position[i]
                .iter_mut()
                .zip(other.coverage_per_position[i].iter())
            {
                *a += b;
            }
        }

        self.acum_read_size += other.acum_read_size;
        self.max_read_size = self.max_read_size.max(other.max_read_size);
        self.min_read_size = self.min_read_size.min(other.min_read_size);
        self.number_of_secondary_alignments += other.number_of_secondary_alignments;
        self.number_of_reads += other.number_of_reads;
        self.number_of_valid_reads += other.number_of_valid_reads;
        self.num_supplementary_alignments += other.num_supplementary_alignments;
        self.num_mapped_reads += other.num_mapped_reads;
        self.num_paired_reads += other.num_paired_reads;
        self.num_mapped_first_in_pair += other.num_mapped_first_in_pair;
        self.num_mapped_second_in_pair += other.num_mapped_second_in_pair;
        self.num_singletons += other.num_singletons;
        self.num_marked_duplicates += other.num_marked_duplicates;
        self.num_estimated_duplicate_reads += other.num_estimated_duplicate_reads;
        self.num_of_duplicates_skip += other.num_of_duplicates_skip;

        self.num_read_with_correct_insert_sizes += other.num_read_with_correct_insert_sizes;
        self.acum_insert_size += other.acum_insert_size;
        self.number_of_reads_with_start_greater_than_end +=
            other.number_of_reads_with_start_greater_than_end;

        self.number_of_mapped_bases += other.number_of_mapped_bases;
        self.number_of_sequenced_bases += other.number_of_sequenced_bases;
        self.number_of_as += other.number_of_as;
        self.number_of_ts += other.number_of_ts;
        self.number_of_cs += other.number_of_cs;
        self.number_of_gs += other.number_of_gs;
        self.number_of_ns += other.number_of_ns;
        self.number_of_problematic_reads += other.number_of_problematic_reads;

        self.num_insertions += other.num_insertions;
        self.num_read_with_insertion += other.num_read_with_insertion;
        self.num_deletions += other.num_deletions;
        self.num_read_with_deletion += other.num_read_with_deletion;
        self.num_clipped_reads += other.num_clipped_reads;

        self.num_mismatches += other.num_mismatches;
        self.edit_distance += other.edit_distance;

        // todo!() merge insert size array
        // todo!() merge homopolymer
        // todo!() merge gc content
        self.homopolymer_indels[0] += other.homopolymer_indels[0];
        self.homopolymer_indels[1] += other.homopolymer_indels[1];
        self.homopolymer_indels[2] += other.homopolymer_indels[2];
        self.homopolymer_indels[3] += other.homopolymer_indels[3];
        self.homopolymer_indels[4] += other.homopolymer_indels[4];
    }

    pub fn finish(&mut self) {
        let num_homopolymer_indels: usize = self.homopolymer_indels[0..=4].iter().sum();
        self.homopolymer_indels[5] +=
            self.num_insertions + self.num_deletions - num_homopolymer_indels;

        self.finish_metrics_across_reference();

        self.finish_metrics_across_histograms();

        println!("Finish All");
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

        // todo!()
        // if self.collect_intersecting_paired_end_reads_flag {
        //     self.bam_stats_collector.finalize_alignment_info();
        // }

        println!("\nInside of regions...");
        println!("{}", self.report_inside_region_info());

        // todo!()
        // if self.compute_outside_stats_flag {
        //     println!("\nOutside of regions...");
        //     println!("{}", self.outside_bam_stats_collector.report());
        // }
        if self.number_of_reads == 0 {
            panic!("The BAM file is empty or corrupt")
        }

        if self.num_mapped_reads > 0 {
            println!("numberOfMappedBases: {}", self.number_of_mapped_bases);
            println!("referenceSize: {}", self.reference_size);
            println!("numberOfSequencedBases: {}", self.number_of_sequenced_bases);
            println!("numberOfAs: {}", self.number_of_as);
        } else {
            println!("\nWARNING: number of mapped reads equals zero");
            println!(
                "Total number of mapped reads or mapped reads in region equals zero.\n
                        For more details, check the number of Unmapped reads."
            );
        }
    }

    pub fn finish_metrics_across_histograms(&mut self) {
        // update from reference array

        // todo!() region analyze
        // if self.num_selected_regions > 0 {
        //     self.mean_mapping_quality_per_window =
        //         Self::compute_mean_val_from_histogram(&self.mapping_quality_histogram);
        // }

        self.compute_mapping_quality_hisitogram();
        self.compute_insert_size_histogram();
        self.compute_coverage_histogram();
        self.compute_unique_read_starts_histogram();
        self.compute_gc_content_histogram();
        self.compute_reads_content_histogram();
        self.compute_reads_clipping_profile_histogram();
    }

    fn compute_mapping_quality_hisitogram(&mut self) {
        // update cache and map
        self.update_mapping_histogram();

        self.mapping_quality_histogram =
            Self::compute_vector_histogram(&self.mapping_quality_histogram_map);
    }
    fn compute_reads_clipping_profile_histogram(&mut self) {
        if self.max_read_size == 0 || !self.clipping_is_present() {
            return;
        }
        Self::ensure_list_size(&mut self.reads_clipping_content, self.max_read_size);

        let mut total_bases_clipped = 0.0;
        for val in &self.reads_clipping_content {
            total_bases_clipped += *val as f64;
        }

        for pos in 0..self.max_read_size {
            let val =
                (self.reads_clipping_content[pos as usize] as f64 / total_bases_clipped) * 100.0;
            self.reads_clipping_profile_histogram
                .add_item(XYItem::new(pos as f64, val));
        }
    }

    pub fn clipping_is_present(&self) -> bool {
        for val in &self.reads_clipping_content {
            if *val > 0 {
                return true;
            }
        }
        return false;
    }

    // easy
    pub fn compute_reads_content_histogram(&mut self) {
        // let total_size = self.reads_as_data.len()
        //     + self.reads_ts_data.len()
        //     + self.reads_cs_data.len()
        //     + self.reads_gs_data.len()
        //     + self.reads_ns_data.len();

        // make sure that we have enough data
        Self::ensure_list_size(&mut self.reads_a_content, self.max_read_size);
        Self::ensure_list_size(&mut self.reads_t_content, self.max_read_size);
        Self::ensure_list_size(&mut self.reads_c_content, self.max_read_size);
        Self::ensure_list_size(&mut self.reads_g_content, self.max_read_size);
        Self::ensure_list_size(&mut self.reads_n_content, self.max_read_size);

        for i in 0..self.max_read_size {
            let num_a = self.reads_a_content[i as usize] as f64;
            let num_t = self.reads_t_content[i as usize] as f64;
            let num_c = self.reads_c_content[i as usize] as f64;
            let num_g = self.reads_g_content[i as usize] as f64;
            let num_n = self.reads_n_content[i as usize] as f64;
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

    fn ensure_list_size(array: &mut Vec<usize>, expected_size: usize) {
        let size = array.len();
        if size < expected_size {
            for _ in 0..expected_size - size + 1 {
                array.push(0);
            }
        }
    }

    pub fn compute_gc_content_histogram(&mut self) {
        for i in 0..self.reads_gc_content.len() {
            let index = (self.reads_gc_content[i] * Self::NUM_BINS as f32).floor();
            self.gc_content_histogram[index as usize] += 1.0;
        }

        self.sample_count = self.reads_gc_content.len();

        //normalize
        //double normalizer = sampleCount - gcContentHistogram[0];
        for i in 0..Self::NUM_BINS + 1 {
            self.gc_content_histogram[i as usize] /= self.sample_count as f64;
        }

        //smooth
        for i in Self::SMOOTH_DISTANCE..Self::NUM_BINS as i32 - Self::SMOOTH_DISTANCE {
            let mut res = 0.0;
            for j in 0..2 * Self::SMOOTH_DISTANCE + 1 {
                let index = (i + j - Self::SMOOTH_DISTANCE) as usize;
                res += self.gc_content_histogram[index];
            }
            self.gc_content_histogram[i as usize] = res / (Self::SMOOTH_DISTANCE * 2 + 1) as f64;
        }
    }

    // easy
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

    pub fn compute_coverage_histogram(&mut self) {
        // calculate cache and map
        self.update_coverage_histogram();

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

    // easy to do
    pub fn compute_insert_size_histogram(&mut self) {
        if self.insert_size_array.is_empty() {
            return;
        }

        self.insert_size_array.sort();
        let size = self.insert_size_array.len();
        let median_index = size / 2;
        let percentile_25_index = size / 4;
        let percentile_75_index = percentile_25_index * 3;

        self.p25_insert_size = self.insert_size_array[percentile_25_index] as i32;
        self.median_insert_size = self.insert_size_array[median_index] as i32;
        self.p75_insert_size = self.insert_size_array[percentile_75_index] as i32;

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

    fn compute_vector_histogram(map: &HashMap<i64, i64>) -> XYVector {
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

    fn update_coverage_histogram(&mut self) {
        // insert val in cache vector range from 0 to Self::CACHE_SIZE-1 into  map
        for i in 0..Self::CACHE_SIZE as usize {
            let val = self.coverage_histogram_cache[i];
            if val > 0 {
                self.coverage_histogram_map.insert(i as i64, val);
            }
        }
    }

    fn update_mapping_histogram(&mut self) {
        // insert val in cache vector range from 0 to Self::CACHE_SIZE-1 into  map
        for i in 0..Self::CACHE_SIZE as usize {
            let val = self.mapping_quality_histogram_cache[i];
            if val > 0 {
                self.mapping_quality_histogram_map.insert(i as i64, val);
            }
        }
    }

    fn allocate_space_in_window(&mut self, window_num: usize) {
        // coverage
        self.coverage_across_reference = vec![0.0; window_num];
        self.std_coverage_across_reference = vec![0.0; window_num];
        self.coverage_per_window = vec![0; window_num];
        self.coverage_squared_per_window = vec![0; window_num];

        // mapping quality
        self.mapping_quality_across_reference = vec![0.0; window_num];

        // insert size
        self.insert_size_across_reference = vec![0.0; window_num];

        // base content
        self.a_content_across_reference = vec![0.0; window_num];
        self.t_content_across_reference = vec![0.0; window_num];
        self.c_content_across_reference = vec![0.0; window_num];
        self.g_content_across_reference = vec![0.0; window_num];
        self.n_content_across_reference = vec![0.0; window_num];
        self.a_relative_content_across_reference = vec![0.0; window_num];
        self.t_relative_content_across_reference = vec![0.0; window_num];
        self.c_relative_content_across_reference = vec![0.0; window_num];
        self.g_relative_content_across_reference = vec![0.0; window_num];
        self.n_relative_content_across_reference = vec![0.0; window_num];

        // gc content
        self.gc_content_across_reference = vec![0.0; window_num];
        self.gc_relative_content_across_reference = vec![0.0; window_num];
    }

    pub fn finish_metrics_across_reference(&mut self) {
        // todo!() region analyze
        // if self.selected_regions_available_flag {
        //     self.effective_window_length =
        //         self.selected_regions.iter().filter(|x| *x).count() as i64
        // }

        // allocate vertor space to speed up push
        let window_num = self.window_starts.len();
        self.allocate_space_in_window(window_num);

        for i in 0..window_num {
            let mut window = &mut self.windows[i];
            let effective_window_length = window.end - window.start + 1;
            // todo! region
            // if self.selected_regions_available_flag {
            //     self.effective_window_length =
            //         self.selected_regions.iter().filter(|x| *x).count() as i64
            // }

            let mut mean_coverage = 0.0;
            if effective_window_length != 0 {
                mean_coverage =
                    window.number_of_mapped_bases as f64 / effective_window_length as f64;

                self.coverage_across_reference[i] = mean_coverage;
            }

            if window.correct_insert_sizes > 0 {
                let mean_insert_size =
                    window.acum_insert_sizes as f64 / window.correct_insert_sizes as f64;

                self.insert_size_across_reference[i] = mean_insert_size;
            }

            // ACTG absolute content
            if mean_coverage > 0.0 {
                let mean_mapping_quality =
                    window.acum_mapping_quality / window.number_of_mapped_bases as f64;

                self.mapping_quality_across_reference[i] = mean_mapping_quality;

                let mean_a_content = window.number_of_as as f64 / mean_coverage;
                let mean_t_content = window.number_of_ts as f64 / mean_coverage;
                let mean_c_content = window.number_of_cs as f64 / mean_coverage;
                let mean_g_content = window.number_of_gs as f64 / mean_coverage;
                let mean_n_content = window.number_of_ns as f64 / mean_coverage;
                let mean_gc_content = mean_g_content + mean_c_content;

                // ACTG relative content
                let acum_mean_content = mean_a_content
                    + mean_t_content
                    + mean_c_content
                    + mean_g_content
                    + mean_n_content;
                let mean_a_relative_content = (mean_a_content / acum_mean_content) * 100.0;
                let mean_t_relative_content = (mean_t_content / acum_mean_content) * 100.0;
                let mean_c_relative_content = (mean_c_content / acum_mean_content) * 100.0;
                let mean_g_relative_content = (mean_g_content / acum_mean_content) * 100.0;
                let mean_n_relative_content = (mean_n_content / acum_mean_content) * 100.0;
                let mean_gc_relative_content = mean_g_relative_content + mean_c_relative_content;

                self.a_content_across_reference[i] = mean_a_content;
                self.t_content_across_reference[i] = mean_t_content;
                self.c_content_across_reference[i] = mean_c_content;
                self.g_content_across_reference[i] = mean_g_content;
                self.n_content_across_reference[i] = mean_n_content;
                self.gc_content_across_reference[i] = mean_gc_content;

                self.a_relative_content_across_reference[i] = mean_a_relative_content;
                self.t_relative_content_across_reference[i] = mean_t_relative_content;
                self.c_relative_content_across_reference[i] = mean_c_relative_content;
                self.g_relative_content_across_reference[i] = mean_g_relative_content;
                self.n_relative_content_across_reference[i] = mean_n_relative_content;
                self.gc_relative_content_across_reference[i] = mean_gc_relative_content;
            }

            let mut sum_coverage = 0;
            let mut sum_coverage_squared = 0;
            for j in 0..self.coverage_per_position[i].len() {
                let coverage_at_postion = self.coverage_per_position[i][j];
                if coverage_at_postion > 0 {
                    // quality
                    self.mapping_quality_per_position[i][j] =
                        self.mapping_quality_per_position[i][j] / coverage_at_postion as i64;

                    sum_coverage_squared += coverage_at_postion as i64 * coverage_at_postion as i64;
                    sum_coverage += coverage_at_postion as i64;
                } else {
                    // make it invalid for histogram
                    self.mapping_quality_per_position[i][j] = -1;
                }
            }

            // compute std coverageData
            let mean_coverage =
                (sum_coverage as f64 / effective_window_length as f64).floor() as usize;
            let std_coverage = (sum_coverage_squared as f64 / effective_window_length as f64
                - mean_coverage as f64 * mean_coverage as f64)
                .sqrt();

            self.sum_coverage += sum_coverage as usize;
            self.sum_coverage_squared += sum_coverage_squared as usize;
            self.coverage_per_window[i] = sum_coverage as usize;
            self.coverage_squared_per_window[i] = sum_coverage_squared as usize;
            self.std_coverage_across_reference[i] = std_coverage;

            // update map and cache
            for j in 0..self.coverage_per_position[i].len() {
                // todo!() region analyze
                // if window.selected_regions_available_flag && !window.get_region().get(i).unwrap() {
                //     continue;
                // }

                // coverageData
                let coverage = self.coverage_per_position[i][j];
                if coverage < Self::CACHE_SIZE as i64 {
                    self.coverage_histogram_cache[coverage as usize] += 1;
                } else if !self.coverage_histogram_map.contains_key(&coverage) {
                    self.coverage_histogram_map.insert(coverage, 1);
                } else {
                    let value = self.coverage_histogram_map.get(&coverage).unwrap();
                    let update_value = *value + 1;
                    self.coverage_histogram_map.insert(coverage, update_value);
                }

                // quality data
                let quality = self.mapping_quality_per_position[i][j];
                if quality != -1 {
                    if quality < Self::CACHE_SIZE as i64 {
                        self.mapping_quality_histogram_cache[quality as usize] += 1;
                    } else if !self.mapping_quality_histogram_map.contains_key(&quality) {
                        self.mapping_quality_histogram_map.insert(quality, 1);
                    } else {
                        let value = self.mapping_quality_histogram_map.get(&quality).unwrap();
                        let update_value = *value + 1;
                        self.mapping_quality_histogram_map
                            .insert(quality, update_value);
                    }
                }
            }
        }

        println!("hello world");
    }

    fn compute_window_size(&self, reference_size: usize, number_of_windows: usize) -> usize {
        let window_size = (reference_size as f64 / number_of_windows as f64).floor() as usize;
        if (reference_size as f64 / number_of_windows as f64) > window_size as f64 {
            return window_size + 1;
        }
        window_size
    }

    fn compute_window_positions(&mut self, window_size: usize) {
        let contigs = self.locator.get_contigs();

        let mut start_pos = 1;
        let mut i = 0;
        let num_configs = contigs.len();
        while start_pos < self.reference_size {
            self.window_starts.push(start_pos);
            start_pos += window_size;
            while i < num_configs {
                let next_contig_start = contigs[i].end() as usize + 1;
                if start_pos >= next_contig_start {
                    if next_contig_start < self.reference_size && start_pos > next_contig_start {
                        self.window_starts.push(next_contig_start);
                    }
                    i += 1;
                } else {
                    break;
                }
            }
        }

        // calculate window end position and windows
        self.window_ends = vec![0; self.window_starts.len()];
        self.windows = Vec::with_capacity(self.window_starts.len() * 2);
        self.coverage_per_position = vec![vec![]; self.window_starts.len()];
        self.mapping_quality_per_position = vec![vec![]; self.window_starts.len()];
        for i in 0..self.window_starts.len() {
            if i == self.window_starts.len() - 1 {
                self.window_ends[i] = self.reference_size;
            } else {
                self.window_ends[i] = self.window_starts[i + 1] - 1;
            }
            self.windows.push(BamGenomeWindow::new(
                self.window_starts[i],
                self.window_ends[i],
            ));

            // allocate
            self.coverage_per_position[i] =
                vec![0; self.window_ends[i] - self.window_starts[i] + 1];
            self.mapping_quality_per_position[i] =
                vec![0; self.window_ends[i] - self.window_starts[i] + 1];
        }

        println!("finish init windows");
    }

    fn update_read_start_histogram_and_judge_is_detected_dup(&mut self, position: i64) -> bool {
        let duplicate = self.read_starts_histogram.update(position);
        if duplicate {
            self.num_estimated_duplicate_reads += 1;
        }
        duplicate
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
                self.num_insertions += 1;
                if !read_has_insertion {
                    self.num_read_with_insertion += 1;
                }
                read_has_insertion = true;
            } else if c.char() == 'D' {
                self.num_deletions += 1;
                if !read_has_deletion {
                    self.num_read_with_deletion += 1;
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
        self.reset_counters();
        for pos in 0..extended_cigar_vector.len() {
            let cigar_char = extended_cigar_vector[pos];

            if cigar_char == 'M' || cigar_char == '=' {
                let base = read_bases[read_pos];
                self.collect_base(read_pos, base, false);
                read_pos += 1;
                alignment_vector[alignment_pos as usize] = base as char;
                alignment_pos += 1;
            } else if cigar_char == 'I' {
                let base = read_bases[read_pos];
                self.collect_base(read_pos, base, true);
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
                    self.collect_deleted_base(next_base);
                }
                alignment_vector[alignment_pos as usize] = '-';
                alignment_pos += 1;
            } else if cigar_char == 'N' {
                alignment_vector[alignment_pos as usize] = 'N';
                alignment_pos += 1;
            } else if cigar_char == 'S' {
                self.inc_clipping_content(read_pos);
                read_pos += 1;
            } else if cigar_char == 'H' {
                self.inc_clipping_content(read_pos);
            } else if cigar_char == 'P' {
                alignment_vector[alignment_pos as usize] = '-';
                alignment_pos += 1;
            }
        }

        if read_is_clipped {
            self.num_clipped_reads += 1
        }

        let num_mismatches = self.compute_read_mismatches(read);

        self.num_mismatches += num_mismatches;

        match read.aux(b"NM") {
            Ok(value) => {
                if let Aux::U8(v) = value {
                    self.edit_distance += v as usize;
                }
            }
            Err(e) => {}
        }

        alignment_vector
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

    fn save_gc(&mut self) {
        if self.num_gc != 0 {
            let gc_content = self.num_gc as f32 / self.num_bases as f32;
            self.reads_gc_content.push(gc_content);
        }

        self.num_bases = 0;
        self.num_gc = 0;
    }

    fn update_homopolymer_stats(&mut self, base: u8, inside_indel_region: bool) {
        if base == self.prev_base {
            self.homopolymer_size += 1;
        } else {
            // judge boundaries
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

    fn reset_counters(&mut self) {
        self.prev_base = 0;
        self.homopolymer_size = 1;
        self.prev_base_inside_indel_region_flag = false;
        self.homopolymer_starts_inside_indel_region_flag = false;
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

    fn ensure_array_size(array: &mut Vec<usize>, pos: usize) {
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
            self.reference_size = self.locator.get_total_size() as usize;
            self.number_of_reference_config = self.locator.get_contigs().len();
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
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct FastBam {}

impl FastBam {
    pub fn new() -> Self {
        Self {}
    }

    pub fn run(&mut self, bam_file: &str, qualimap_config: &QualimapConfig) {
        let start_run_time = Instant::now();
        let thread_num = qualimap_config.get_thread_num();
        let mut time_to_init = Duration::seconds(0);

        let mut fast_bam_analysis = FastBamAnalysis::new(bam_file.to_string(), qualimap_config);
        let start_time = Instant::now();
        fast_bam_analysis.load_and_init();
        time_to_init += Instant::now() - start_time;
        // fast_bam_analysis.run();
        println!("total_record:{}", fast_bam_analysis.total_num_of_records);
        let intervals =
            Self::divide_interval(fast_bam_analysis.total_num_of_records, thread_num + 1);

        let arc_bam_file = Arc::new(bam_file.to_string());
        let arc_qualimap_config = Arc::new(qualimap_config.clone());
        let arc_intervals = Arc::new(intervals);
        let arc_locator = Arc::new(fast_bam_analysis.locator.clone());

        let mut handles = vec![];
        for i in 1..=fast_bam_analysis.thread_num {
            let bam_file = arc_bam_file.clone();
            let qualimap_config = arc_qualimap_config.clone();
            let intervals = arc_intervals.clone();
            let locator = arc_locator.clone();

            handles.push(thread::spawn(move || {
                let mut fast_bam_analysis_instance =
                    FastBamAnalysis::new(bam_file.as_ref().to_owned(), qualimap_config.as_ref());
                let start_time = Instant::now();

                fast_bam_analysis_instance.simple_load_and_init(locator.as_ref().to_owned());
                time_to_init += Instant::now() - start_time;
                fast_bam_analysis_instance.analyze_intervals(intervals.as_ref()[i]);
                fast_bam_analysis_instance
            }))
        }

        fast_bam_analysis.analyze_intervals(arc_intervals.as_ref()[0]);
        let mut results = vec![];
        for handle in handles {
            let t = handle.join().unwrap();
            results.push(t);
        }

        println!("Time to init: {:?}", time_to_init);
        let start_time_to_merge = Instant::now();
        for i in 0..results.len() {
            fast_bam_analysis.merge(&results[i]);
            println!("merge succeeded");
        }
        let time_to_merge = Instant::now() - start_time_to_merge;
        println!("Time to merge: {:?}", time_to_merge);

        let start_time_to_finish = Instant::now();
        fast_bam_analysis.finish();

        let time_to_finish = Instant::now() - start_time_to_finish;
        println!("Time to finish: {:?}", time_to_finish);

        let time_to_overall_analysis = Instant::now() - start_run_time;
        println!("Time to overall analysis: {:?}", time_to_overall_analysis);

        println!(
            "Time to process sequence: {:?}",
            time_to_overall_analysis - time_to_finish - time_to_merge - time_to_init
        );
    }

    fn divide_interval(m: usize, n: usize) -> Vec<(usize, usize)> {
        let interval_size = m / n;
        let mut remainder = m % n;
        let mut intervals = Vec::with_capacity(n as usize);

        let mut start = 0;
        for _ in 0..n {
            let end = start + interval_size + if remainder > 0 { 1 } else { 0 } - 1;
            intervals.push((start, end));
            start = end + 1;
            remainder = remainder.saturating_sub(1);
        }

        intervals
    }
}
