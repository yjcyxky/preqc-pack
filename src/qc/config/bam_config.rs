use crate::qc::bam::qualimap::Qualimap;
use serde::{Deserialize, Serialize};

#[allow(dead_code)]
pub struct Constants {}
impl Constants {
    pub const DEFAULT_NUMBER_OF_WINDOWS: i32 = 400;
    pub const DEFAULT_CHUNK_SIZE: i32 = 1000;
    pub const DEFAULT_HOMOPOLYMER_SIZE: i32 = 3;
    pub const DEFAULT_STABLIZED_WINDOW_PROPORTION: i32 = 500;
    pub const DEFAULT_LIB_PROTOCOL: &'static str = "non-strand-specific";
    pub const DEFAULT_SKIP_DUPLICATE_MODE: &'static str = "flagged";

    pub const GRAPHIC_TO_SAVE_WIDTH: i32 = 1024;
    pub const GRAPHIC_TO_SAVE_HEIGHT: i32 = 768;

    pub const SAM_FLAG_SUPP_ALIGNMENT: i32 = 0x800;

    /** Path to locate the images when the application is running in a jar file */
    pub const PATH_IMAGES: &str = "/org/bioinfo/ngs/qc/qualimap/gui/images/";

    /** Path to locate the resources of the application */
    pub const PATH_RESOURCES: &str = "/org/bioinfo/ngs/qc/qualimap/";

    //******************************************************************************************
    //******************************* FILE EXTENSION CONSTANTS *********************************
    //******************************************************************************************
    // Extension for the data Input File
    pub const FILE_EXTENSION_BAM: &'static str = "BAM";
    pub const FILE_EXTENSION_SAM: &'static str = "SAM";

    // Extension for the Region Input File
    // pub const FILE_EXTENSION_REGION: HashMap<&'static str, &'static str> = hashmap! {
    //     "GFF" => "GFF",
    //     "GTF" => "GTF",
    //     "BED" => "BED",
    // };

    // Extension for the PDF File
    pub const FILE_EXTENSION_PDF_FILE: &'static str = "PDF";

    //******************************************************************************************
    //******************************* GRAPHICS NAMES CONSTANTS *********************************
    //******************************************************************************************

    pub const PLOT_TITLE_COVERAGE_ACROSS_REFERENCE: &'static str = "Coverage Across Reference";
    pub const PLOT_TITLE_COVERAGE_HISTOGRAM: &'static str = "Coverage Histogram";
    pub const PLOT_TITLE_COVERAGE_HISTOGRAM_0_50: &'static str = "Coverage Histogram (0-50X)";
    pub const PLOT_TITLE_MAPPING_QUALITY_ACROSS_REFERENCE: &'static str =
        "Mapping Quality Across Reference";
    pub const PLOT_TITLE_MAPPING_QUALITY_HISTOGRAM: &'static str = "Mapping Quality Histogram";
    pub const PLOT_TITLE_INSERT_SIZE_ACROSS_REFERENCE: &'static str =
        "Insert Size Across Reference";
    pub const PLOT_TITLE_INSERT_SIZE_HISTOGRAM: &'static str = "Insert Size Histogram";
    pub const PLOT_TITLE_READS_NUCLEOTIDE_CONTENT: &'static str = "Mapped Reads Nucleotide Content";
    pub const PLOT_TITLE_READS_CLIPPING_PROFILE: &'static str = "Mapped Reads Clipping Profile";
    pub const PLOT_TITLE_GENOME_FRACTION_COVERAGE: &'static str = "Genome Fraction Coverage";
    pub const PLOT_TITLE_READS_GC_CONTENT: &'static str = "Mapped Reads GC-content Distribution";
    pub const PLOT_TITLE_DUPLICATION_RATE_HISTOGRAM: &'static str = "Duplication Rate Histogram";
    pub const PLOT_TITLE_HOMOPOLYMER_INDELS: &'static str = "Homopolymer Indels";

    //******************************************************************************************
    //*********************************** TYPES OF SPECIES *************************************
    //******************************************************************************************
    pub const TYPE_COMBO_SPECIES_HUMAN: &'static str = "HUMAN.ENS68";
    pub const TYPE_COMBO_SPECIES_MOUSE: &'static str = "MOUSE.ENS68";

    //******************************************************************************************
    //*********************************** FILES OF SPECIES *************************************
    //******************************************************************************************
    pub const FILE_SPECIES_INFO_HUMAN: &'static str = "human.64.genes.biotypes.txt";
    pub const FILE_SPECIES_GROUPS_HUMAN: &'static str = "human.biotypes.groups.txt";
    pub const FILE_SPECIES_INFO_MOUSE: &'static str = "mouse.64.genes.biotypes.txt";
    pub const FILE_SPECIES_GROUPS_MOUSE: &'static str = "mouse.biotypes.groups.txt";
    pub const FILE_SPECIES_INFO_HUMAN_ENS68: &'static str = "human.ens68.txt";
    pub const FILE_SPECIES_INFO_MOUSE_ENS68: &'static str = "mouse.ens68.txt";
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct QualimapConfig {
    // advanced options
    thread_num: usize,
    window_num: usize,
    bunch_size: usize,
    min_homopolymer_size: usize,
    // region analyze
    feature_file: String,
    outside_stats_flag: bool,
    lib_protocol: String,
    // overlap detection
    collect_overlap_flag: bool,
    // skip duplicates option
    skip_dup_flag: bool,
    skip_dup_mode: String,
    // reference gc content
    gc_genome: String,
}

impl QualimapConfig {
    pub fn new() -> Self {
        Self {
            // advanced options
            thread_num: 1,
            window_num: Constants::DEFAULT_NUMBER_OF_WINDOWS as usize,
            bunch_size: Constants::DEFAULT_CHUNK_SIZE as usize,
            min_homopolymer_size: Constants::DEFAULT_HOMOPOLYMER_SIZE as usize,
            // region analyze
            feature_file: "".to_string(),
            outside_stats_flag: false,
            lib_protocol: Constants::DEFAULT_LIB_PROTOCOL.to_string(),
            // overlap detection
            collect_overlap_flag: false,
            // skip duplicates option
            skip_dup_flag: false,
            skip_dup_mode: Constants::DEFAULT_SKIP_DUPLICATE_MODE.to_string(),
            // reference gc content
            gc_genome: "".to_string(),
        }
    }

    pub fn set_thread_num(&mut self, nthreads: usize) {
        self.thread_num = nthreads;
    }

    pub fn get_thread_num(&self) -> usize {
        self.thread_num
    }

    pub fn set_window_num(&mut self, nwindows: usize) {
        self.window_num = nwindows;
    }

    pub fn get_window_num(&self) -> usize {
        self.window_num
    }

    pub fn set_bunch_size(&mut self, bunch_size: usize) {
        self.bunch_size = bunch_size;
    }

    pub fn get_bunch_size(&self) -> usize {
        self.bunch_size
    }

    pub fn set_min_homopolymer_size(&mut self, min_homopoly: usize) {
        self.min_homopolymer_size = min_homopoly;
    }

    pub fn get_min_homopolymer_size(&self) -> usize {
        self.min_homopolymer_size
    }

    pub fn set_feature_file(&mut self, feature_file: String) {
        self.feature_file = feature_file;
    }

    pub fn get_feature_file(&self) -> String {
        self.feature_file.clone()
    }

    pub fn set_outside_region_analyze_flag(&mut self, flag: bool) {
        self.outside_stats_flag = flag;
    }

    pub fn get_outside_region_analyze_flag(&self) -> bool {
        self.outside_stats_flag
    }

    pub fn set_lib_protocol(&mut self, protocol: String) {
        self.lib_protocol = protocol.to_string();
    }

    pub fn get_lib_protocol(&self) -> &str {
        &self.lib_protocol
    }

    pub fn set_collect_overlap_flag(&mut self, flag: bool) {
        self.collect_overlap_flag = flag;
    }

    pub fn get_collect_overlap_flag(&self) -> bool {
        self.collect_overlap_flag
    }

    pub fn set_skip_duplicate_flag(&mut self, flag: bool) {
        self.skip_dup_flag = flag;
    }

    pub fn get_skip_duplicate_flag(&self) -> bool {
        self.skip_dup_flag
    }

    pub fn set_skip_duplicate_mode(&mut self, mode: String) {
        self.skip_dup_mode = mode;
    }

    pub fn get_skip_duplicate_mode(&self) -> &str {
        &self.skip_dup_mode
    }

    pub fn set_gc_genome(&mut self, genome: String) {
        self.gc_genome = genome;
    }

    pub fn get_gc_genome(&self) -> String {
        self.gc_genome.clone()
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct QCResults {
    qualimap: Option<Qualimap>,
    fastqscreen: Option<String>,
}

impl QCResults {
    pub fn new() -> Self {
        Self {
            qualimap: None,
            fastqscreen: None,
        }
    }
    pub fn run_qc(bam_path: &str, which: &str, qualimap_config: &QualimapConfig) -> QCResults {
        let mut result = QCResults::new();

        let mut qc = Qualimap::new();
        if which == "qualimap" || which == "all" {
            qc.run(bam_path, qualimap_config);
            result.qualimap = Some(qc);
        }
        if which == "fastqscreen" || which == "all" {
            // run fastqscreen
            // set results
        }

        result
    }
}
