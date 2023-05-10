use crate::qc::bam::qualimap::Qualimap;
use serde::{Deserialize, Serialize};

#[allow(dead_code)]
pub struct Constants {}
impl Constants {
    pub const DEFAULT_NUMBER_OF_WINDOWS: i32 = 400;
    pub const DEFAULT_CHUNK_SIZE: i32 = 1000;
    pub const DEFAULT_HOMOPOLYMER_SIZE: i32 = 3;
    pub const DEFAULT_STABLIZED_WINDOW_PROPORTION: i32 = 500;

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
    thread_num: usize,
    window_num: usize,
    bunch_size: usize,
}

impl QualimapConfig {
    pub fn new(_thread_num: usize, _window_num: usize, _bunch_size: usize) -> Self {
        Self {
            thread_num: _thread_num,
            window_num: _window_num,
            bunch_size: _bunch_size,
        }
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
            qc.run(bam_path.to_string());
            result.qualimap = Some(qc);
        }
        if which == "fastqscreen" || which == "all" {
            // run fastqscreen
            // set results
        }

        result
    }
}
