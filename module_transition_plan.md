# Module-by-Module Transition Plan

## Module 1: `args.py` → `cli/args.rs`

### Current Python Implementation
- Uses custom `RichParser` extending `argparse.ArgumentParser`
- Custom `FlexibleArgFormatter` for help formatting
- File extension and existence validation
- Default thread count detection

### Rust Implementation Strategy

#### Crate: `clap` (v4.x with derive macros)

```rust
use clap::{Parser, ValueEnum};
use std::path::PathBuf;

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum AmpliconType {
    EndToEnd,
    EndToMid,
    Fragmented,
}

#[derive(Parser, Debug)]
#[command(name = "ampligone")]
#[command(author, version, about = "Accurate primer removal from NGS reads")]
pub struct Args {
    /// Input file (FASTQ, FASTQ.GZ, or BAM)
    #[arg(short, long, value_name = "FILE")]
    pub input: PathBuf,

    /// Second input file for paired-end reads (R2)
    #[arg(long, value_name = "FILE")]
    pub input2: Option<PathBuf>,

    /// Output file (FASTQ or FASTQ.GZ)
    #[arg(short, long, value_name = "FILE")]
    pub output: PathBuf,

    /// Second output file for paired-end reads (R2)
    #[arg(long, value_name = "FILE")]
    pub output2: Option<PathBuf>,

    /// Reference genome in FASTA format
    #[arg(short = 'r', long, value_name = "FILE")]
    pub reference: PathBuf,

    /// Primer file (FASTA or BED format)
    #[arg(short = 'p', long, value_name = "FILE")]
    pub primers: PathBuf,

    /// Amplicon type
    #[arg(short = 'a', long, value_enum, default_value = "end-to-end")]
    pub amplicon_type: AmpliconType,

    /// Number of threads
    #[arg(short, long, default_value_t = num_cpus::get().min(2))]
    pub threads: usize,

    /// Primer mismatch error rate (0.0-1.0)
    #[arg(short, long, default_value = "0.1")]
    pub error_rate: f64,

    /// Fragment lookaround size
    #[arg(long, default_value = "10")]
    pub fragment_lookaround_size: usize,

    /// Export found primers to BED file
    #[arg(long, value_name = "FILE")]
    pub export_primers: Option<PathBuf>,

    /// Enable virtual primer binding
    #[arg(long, default_value = "true")]
    pub virtual_primers: bool,

    /// Quiet mode (suppress progress output)
    #[arg(short, long)]
    pub quiet: bool,

    /// Verbose/debug mode
    #[arg(short, long)]
    pub verbose: bool,
}

impl Args {
    pub fn validate(&self) -> Result<(), String> {
        // Validate file extensions
        self.validate_input_extension()?;
        self.validate_output_extension()?;
        self.validate_files_exist()?;
        self.validate_paired_end_consistency()?;
        Ok(())
    }

    fn validate_paired_end_consistency(&self) -> Result<(), String> {
        match (&self.input2, &self.output2) {
            (Some(_), None) | (None, Some(_)) => {
                Err("Both --input2 and --output2 must be provided for paired-end mode".into())
            }
            _ => Ok(())
        }
    }
}
```

### Migration Checklist
- [ ] Port all argument definitions
- [ ] Implement file extension validation
- [ ] Implement file existence checks
- [ ] Add thread count auto-detection
- [ ] Add paired-end argument validation
- [ ] Port help text and descriptions
- [ ] Add shell completion generation

---

## Module 2: `io_ops.py` → `io/` module

### Current Python Implementation

#### `SequenceReads` class
- Reads FASTQ or BAM files
- Stores tuples of (name, sequence, quality)
- Creates pandas DataFrame for processing
- Handles gzip decompression
- Implements `_flip_strand` for reverse complement

#### `read_bed` function
- Reads BED files with pandas
- Filters browser/track lines

#### `write_output` function
- Writes FASTQ output
- Supports gzip compression with pgzip

### Rust Implementation Strategy

#### File: `io/fastq.rs`

```rust
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct FastqRecord {
    pub name: String,
    pub sequence: Vec<u8>,
    pub quality: Vec<u8>,
}

impl FastqRecord {
    pub fn reverse_complement(&self) -> Self {
        let rc_seq = reverse_complement(&self.sequence);
        let rev_qual: Vec<u8> = self.quality.iter().rev().copied().collect();
        Self {
            name: self.name.clone(),
            sequence: rc_seq,
            quality: rev_qual,
        }
    }
}

pub struct FastqReader<R: BufRead> {
    reader: R,
    buffer: String,
}

impl<R: BufRead> FastqReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            buffer: String::with_capacity(1024),
        }
    }

    pub fn from_path<P: AsRef<Path>>(path: P) -> std::io::Result<FastqReader<Box<dyn BufRead>>> {
        let path = path.as_ref();
        let file = File::open(path)?;
        
        let reader: Box<dyn BufRead> = if path.extension().map_or(false, |e| e == "gz") {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };
        
        Ok(FastqReader::new(reader))
    }
}

impl<R: BufRead> Iterator for FastqReader<R> {
    type Item = std::io::Result<FastqRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        // Read name line
        self.buffer.clear();
        match self.reader.read_line(&mut self.buffer) {
            Ok(0) => return None,
            Ok(_) => {}
            Err(e) => return Some(Err(e)),
        }
        
        let name = self.buffer.trim_start_matches('@')
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();
        
        // Read sequence
        self.buffer.clear();
        if let Err(e) = self.reader.read_line(&mut self.buffer) {
            return Some(Err(e));
        }
        let sequence = self.buffer.trim().as_bytes().to_vec();
        
        // Skip + line
        self.buffer.clear();
        if let Err(e) = self.reader.read_line(&mut self.buffer) {
            return Some(Err(e));
        }
        
        // Read quality
        self.buffer.clear();
        if let Err(e) = self.reader.read_line(&mut self.buffer) {
            return Some(Err(e));
        }
        let quality = self.buffer.trim().as_bytes().to_vec();
        
        Some(Ok(FastqRecord { name, sequence, quality }))
    }
}

pub struct FastqWriter<W: Write> {
    writer: W,
}

impl<W: Write> FastqWriter<W> {
    pub fn new(writer: W) -> Self {
        Self { writer }
    }

    pub fn from_path<P: AsRef<Path>>(path: P, threads: usize) -> std::io::Result<FastqWriter<Box<dyn Write>>> {
        let path = path.as_ref();
        let file = File::create(path)?;
        
        let writer: Box<dyn Write> = if path.extension().map_or(false, |e| e == "gz") {
            // Use pigz-like parallel compression
            Box::new(GzEncoder::new(BufWriter::new(file), flate2::Compression::default()))
        } else {
            Box::new(BufWriter::new(file))
        };
        
        Ok(FastqWriter::new(writer))
    }

    pub fn write_record(&mut self, record: &FastqRecord) -> std::io::Result<()> {
        writeln!(self.writer, "@{}", record.name)?;
        self.writer.write_all(&record.sequence)?;
        writeln!(self.writer)?;
        writeln!(self.writer, "+")?;
        self.writer.write_all(&record.quality)?;
        writeln!(self.writer)?;
        Ok(())
    }
}
```

#### File: `io/bam.rs`

```rust
use rust_htslib::bam::{self, Read, Record};
use std::path::Path;

pub struct BamReader {
    reader: bam::Reader,
}

impl BamReader {
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, rust_htslib::errors::Error> {
        let reader = bam::Reader::from_path(path)?;
        Ok(Self { reader })
    }
}

impl Iterator for BamReader {
    type Item = Result<FastqRecord, rust_htslib::errors::Error>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = Record::new();
        match self.reader.read(&mut record) {
            Some(Ok(())) => {
                if record.is_unmapped() {
                    return self.next(); // Skip unmapped reads
                }

                let name = String::from_utf8_lossy(record.qname()).to_string();
                let mut sequence = record.seq().as_bytes();
                let mut quality: Vec<u8> = record.qual().iter().map(|q| q + 33).collect();

                // Reverse complement if on reverse strand
                if record.is_reverse() {
                    sequence = reverse_complement(&sequence);
                    quality.reverse();
                }

                Some(Ok(FastqRecord {
                    name,
                    sequence,
                    quality,
                }))
            }
            Some(Err(e)) => Some(Err(e)),
            None => None,
        }
    }
}
```

#### File: `io/bed.rs`

```rust
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct BedRecord {
    pub reference: String,
    pub start: u64,
    pub end: u64,
    pub name: String,
    pub score: String,
    pub strand: Strand,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

impl Strand {
    pub fn from_char(c: char) -> Option<Self> {
        match c {
            '+' => Some(Strand::Forward),
            '-' => Some(Strand::Reverse),
            _ => None,
        }
    }

    pub fn to_char(self) -> char {
        match self {
            Strand::Forward => '+',
            Strand::Reverse => '-',
        }
    }
}

pub fn read_bed<P: AsRef<Path>>(path: P) -> std::io::Result<Vec<BedRecord>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut records = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        
        // Skip browser/track lines
        if trimmed.starts_with("browser ") || trimmed.starts_with("track ") {
            continue;
        }
        
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 6 {
            continue;
        }

        let record = BedRecord {
            reference: fields[0].to_string(),
            start: fields[1].parse().unwrap_or(0),
            end: fields[2].parse().unwrap_or(0),
            name: fields[3].to_string(),
            score: fields[4].to_string(),
            strand: Strand::from_char(fields[5].chars().next().unwrap_or('+')).unwrap_or(Strand::Forward),
        };
        records.push(record);
    }

    Ok(records)
}

pub fn write_bed<P: AsRef<Path>>(path: P, records: &[BedRecord]) -> std::io::Result<()> {
    let mut file = File::create(path)?;
    for record in records {
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}",
            record.reference,
            record.start,
            record.end,
            record.name,
            record.score,
            record.strand.to_char()
        )?;
    }
    Ok(())
}
```

### Migration Checklist
- [ ] Implement FastqReader with iterator pattern
- [ ] Implement FastqWriter with buffered writing
- [ ] Add gzip support (flate2 or gzp for parallel)
- [ ] Implement BAM reading with rust-htslib
- [ ] Implement BED reading/writing
- [ ] Add FASTA parsing for reference and primers
- [ ] Implement reverse complement function with SIMD
- [ ] Add comprehensive error handling

---

## Module 3: `fasta2bed.py` → `primer/` module

### Current Python Implementation

#### Key Functions
- `find_ambiguous_options(seq)`: Expands IUPAC ambiguous nucleotides
- `parse_cigar_obj(cig_obj)`: Parses parasail CIGAR objects
- `count_cigar_information(cigar)`: Counts matches, mismatches, insertions, deletions
- `get_coords(seq, ref_seq, err_rate)`: Finds primer coordinates via semi-global alignment
- `find_or_read_primers(primerfile, referencefile, err_rate)`: Main primer finding function
- `coord_list_gen(...)`: Generator for primer coordinates

### Rust Implementation Strategy

#### File: `primer/search.rs`

```rust
use parasailors::{Matrix, Profile, semi_global_alignment};
use bio::alphabets::dna::iupac_alphabet;
use std::collections::HashMap;

/// IUPAC ambiguity code expansion
pub fn expand_ambiguous(seq: &[u8]) -> Vec<Vec<u8>> {
    static IUPAC: phf::Map<u8, &[u8]> = phf::phf_map! {
        b'A' => b"A",
        b'C' => b"C",
        b'G' => b"G",
        b'T' => b"T",
        b'R' => b"AG",
        b'Y' => b"CT",
        b'S' => b"GC",
        b'W' => b"AT",
        b'K' => b"GT",
        b'M' => b"AC",
        b'B' => b"CGT",
        b'D' => b"AGT",
        b'H' => b"ACT",
        b'V' => b"ACG",
        b'N' => b"ACGT",
    };

    let mut results = vec![Vec::new()];
    
    for &base in seq {
        let base_upper = base.to_ascii_uppercase();
        let options = IUPAC.get(&base_upper).copied().unwrap_or(&[base_upper]);
        
        let mut new_results = Vec::with_capacity(results.len() * options.len());
        for result in &results {
            for &option in options.iter() {
                let mut new_seq = result.clone();
                new_seq.push(option);
                new_results.push(new_seq);
            }
        }
        results = new_results;
    }
    
    results
}

#[derive(Debug, Clone)]
pub struct PrimerCoordinates {
    pub reference: String,
    pub start: u64,
    pub end: u64,
    pub name: String,
    pub score: i32,
    pub strand: Strand,
    pub sequence: Vec<u8>,
    pub revcomp: Vec<u8>,
}

pub struct PrimerSearcher {
    error_rate: f64,
    scoring_matrix: Matrix,
    gap_open: i32,
    gap_extend: i32,
}

impl PrimerSearcher {
    pub fn new(error_rate: f64) -> Self {
        Self {
            error_rate,
            scoring_matrix: Matrix::create("ACGT", 2, -1), // nuc44 equivalent
            gap_open: 8,
            gap_extend: 30,
        }
    }

    pub fn find_primer_coordinates(
        &self,
        primer_seq: &[u8],
        reference_seq: &[u8],
    ) -> Option<(u64, u64, i32, i32)> {
        let max_errors = (primer_seq.len() as f64 * self.error_rate).floor() as usize;
        let options = expand_ambiguous(primer_seq);
        
        let mut best_result: Option<(Vec<u8>, u64, u64, i32, i32)> = None;
        
        for option in options {
            let profile = Profile::new(&option, &self.scoring_matrix);
            let result = semi_global_alignment(
                &profile,
                reference_seq,
                self.gap_open,
                self.gap_extend,
            );
            
            let (cigar, cleaned_cigar) = parse_cigar(&result.cigar);
            let errors = count_errors(&cleaned_cigar);
            
            if errors <= max_errors {
                let score = result.score;
                if best_result.is_none() || score > best_result.as_ref().unwrap().3 {
                    best_result = Some((
                        option,
                        result.ref_begin as u64,
                        result.ref_end as u64 + 1,
                        score,
                        calculate_percentage(option.len(), score),
                    ));
                }
            }
        }
        
        best_result.map(|(_, start, end, score, pct)| (start, end, score, pct))
    }

    pub fn search_all_primers(
        &self,
        primers: &[FastaRecord],
        references: &[FastaRecord],
    ) -> Vec<PrimerCoordinates> {
        let mut results = Vec::new();
        
        for reference in references {
            for primer in primers {
                let strand = detect_strand(&primer.name)?;
                let revcomp = reverse_complement(&primer.sequence);
                
                // Try forward orientation
                let fw_coords = self.find_primer_coordinates(&primer.sequence, &reference.sequence);
                let rv_coords = self.find_primer_coordinates(&revcomp, &reference.sequence);
                
                let best = choose_best_fitting(fw_coords, rv_coords)?;
                
                results.push(PrimerCoordinates {
                    reference: reference.name.clone(),
                    start: best.0,
                    end: best.1,
                    name: primer.name.clone(),
                    score: best.2,
                    strand,
                    sequence: primer.sequence.clone(),
                    revcomp,
                });
            }
        }
        
        results
    }
}

fn detect_strand(primer_name: &str) -> Option<Strand> {
    let name_upper = primer_name.to_uppercase();
    let forward_keys = ["LEFT", "PLUS", "POSITIVE", "FORWARD"];
    let reverse_keys = ["RIGHT", "MINUS", "NEGATIVE", "REVERSE"];
    
    for key in forward_keys {
        if name_upper.contains(key) {
            return Some(Strand::Forward);
        }
    }
    for key in reverse_keys {
        if name_upper.contains(key) {
            return Some(Strand::Reverse);
        }
    }
    None
}
```

#### File: `primer/index.rs`

```rust
use std::collections::{HashMap, HashSet};
use std::ops::Range;

/// Index structure for fast primer coordinate lookup
#[derive(Debug, Default)]
pub struct PrimerIndex {
    forward: HashMap<String, HashSet<u64>>,
    reverse: HashMap<String, HashSet<u64>>,
}

impl PrimerIndex {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_coordinates(
        coords: &[PrimerCoordinates],
        bind_virtual_primers: bool,
    ) -> Self {
        let mut index = Self::new();
        
        for coord in coords {
            let range = if bind_virtual_primers {
                // Extend range to include nearby primers
                expand_range_for_virtual_primers(coord, coords)
            } else {
                coord.start..coord.end
            };
            
            let target = match coord.strand {
                Strand::Forward => &mut index.forward,
                Strand::Reverse => &mut index.reverse,
            };
            
            target.entry(coord.reference.clone())
                .or_default()
                .extend(range);
        }
        
        index
    }

    pub fn get_forward(&self, reference: &str) -> Option<&HashSet<u64>> {
        self.forward.get(reference)
    }

    pub fn get_reverse(&self, reference: &str) -> Option<&HashSet<u64>> {
        self.reverse.get(reference)
    }

    pub fn contains_forward(&self, reference: &str, position: u64) -> bool {
        self.forward.get(reference)
            .map_or(false, |set| set.contains(&position))
    }

    pub fn contains_reverse(&self, reference: &str, position: u64) -> bool {
        self.reverse.get(reference)
            .map_or(false, |set| set.contains(&position))
    }
}

fn expand_range_for_virtual_primers(
    coord: &PrimerCoordinates,
    all_coords: &[PrimerCoordinates],
) -> Range<u64> {
    let length = coord.end - coord.start;
    let mut start = coord.start;
    let mut end = coord.end;
    
    for other in all_coords {
        if other.reference == coord.reference 
            && other.strand == coord.strand
            && other.start <= end + length
            && other.end >= start.saturating_sub(length)
        {
            start = start.min(other.start);
            end = end.max(other.end);
        }
    }
    
    (start + 1)..end  // Match Python's range(start + 1, end)
}
```

### Migration Checklist
- [ ] Implement IUPAC ambiguity expansion
- [ ] Port CIGAR parsing logic
- [ ] Implement semi-global alignment with parasailors or rust-bio
- [ ] Create PrimerCoordinates struct
- [ ] Implement PrimerIndex with virtual primer support
- [ ] Add strand detection from primer names
- [ ] Implement coordinate comparison and best-fit selection

---

## Module 4: `alignmentpreset.py` → `alignment/preset.rs`

### Current Python Implementation

The module automatically determines the appropriate minimap2 preset based on read characteristics:
- Analyzes read lengths and quality scores
- Determines sequence variance and stability
- Returns `sr` (short read), `map-ont` (ONT), or other presets

### Rust Implementation Strategy

```rust
use rayon::prelude::*;
use statistical::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentPreset {
    ShortRead,  // sr
    MapOnt,     // map-ont
    MapPb,      // map-pb (PacBio)
}

impl AlignmentPreset {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::ShortRead => "sr",
            Self::MapOnt => "map-ont",
            Self::MapPb => "map-pb",
        }
    }
}

pub struct PresetDetector {
    sample_size: usize,
}

impl PresetDetector {
    pub fn new(sample_size: usize) -> Self {
        Self { sample_size }
    }

    pub fn detect(&self, reads: &[FastqRecord]) -> AlignmentPreset {
        let sample: Vec<_> = reads.iter()
            .take(self.sample_size)
            .collect();
        
        if sample.is_empty() {
            return AlignmentPreset::ShortRead;
        }

        let stats = self.calculate_statistics(&sample);
        self.determine_preset(&stats)
    }

    fn calculate_statistics(&self, reads: &[&FastqRecord]) -> ReadStatistics {
        // Parallel calculation using rayon
        let lengths: Vec<f64> = reads.par_iter()
            .map(|r| r.sequence.len() as f64)
            .collect();
        
        let qualities: Vec<f64> = reads.par_iter()
            .flat_map(|r| r.quality.iter().map(|&q| (q - 33) as f64))
            .collect();
        
        ReadStatistics {
            avg_length: mean(&lengths),
            length_stddev: standard_deviation(&lengths, None),
            avg_quality: mean(&qualities),
            quality_stddev: standard_deviation(&qualities, None),
            min_length: lengths.iter().cloned().fold(f64::INFINITY, f64::min) as usize,
            max_length: lengths.iter().cloned().fold(0.0, f64::max) as usize,
        }
    }

    fn determine_preset(&self, stats: &ReadStatistics) -> AlignmentPreset {
        let length_range = stats.max_length - stats.min_length;
        let quality_range = (stats.quality_stddev * 2.0) as usize;
        
        let stability = self.calculate_stability(stats);
        
        if stability > 97.0 {
            if stats.avg_length < 500.0 {
                AlignmentPreset::ShortRead
            } else {
                // Long Illumina reads (MiSeq)
                AlignmentPreset::ShortRead
            }
        } else if stats.avg_quality < 15.0 {
            AlignmentPreset::MapOnt
        } else {
            AlignmentPreset::MapOnt
        }
    }

    fn calculate_stability(&self, stats: &ReadStatistics) -> f64 {
        let length_variance = stats.length_stddev / stats.avg_length * 100.0;
        let quality_variance = stats.quality_stddev / stats.avg_quality * 100.0;
        
        100.0 - (length_variance + quality_variance) / 2.0
    }
}

#[derive(Debug)]
struct ReadStatistics {
    avg_length: f64,
    length_stddev: f64,
    avg_quality: f64,
    quality_stddev: f64,
    min_length: usize,
    max_length: usize,
}
```

### Migration Checklist
- [ ] Implement read statistics calculation
- [ ] Port preset determination logic
- [ ] Add parallel statistics computation with rayon
- [ ] Implement stability calculation
- [ ] Add unit tests for various read types

---

## Module 5: `cut_reads.py` → `cutting/` module

### Current Python Implementation

This is the core module containing:
- `CuttingParameters` dataclass
- `Read` dataclass
- `cut_read()` function for single read cutting
- `cut_reads()` function for batch processing with minimap2 alignment
- Position checking functions (`position_in_or_before_primer`, `position_in_or_after_primer`)

### Rust Implementation Strategy

#### File: `cutting/types.rs`

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CutDirection {
    Forward,  // 1 in Python
    Reverse,  // -1 in Python
}

#[derive(Debug, Clone)]
pub struct Read {
    pub name: String,
    pub sequence: Vec<u8>,
    pub quality: Vec<u8>,
}

impl Read {
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }
}

#[derive(Debug)]
pub struct CuttingParameters<'a> {
    pub primer_positions: &'a HashSet<u64>,
    pub position_on_reference: u64,
    pub cut_direction: CutDirection,
    pub read_direction: i32,  // 1 for forward, -1 for reverse strand
    pub cigar: Vec<CigarOp>,
    pub query_start: usize,
    pub query_end: usize,
    pub fragment_lookaround_size: u64,
}

#[derive(Debug, Clone, Copy)]
pub struct CigarOp {
    pub op: CigarType,
    pub len: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarType {
    Match,      // M or = (0, 7)
    Insertion,  // I (1)
    Deletion,   // D (2)
    Mismatch,   // X (8)
    SoftClip,   // S (4)
    HardClip,   // H (5)
}

impl CigarType {
    pub fn consumes_query(&self) -> bool {
        matches!(self, Self::Match | Self::Insertion | Self::Mismatch | Self::SoftClip)
    }

    pub fn consumes_reference(&self) -> bool {
        matches!(self, Self::Match | Self::Deletion | Self::Mismatch)
    }
}

#[derive(Debug)]
pub struct CuttingResult {
    pub sequence: Vec<u8>,
    pub quality: Vec<u8>,
    pub removed_coords: Vec<u64>,
    pub new_query_start: usize,
    pub new_query_end: usize,
}
```

#### File: `cutting/read.rs`

```rust
use super::types::*;

pub fn position_in_or_before_primer(
    position: u64,
    primer_positions: &HashSet<u64>,
    lookaround: u64,
) -> bool {
    let search_range = position.saturating_sub(lookaround)..=position;
    search_range.into_iter().any(|p| primer_positions.contains(&p))
}

pub fn position_in_or_after_primer(
    position: u64,
    primer_positions: &HashSet<u64>,
    lookaround: u64,
) -> bool {
    let search_range = position..=position + lookaround;
    search_range.into_iter().any(|p| primer_positions.contains(&p))
}

pub fn cut_read(read: &mut Read, params: &mut CuttingParameters) -> CuttingResult {
    let mut removed_coords = Vec::new();
    
    // Determine starting position based on read and cut direction
    let mut position_on_sequence = if params.read_direction as i8 == match params.cut_direction {
        CutDirection::Forward => 1,
        CutDirection::Reverse => -1,
    } {
        params.query_start
    } else {
        params.query_end
    };

    let position_checker: fn(u64, &HashSet<u64>, u64) -> bool = match params.cut_direction {
        CutDirection::Forward => position_in_or_before_primer,
        CutDirection::Reverse => position_in_or_after_primer,
    };

    let cigar_iter: Box<dyn Iterator<Item = &CigarOp>> = match params.cut_direction {
        CutDirection::Forward => Box::new(params.cigar.iter()),
        CutDirection::Reverse => Box::new(params.cigar.iter().rev()),
    };

    for op in cigar_iter {
        let mut remaining = op.len as usize;
        
        while remaining > 0 {
            let needs_cutting = position_checker(
                params.position_on_reference,
                params.primer_positions,
                params.fragment_lookaround_size,
            );
            
            // Check if we should continue cutting
            if !needs_cutting && matches!(op.op, CigarType::Match) {
                break;
            }
            
            if needs_cutting || !matches!(op.op, CigarType::Match) {
                remaining -= 1;
                removed_coords.push(params.position_on_reference);
                
                // Update sequence position for consuming operations
                if op.op.consumes_query() {
                    let delta = match params.cut_direction {
                        CutDirection::Forward => params.read_direction,
                        CutDirection::Reverse => -params.read_direction,
                    };
                    position_on_sequence = (position_on_sequence as i64 + delta as i64) as usize;
                }
                
                // Update reference position
                if op.op.consumes_reference() {
                    let delta = match params.cut_direction {
                        CutDirection::Forward => 1i64,
                        CutDirection::Reverse => -1i64,
                    };
                    params.position_on_reference = 
                        (params.position_on_reference as i64 + delta) as u64;
                }
            }
        }
    }

    // Calculate new sequence bounds and extract cut sequence
    let (new_start, new_end) = calculate_new_bounds(
        params.query_start,
        params.query_end,
        position_on_sequence,
        params.cut_direction,
    );

    CuttingResult {
        sequence: read.sequence[new_start..new_end].to_vec(),
        quality: read.quality[new_start..new_end].to_vec(),
        removed_coords,
        new_query_start: new_start,
        new_query_end: new_end,
    }
}
```

#### File: `cutting/parallel.rs`

```rust
use rayon::prelude::*;
use crossbeam_channel::{bounded, Sender, Receiver};
use std::thread;

pub struct ParallelCutter {
    threads: usize,
    reference: PathBuf,
    preset: AlignmentPreset,
    scoring: ScoringMatrix,
    primer_index: PrimerIndex,
    amplicon_type: AmpliconType,
    fragment_lookaround: u64,
}

impl ParallelCutter {
    pub fn process_reads(
        &self,
        reads: Vec<FastqRecord>,
    ) -> Vec<ProcessedRead> {
        reads.into_par_iter()
            .filter_map(|record| self.process_single_read(record))
            .collect()
    }

    pub fn process_reads_streaming<I, W>(
        &self,
        input: I,
        mut output: W,
        chunk_size: usize,
    ) -> std::io::Result<ProcessingStats>
    where
        I: Iterator<Item = std::io::Result<FastqRecord>>,
        W: FastqWriter,
    {
        let mut stats = ProcessingStats::default();
        let mut chunk = Vec::with_capacity(chunk_size);
        
        for result in input {
            let record = result?;
            chunk.push(record);
            
            if chunk.len() >= chunk_size {
                let processed = self.process_reads(std::mem::take(&mut chunk));
                for read in processed {
                    stats.add(&read);
                    output.write_record(&read.to_fastq_record())?;
                }
                chunk = Vec::with_capacity(chunk_size);
            }
        }
        
        // Process remaining reads
        if !chunk.is_empty() {
            let processed = self.process_reads(chunk);
            for read in processed {
                stats.add(&read);
                output.write_record(&read.to_fastq_record())?;
            }
        }
        
        Ok(stats)
    }

    fn process_single_read(&self, record: FastqRecord) -> Option<ProcessedRead> {
        if record.sequence.len() < 42 {
            return None;  // Too short for alignment
        }

        let mut read = Read::from(record);
        let mut removed_coords = Vec::new();
        let max_iterations = 10;
        
        // Thread-local aligner (minimap2)
        thread_local! {
            static ALIGNER: RefCell<Option<Aligner>> = RefCell::new(None);
        }
        
        ALIGNER.with(|aligner| {
            let mut aligner = aligner.borrow_mut();
            if aligner.is_none() {
                *aligner = Some(Aligner::new(&self.reference, &self.preset, &self.scoring));
            }
            let aligner = aligner.as_mut().unwrap();
            
            for _ in 0..max_iterations {
                let prev_seq = read.sequence.clone();
                
                if let Some(hit) = aligner.map(&read.sequence).next() {
                    let (fw_coords, rv_coords) = self.cut_for_hit(&mut read, &hit);
                    removed_coords.extend(fw_coords);
                    removed_coords.extend(rv_coords);
                    
                    if read.sequence == prev_seq {
                        break;  // No more cutting needed
                    }
                } else {
                    break;  // No alignment found
                }
            }
        });

        Some(ProcessedRead {
            name: read.name,
            sequence: read.sequence,
            quality: read.quality,
            removed_coords,
        })
    }
}

#[derive(Debug, Default)]
pub struct ProcessingStats {
    pub total_reads: usize,
    pub processed_reads: usize,
    pub skipped_reads: usize,
    pub total_removed_bases: usize,
    pub unique_removed_coords: HashSet<u64>,
}
```

### Migration Checklist
- [ ] Define CuttingParameters and related types
- [ ] Implement position checking functions
- [ ] Port cut_read function with CIGAR handling
- [ ] Implement parallel processing with rayon
- [ ] Add streaming/chunked processing
- [ ] Implement thread-local aligner caching
- [ ] Add iteration limit and convergence detection
- [ ] Create comprehensive tests for edge cases

---

## Module 6: `log.py` → Rust logging

### Current Python Implementation
Uses custom Rich-based logging with colored output.

### Rust Implementation Strategy

Use the `tracing` crate with `tracing-subscriber`:

```rust
use tracing::{info, warn, error, debug, Level};
use tracing_subscriber::{fmt, EnvFilter};
use indicatif::{ProgressBar, ProgressStyle};

pub fn init_logging(verbose: bool, quiet: bool) {
    let level = if verbose {
        Level::DEBUG
    } else if quiet {
        Level::ERROR
    } else {
        Level::INFO
    };

    tracing_subscriber::fmt()
        .with_max_level(level)
        .with_target(false)
        .with_ansi(true)
        .init();
}

pub fn create_progress_bar(total: u64, message: &str) -> ProgressBar {
    let pb = ProgressBar::new(total);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("#>-")
    );
    pb.set_message(message.to_string());
    pb
}
```

### Migration Checklist
- [ ] Set up tracing with colored output
- [ ] Implement progress bars with indicatif
- [ ] Add log level configuration
- [ ] Port log messages to match Python output format