# Paired-End Read Support Implementation

## Overview

Illumina sequencing produces paired-end reads (R1 and R2) that need to be processed together. This document describes the implementation strategy for proper paired-end support in AmpliGone.

## Current State

The Python version of AmpliGone does not have native paired-end support. Reads are processed independently, which can lead to:
- Orphaned reads (one of the pair is discarded)
- Inconsistent primer removal between pairs
- Loss of pairing information

## Requirements

### Input/Output
- Accept two input files: `--input` (R1) and `--input2` (R2)
- Produce two output files: `--output` (R1) and `--output2` (R2)
- Maintain read pairing throughout processing
- Handle unpaired/orphan reads appropriately

### Processing Rules
1. Both reads in a pair must be processed together
2. If one read is discarded (too short after cutting), both should be handled
3. Primer orientation awareness: R1 typically starts with forward primer, R2 with reverse
4. Support for different handling modes:
   - `paired`: Both reads must pass (default)
   - `relaxed`: Keep pairs where at least one read passes

## Data Structures

```rust
// filepath: src/paired/types.rs
use crate::io::fastq::FastqRecord;

/// A paired-end read consisting of R1 and R2
#[derive(Debug, Clone)]
pub struct ReadPair {
    /// Read 1 (forward)
    pub r1: FastqRecord,
    /// Read 2 (reverse)
    pub r2: FastqRecord,
}

impl ReadPair {
    /// Create a new read pair
    pub fn new(r1: FastqRecord, r2: FastqRecord) -> Self {
        Self { r1, r2 }
    }

    /// Check if both reads have the same base name (ignoring /1, /2 suffixes)
    pub fn is_properly_paired(&self) -> bool {
        let name1 = Self::base_name(&self.r1.name);
        let name2 = Self::base_name(&self.r2.name);
        name1 == name2
    }

    /// Extract base name without /1, /2 or _1, _2 suffix
    fn base_name(name: &str) -> &str {
        name.trim_end_matches("/1")
            .trim_end_matches("/2")
            .trim_end_matches(" 1")
            .trim_end_matches(" 2")
            .split_whitespace()
            .next()
            .unwrap_or(name)
    }

    /// Get the common base name for the pair
    pub fn base_name_owned(&self) -> String {
        Self::base_name(&self.r1.name).to_string()
    }
}

/// Result of processing a read pair
#[derive(Debug)]
pub enum PairProcessingResult {
    /// Both reads processed successfully
    BothPassed {
        r1: ProcessedRead,
        r2: ProcessedRead,
    },
    /// Only R1 passed (R2 too short or failed)
    OnlyR1Passed {
        r1: ProcessedRead,
        r2_reason: DiscardReason,
    },
    /// Only R2 passed (R1 too short or failed)
    OnlyR2Passed {
        r2: ProcessedRead,
        r1_reason: DiscardReason,
    },
    /// Both reads failed
    BothFailed {
        r1_reason: DiscardReason,
        r2_reason: DiscardReason,
    },
}

/// Reason why a read was discarded
#[derive(Debug, Clone)]
pub enum DiscardReason {
    /// Read too short after primer removal
    TooShort { remaining_length: usize },
    /// No alignment to reference found
    NoAlignment,
    /// Quality too low
    LowQuality,
    /// Read was entirely primer sequence
    EntirelyPrimer,
}

/// Mode for handling paired-end reads
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum PairedMode {
    /// Both reads must pass for pair to be kept
    Paired,
    /// Keep pairs where at least one read passes
    Relaxed,
    /// Process reads independently (legacy mode)
    Independent,
}
```

## Paired FASTQ Reader

```rust
// filepath: src/paired/reader.rs
use crate::io::fastq::{FastqReader, FastqRecord};
use std::io::BufRead;
use std::path::Path;

/// Iterator over paired-end reads from two FASTQ files
pub struct PairedFastqReader<R1: BufRead, R2: BufRead> {
    reader1: FastqReader<R1>,
    reader2: FastqReader<R2>,
    strict_pairing: bool,
}

impl<R1: BufRead, R2: BufRead> PairedFastqReader<R1, R2> {
    pub fn new(reader1: FastqReader<R1>, reader2: FastqReader<R2>) -> Self {
        Self {
            reader1,
            reader2,
            strict_pairing: true,
        }
    }

    /// Disable strict pairing checks (for performance)
    pub fn relaxed_pairing(mut self) -> Self {
        self.strict_pairing = false;
        self
    }
}

impl PairedFastqReader<Box<dyn BufRead>, Box<dyn BufRead>> {
    /// Open paired FASTQ files from paths
    pub fn from_paths<P1: AsRef<Path>, P2: AsRef<Path>>(
        path1: P1,
        path2: P2,
    ) -> std::io::Result<Self> {
        let reader1 = FastqReader::from_path(path1)?;
        let reader2 = FastqReader::from_path(path2)?;
        Ok(Self::new(reader1, reader2))
    }
}

impl<R1: BufRead, R2: BufRead> Iterator for PairedFastqReader<R1, R2> {
    type Item = Result<ReadPair, PairedReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        match (self.reader1.next(), self.reader2.next()) {
            (Some(Ok(r1)), Some(Ok(r2))) => {
                let pair = ReadPair::new(r1, r2);
                
                if self.strict_pairing && !pair.is_properly_paired() {
                    Some(Err(PairedReadError::MismatchedNames {
                        r1_name: pair.r1.name.clone(),
                        r2_name: pair.r2.name.clone(),
                    }))
                } else {
                    Some(Ok(pair))
                }
            }
            (None, None) => None,
            (Some(_), None) => Some(Err(PairedReadError::UnequalLength {
                file: "R2".to_string(),
            })),
            (None, Some(_)) => Some(Err(PairedReadError::UnequalLength {
                file: "R1".to_string(),
            })),
            (Some(Err(e)), _) | (_, Some(Err(e))) => {
                Some(Err(PairedReadError::IoError(e)))
            }
        }
    }
}

#[derive(Debug, thiserror::Error)]
pub enum PairedReadError {
    #[error("Read names do not match: R1='{r1_name}', R2='{r2_name}'")]
    MismatchedNames { r1_name: String, r2_name: String },
    
    #[error("Unequal file lengths: {file} ended early")]
    UnequalLength { file: String },
    
    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),
}
```

## Paired FASTQ Writer

```rust
// filepath: src/paired/writer.rs
use crate::io::fastq::{FastqWriter, FastqRecord};
use std::io::Write;
use std::path::Path;

/// Writer for paired-end FASTQ output
pub struct PairedFastqWriter<W1: Write, W2: Write> {
    writer1: FastqWriter<W1>,
    writer2: FastqWriter<W2>,
    written_pairs: usize,
    orphan_r1: usize,
    orphan_r2: usize,
}

impl<W1: Write, W2: Write> PairedFastqWriter<W1, W2> {
    pub fn new(writer1: FastqWriter<W1>, writer2: FastqWriter<W2>) -> Self {
        Self {
            writer1,
            writer2,
            written_pairs: 0,
            orphan_r1: 0,
            orphan_r2: 0,
        }
    }

    /// Write a complete pair
    pub fn write_pair(&mut self, r1: &FastqRecord, r2: &FastqRecord) -> std::io::Result<()> {
        self.writer1.write_record(r1)?;
        self.writer2.write_record(r2)?;
        self.written_pairs += 1;
        Ok(())
    }

    /// Write only R1 (orphan mode)
    pub fn write_r1_only(&mut self, r1: &FastqRecord) -> std::io::Result<()> {
        self.writer1.write_record(r1)?;
        self.orphan_r1 += 1;
        Ok(())
    }

    /// Write only R2 (orphan mode)  
    pub fn write_r2_only(&mut self, r2: &FastqRecord) -> std::io::Result<()> {
        self.writer2.write_record(r2)?;
        self.orphan_r2 += 1;
        Ok(())
    }

    /// Get writing statistics
    pub fn stats(&self) -> PairedWriterStats {
        PairedWriterStats {
            written_pairs: self.written_pairs,
            orphan_r1: self.orphan_r1,
            orphan_r2: self.orphan_r2,
        }
    }
}

impl PairedFastqWriter<Box<dyn Write>, Box<dyn Write>> {
    pub fn from_paths<P1: AsRef<Path>, P2: AsRef<Path>>(
        path1: P1,
        path2: P2,
        threads: usize,
    ) -> std::io::Result<Self> {
        let writer1 = FastqWriter::from_path(path1, threads)?;
        let writer2 = FastqWriter::from_path(path2, threads)?;
        Ok(Self::new(writer1, writer2))
    }
}

#[derive(Debug, Clone)]
pub struct PairedWriterStats {
    pub written_pairs: usize,
    pub orphan_r1: usize,
    pub orphan_r2: usize,
}
```

## Paired-End Processor

```rust
// filepath: src/paired/processor.rs
use super::types::*;
use crate::cutting::ParallelCutter;
use rayon::prelude::*;

pub struct PairedEndProcessor {
    cutter: ParallelCutter,
    mode: PairedMode,
    min_length: usize,
}

impl PairedEndProcessor {
    pub fn new(cutter: ParallelCutter, mode: PairedMode, min_length: usize) -> Self {
        Self {
            cutter,
            mode,
            min_length,
        }
    }

    /// Process a batch of read pairs in parallel
    pub fn process_pairs(&self, pairs: Vec<ReadPair>) -> Vec<PairProcessingResult> {
        pairs
            .into_par_iter()
            .map(|pair| self.process_single_pair(pair))
            .collect()
    }

    /// Process a single read pair
    fn process_single_pair(&self, pair: ReadPair) -> PairProcessingResult {
        // Process R1
        let r1_result = self.cutter.process_single_read(pair.r1.clone());
        
        // Process R2
        let r2_result = self.cutter.process_single_read(pair.r2.clone());

        // Evaluate results
        match (r1_result, r2_result) {
            (Some(r1), Some(r2)) => {
                // Check minimum length
                let r1_ok = r1.sequence.len() >= self.min_length;
                let r2_ok = r2.sequence.len() >= self.min_length;

                match (r1_ok, r2_ok) {
                    (true, true) => PairProcessingResult::BothPassed { r1, r2 },
                    (true, false) => PairProcessingResult::OnlyR1Passed {
                        r1,
                        r2_reason: DiscardReason::TooShort {
                            remaining_length: r2.sequence.len(),
                        },
                    },
                    (false, true) => PairProcessingResult::OnlyR2Passed {
                        r2,
                        r1_reason: DiscardReason::TooShort {
                            remaining_length: r1.sequence.len(),
                        },
                    },
                    (false, false) => PairProcessingResult::BothFailed {
                        r1_reason: DiscardReason::TooShort {
                            remaining_length: r1.sequence.len(),
                        },
                        r2_reason: DiscardReason::TooShort {
                            remaining_length: r2.sequence.len(),
                        },
                    },
                }
            }
            (Some(r1), None) => {
                if r1.sequence.len() >= self.min_length {
                    PairProcessingResult::OnlyR1Passed {
                        r1,
                        r2_reason: DiscardReason::NoAlignment,
                    }
                } else {
                    PairProcessingResult::BothFailed {
                        r1_reason: DiscardReason::TooShort {
                            remaining_length: r1.sequence.len(),
                        },
                        r2_reason: DiscardReason::NoAlignment,
                    }
                }
            }
            (None, Some(r2)) => {
                if r2.sequence.len() >= self.min_length {
                    PairProcessingResult::OnlyR2Passed {
                        r2,
                        r1_reason: DiscardReason::NoAlignment,
                    }
                } else {
                    PairProcessingResult::BothFailed {
                        r1_reason: DiscardReason::NoAlignment,
                        r2_reason: DiscardReason::TooShort {
                            remaining_length: r2.sequence.len(),
                        },
                    }
                }
            }
            (None, None) => PairProcessingResult::BothFailed {
                r1_reason: DiscardReason::NoAlignment,
                r2_reason: DiscardReason::NoAlignment,
            },
        }
    }

    /// Process pairs and write output based on mode
    pub fn process_and_write<W1: Write, W2: Write>(
        &self,
        pairs: impl Iterator<Item = Result<ReadPair, PairedReadError>>,
        writer: &mut PairedFastqWriter<W1, W2>,
        chunk_size: usize,
    ) -> Result<ProcessingStats, Box<dyn std::error::Error>> {
        let mut stats = ProcessingStats::default();
        let mut chunk = Vec::with_capacity(chunk_size);

        for pair_result in pairs {
            let pair = pair_result?;
            chunk.push(pair);

            if chunk.len() >= chunk_size {
                self.process_chunk(&mut chunk, writer, &mut stats)?;
            }
        }

        // Process remaining
        if !chunk.is_empty() {
            self.process_chunk(&mut chunk, writer, &mut stats)?;
        }

        Ok(stats)
    }

    fn process_chunk<W1: Write, W2: Write>(
        &self,
        chunk: &mut Vec<ReadPair>,
        writer: &mut PairedFastqWriter<W1, W2>,
        stats: &mut ProcessingStats,
    ) -> std::io::Result<()> {
        let results = self.process_pairs(std::mem::take(chunk));

        for result in results {
            match (&self.mode, result) {
                (_, PairProcessingResult::BothPassed { r1, r2 }) => {
                    writer.write_pair(&r1.to_fastq_record(), &r2.to_fastq_record())?;
                    stats.pairs_passed += 1;
                    stats.removed_bases_r1 += r1.removed_coords.len();
                    stats.removed_bases_r2 += r2.removed_coords.len();
                }
                (PairedMode::Relaxed, PairProcessingResult::OnlyR1Passed { r1, .. }) => {
                    writer.write_r1_only(&r1.to_fastq_record())?;
                    stats.orphan_r1 += 1;
                    stats.removed_bases_r1 += r1.removed_coords.len();
                }
                (PairedMode::Relaxed, PairProcessingResult::OnlyR2Passed { r2, .. }) => {
                    writer.write_r2_only(&r2.to_fastq_record())?;
                    stats.orphan_r2 += 1;
                    stats.removed_bases_r2 += r2.removed_coords.len();
                }
                (PairedMode::Paired, PairProcessingResult::OnlyR1Passed { .. })
                | (PairedMode::Paired, PairProcessingResult::OnlyR2Passed { .. }) => {
                    stats.pairs_discarded += 1;
                }
                (_, PairProcessingResult::BothFailed { .. }) => {
                    stats.pairs_discarded += 1;
                }
                (PairedMode::Independent, _) => {
                    // Handle independently - not used in paired mode
                }
            }
        }

        *chunk = Vec::with_capacity(chunk.capacity());
        Ok(())
    }
}

#[derive(Debug, Default)]
pub struct ProcessingStats {
    pub total_pairs: usize,
    pub pairs_passed: usize,
    pub pairs_discarded: usize,
    pub orphan_r1: usize,
    pub orphan_r2: usize,
    pub removed_bases_r1: usize,
    pub removed_bases_r2: usize,
}
```

## CLI Integration

```rust
// filepath: src/cli/args.rs (additions)
use super::*;

#[derive(Parser, Debug)]
pub struct Args {
    // ... existing args ...

    /// Second input file for paired-end reads (R2)
    #[arg(long = "input2", value_name = "FILE")]
    pub input2: Option<PathBuf>,

    /// Second output file for paired-end reads (R2)
    #[arg(long = "output2", value_name = "FILE")]
    pub output2: Option<PathBuf>,

    /// Paired-end processing mode
    #[arg(long, value_enum, default_value = "paired")]
    pub paired_mode: PairedMode,

    /// Minimum read length after primer removal
    #[arg(long, default_value = "20")]
    pub min_length: usize,
}

impl Args {
    pub fn is_paired_end(&self) -> bool {
        self.input2.is_some()
    }

    pub fn validate_paired_end(&self) -> Result<(), String> {
        match (&self.input2, &self.output2) {
            (Some(_), None) => Err("--output2 is required when using --input2".into()),
            (None, Some(_)) => Err("--input2 is required when using --output2".into()),
            _ => Ok(()),
        }
    }
}
```

## Main Entry Point

```rust
// filepath: src/main.rs (paired-end handling)
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    args.validate()?;
    
    init_logging(args.verbose, args.quiet);

    if args.is_paired_end() {
        run_paired_end(&args)?;
    } else {
        run_single_end(&args)?;
    }

    Ok(())
}

fn run_paired_end(args: &Args) -> Result<(), Box<dyn std::error::Error>> {
    let input2 = args.input2.as_ref().unwrap();
    let output2 = args.output2.as_ref().unwrap();

    info!("Running in paired-end mode");
    info!("R1 input: {}", args.input.display());
    info!("R2 input: {}", input2.display());

    // Create paired reader
    let reader = PairedFastqReader::from_paths(&args.input, input2)?;

    // Create paired writer
    let mut writer = PairedFastqWriter::from_paths(&args.output, output2, args.threads)?;

    // Create processor
    let cutter = create_cutter(args)?;
    let processor = PairedEndProcessor::new(cutter, args.paired_mode, args.min_length);

    // Process
    let stats = processor.process_and_write(reader, &mut writer, 10000)?;

    // Report statistics
    info!("Processing complete:");
    info!("  Pairs passed: {}", stats.pairs_passed);
    info!("  Pairs discarded: {}", stats.pairs_discarded);
    info!("  Orphan R1: {}", stats.orphan_r1);
    info!("  Orphan R2: {}", stats.orphan_r2);
    info!("  Bases removed from R1: {}", stats.removed_bases_r1);
    info!("  Bases removed from R2: {}", stats.removed_bases_r2);

    Ok(())
}
```

## Test Cases

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_pair_name_matching() {
        let r1 = FastqRecord {
            name: "read1/1".to_string(),
            sequence: b"ACGT".to_vec(),
            quality: b"IIII".to_vec(),
        };
        let r2 = FastqRecord {
            name: "read1/2".to_string(),
            sequence: b"ACGT".to_vec(),
            quality: b"IIII".to_vec(),
        };
        
        let pair = ReadPair::new(r1, r2);
        assert!(pair.is_properly_paired());
    }

    #[test]
    fn test_read_pair_name_mismatch() {
        let r1 = FastqRecord {
            name: "read1/1".to_string(),
            sequence: b"ACGT".to_vec(),
            quality: b"IIII".to_vec(),
        };
        let r2 = FastqRecord {
            name: "read2/2".to_string(),
            sequence: b"ACGT".to_vec(),
            quality: b"IIII".to_vec(),
        };
        
        let pair = ReadPair::new(r1, r2);
        assert!(!pair.is_properly_paired());
    }

    #[test]
    fn test_illumina_name_formats() {
        // Test various Illumina naming conventions
        let formats = [
            ("read1 1:N:0:ATCACG", "read1 2:N:0:ATCACG"),
            ("read1/1", "read1/2"),
            ("@M00001:1:000:1:1101:1:1 1:N:0", "@M00001:1:000:1:1101:1:1 2:N:0"),
        ];

        for (r1_name, r2_name) in formats {
            let pair = ReadPair::new(
                FastqRecord {
                    name: r1_name.to_string(),
                    sequence: vec![],
                    quality: vec![],
                },
                FastqRecord {
                    name: r2_name.to_string(),
                    sequence: vec![],
                    quality: vec![],
                },
            );
            assert!(pair.is_properly_paired(), "Failed for: {} vs {}", r1_name, r2_name);
        }
    }
}
```