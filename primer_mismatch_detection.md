# Enhanced Primer Mismatch Detection Strategy

## Overview

This document describes the enhanced primer detection and cleanup feature for reads where primers don't align properly to the reference. This addresses cases where reads start or stop too far from expected primer positions.

## Problem Statement

In amplicon sequencing, reads should ideally start/end at primer binding sites. However, several issues can cause misalignment:

1. **Primer mutations**: Primers may have SNPs relative to reference
2. **Incomplete extension**: Reads may start/end within primer regions
3. **Non-specific binding**: Primers may bind at unexpected locations
4. **Adapter contamination**: Residual adapter sequences
5. **Chimeric reads**: Fusion of multiple amplicons

## Detection Algorithm

### Phase 1: Identify Problematic Reads

```rust
// filepath: src/cutting/detection.rs

/// Threshold for considering a read position "too far" from expected primer site
const POSITION_TOLERANCE: u64 = 5;

/// Result of checking read alignment against expected primer positions
#[derive(Debug, Clone)]
pub enum ReadAlignmentStatus {
    /// Read aligns within expected primer boundaries
    Normal,
    /// Read starts too far from expected forward primer site
    SuspiciousStart {
        expected_start: u64,
        actual_start: u64,
        distance: u64,
    },
    /// Read ends too far from expected reverse primer site
    SuspiciousEnd {
        expected_end: u64,
        actual_end: u64,
        distance: u64,
    },
    /// Both start and end are suspicious
    SuspiciousBoth {
        start_info: (u64, u64, u64),  // expected, actual, distance
        end_info: (u64, u64, u64),
    },
}

/// Check if a read's alignment position matches expected primer sites
pub fn check_read_alignment(
    read_start: u64,
    read_end: u64,
    forward_primers: &HashSet<u64>,
    reverse_primers: &HashSet<u64>,
    tolerance: u64,
) -> ReadAlignmentStatus {
    // Find nearest forward primer to read start
    let nearest_forward = find_nearest_position(read_start, forward_primers);
    let start_distance = nearest_forward.map(|p| read_start.abs_diff(p));
    
    // Find nearest reverse primer to read end
    let nearest_reverse = find_nearest_position(read_end, reverse_primers);
    let end_distance = nearest_reverse.map(|p| read_end.abs_diff(p));
    
    let start_suspicious = start_distance.map_or(true, |d| d > tolerance);
    let end_suspicious = end_distance.map_or(true, |d| d > tolerance);
    
    match (start_suspicious, end_suspicious) {
        (false, false) => ReadAlignmentStatus::Normal,
        (true, false) => ReadAlignmentStatus::SuspiciousStart {
            expected_start: nearest_forward.unwrap_or(0),
            actual_start: read_start,
            distance: start_distance.unwrap_or(u64::MAX),
        },
        (false, true) => ReadAlignmentStatus::SuspiciousEnd {
            expected_end: nearest_reverse.unwrap_or(0),
            actual_end: read_end,
            distance: end_distance.unwrap_or(u64::MAX),
        },
        (true, true) => ReadAlignmentStatus::SuspiciousBoth {
            start_info: (
                nearest_forward.unwrap_or(0),
                read_start,
                start_distance.unwrap_or(u64::MAX),
            ),
            end_info: (
                nearest_reverse.unwrap_or(0),
                read_end,
                end_distance.unwrap_or(u64::MAX),
            ),
        },
    }
}

fn find_nearest_position(target: u64, positions: &HashSet<u64>) -> Option<u64> {
    positions
        .iter()
        .min_by_key(|&&p| target.abs_diff(p))
        .copied()
}
```

### Phase 2: Semi-Global Alignment for Suspicious Reads

When a read is flagged as suspicious, we perform semi-global alignment of primer sequences against the read:

```rust
// filepath: src/cutting/primer_alignment.rs

use bio::alignment::pairwise::{Aligner as BioAligner, Scoring};
use bio::alignment::AlignmentOperation;

/// Configuration for primer-to-read alignment
#[derive(Debug, Clone)]
pub struct PrimerAlignmentConfig {
    /// Match score
    pub match_score: i32,
    /// Mismatch penalty (negative)
    pub mismatch_penalty: i32,
    /// Gap open penalty (negative)
    pub gap_open: i32,
    /// Gap extend penalty (negative)
    pub gap_extend: i32,
    /// Minimum alignment score to consider a match (as fraction of perfect score)
    pub min_score_fraction: f64,
    /// Maximum number of mismatches allowed
    pub max_mismatches: usize,
}

impl Default for PrimerAlignmentConfig {
    fn default() -> Self {
        Self {
            match_score: 2,
            mismatch_penalty: -4,
            gap_open: -6,
            gap_extend: -2,
            min_score_fraction: 0.7,
            max_mismatches: 3,
        }
    }
}

/// Result of aligning a primer to a read
#[derive(Debug, Clone)]
pub struct PrimerAlignmentResult {
    /// Primer name/identifier
    pub primer_name: String,
    /// Start position in read (0-based)
    pub read_start: usize,
    /// End position in read (0-based, exclusive)
    pub read_end: usize,
    /// Alignment score
    pub score: i32,
    /// Number of matches
    pub matches: usize,
    /// Number of mismatches
    pub mismatches: usize,
    /// Number of gaps
    pub gaps: usize,
    /// Whether this is a valid primer match
    pub is_valid_match: bool,
}

/// Aligns primer sequences to a suspicious read using semi-global alignment
pub struct PrimerReadAligner {
    config: PrimerAlignmentConfig,
}

impl PrimerReadAligner {
    pub fn new(config: PrimerAlignmentConfig) -> Self {
        Self { config }
    }

    /// Align a single primer to a read using semi-global alignment
    /// 
    /// Semi-global alignment allows free gaps at the start/end of the pattern (primer)
    /// but not the text (read). This finds the best local match of the primer within the read.
    pub fn align_primer_to_read(
        &self,
        primer_name: &str,
        primer_seq: &[u8],
        read_seq: &[u8],
    ) -> Option<PrimerAlignmentResult> {
        // Use rust-bio's semi-global alignment
        let scoring = Scoring::new(
            self.config.gap_open,
            self.config.gap_extend,
            |a: u8, b: u8| {
                if a == b {
                    self.config.match_score
                } else {
                    self.config.mismatch_penalty
                }
            },
        );

        // Semi-global: free gaps at start/end of pattern (primer), not text (read)
        let mut aligner = BioAligner::with_scoring(scoring);
        
        // Perform semi-global alignment (pattern=primer, text=read)
        let alignment = aligner.semiglobal(primer_seq, read_seq);

        // Calculate statistics
        let mut matches = 0;
        let mut mismatches = 0;
        let mut gaps = 0;

        for op in &alignment.operations {
            match op {
                AlignmentOperation::Match => matches += 1,
                AlignmentOperation::Subst => mismatches += 1,
                AlignmentOperation::Ins | AlignmentOperation::Del => gaps += 1,
                AlignmentOperation::Xclip(_) | AlignmentOperation::Yclip(_) => {}
            }
        }

        // Calculate minimum acceptable score
        let perfect_score = primer_seq.len() as i32 * self.config.match_score;
        let min_score = (perfect_score as f64 * self.config.min_score_fraction) as i32;

        let is_valid_match = alignment.score >= min_score
            && mismatches <= self.config.max_mismatches;

        if alignment.score > 0 {
            Some(PrimerAlignmentResult {
                primer_name: primer_name.to_string(),
                read_start: alignment.ystart,
                read_end: alignment.yend,
                score: alignment.score,
                matches,
                mismatches,
                gaps,
                is_valid_match,
            })
        } else {
            None
        }
    }

    /// Search for any matching primer at the start of a read
    pub fn find_primer_at_start(
        &self,
        read_seq: &[u8],
        primers: &[(String, Vec<u8>)],
        search_window: usize,
    ) -> Option<PrimerAlignmentResult> {
        let search_region = &read_seq[..read_seq.len().min(search_window)];
        
        primers
            .iter()
            .filter_map(|(name, seq)| {
                self.align_primer_to_read(name, seq, search_region)
            })
            .filter(|r| r.is_valid_match && r.read_start < 10) // Must start near beginning
            .max_by_key(|r| r.score)
    }

    /// Search for any matching primer at the end of a read
    pub fn find_primer_at_end(
        &self,
        read_seq: &[u8],
        primers: &[(String, Vec<u8>)],
        search_window: usize,
    ) -> Option<PrimerAlignmentResult> {
        let start_pos = read_seq.len().saturating_sub(search_window);
        let search_region = &read_seq[start_pos..];
        
        primers
            .iter()
            .filter_map(|(name, seq)| {
                // For end primers, we often need the reverse complement
                let rc_seq = reverse_complement(seq);
                self.align_primer_to_read(name, &rc_seq, search_region)
            })
            .filter(|r| r.is_valid_match)
            .max_by_key(|r| r.score)
            .map(|mut r| {
                // Adjust positions to be relative to full read
                r.read_start += start_pos;
                r.read_end += start_pos;
                r
            })
    }
}

/// Calculate reverse complement of a DNA sequence
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&base| match base {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'G' | b'g' => b'C',
            b'C' | b'c' => b'G',
            b'N' | b'n' => b'N',
            other => other,
        })
        .collect()
}
```

### Phase 3: Enhanced Cutting Logic

```rust
// filepath: src/cutting/enhanced.rs

use super::detection::{check_read_alignment, ReadAlignmentStatus};
use super::primer_alignment::{PrimerReadAligner, PrimerAlignmentConfig, PrimerAlignmentResult};
use super::types::*;

/// Enhanced read processor that handles both normal and suspicious reads
pub struct EnhancedReadProcessor {
    /// Standard cutter for normal reads
    standard_cutter: ParallelCutter,
    /// Primer-to-read aligner for suspicious reads
    primer_aligner: PrimerReadAligner,
    /// Forward primer sequences (name, sequence)
    forward_primers: Vec<(String, Vec<u8>)>,
    /// Reverse primer sequences (name, sequence)  
    reverse_primers: Vec<(String, Vec<u8>)>,
    /// Position tolerance for suspicious read detection
    position_tolerance: u64,
    /// Search window size for primer alignment
    search_window: usize,
}

impl EnhancedReadProcessor {
    pub fn new(
        standard_cutter: ParallelCutter,
        forward_primers: Vec<(String, Vec<u8>)>,
        reverse_primers: Vec<(String, Vec<u8>)>,
    ) -> Self {
        Self {
            standard_cutter,
            primer_aligner: PrimerReadAligner::new(PrimerAlignmentConfig::default()),
            forward_primers,
            reverse_primers,
            position_tolerance: 5,
            search_window: 50, // Search first/last 50bp for primers
        }
    }

    /// Process a read with enhanced primer detection
    pub fn process_read(&self, read: &FastqRecord) -> Option<ProcessedRead> {
        // First, try standard processing
        let standard_result = self.standard_cutter.process_single_read(read.clone())?;
        
        // Get alignment information (this would come from the aligner)
        // For now, assume we have access to alignment coordinates
        let alignment_status = self.check_alignment_status(&standard_result);
        
        match alignment_status {
            ReadAlignmentStatus::Normal => Some(standard_result),
            
            ReadAlignmentStatus::SuspiciousStart { .. } => {
                self.handle_suspicious_start(read, standard_result)
            }
            
            ReadAlignmentStatus::SuspiciousEnd { .. } => {
                self.handle_suspicious_end(read, standard_result)
            }
            
            ReadAlignmentStatus::SuspiciousBoth { .. } => {
                self.handle_suspicious_both(read, standard_result)
            }
        }
    }

    fn handle_suspicious_start(
        &self,
        original_read: &FastqRecord,
        mut processed: ProcessedRead,
    ) -> Option<ProcessedRead> {
        // Try to find a primer at the start of the read
        if let Some(primer_match) = self.primer_aligner.find_primer_at_start(
            &original_read.sequence,
            &self.forward_primers,
            self.search_window,
        ) {
            tracing::debug!(
                "Found primer '{}' at start of read '{}' (positions {}-{})",
                primer_match.primer_name,
                original_read.name,
                primer_match.read_start,
                primer_match.read_end,
            );
            
            // Cut the primer from the processed read
            if primer_match.read_end < processed.sequence.len() {
                let additional_cut = primer_match.read_end;
                processed.sequence = processed.sequence[additional_cut..].to_vec();
                processed.quality = processed.quality[additional_cut..].to_vec();
                
                // Record the removed coordinates
                for i in 0..additional_cut {
                    processed.removed_coords.push(i as u64);
                }
            }
        }
        
        Some(processed)
    }

    fn handle_suspicious_end(
        &self,
        original_read: &FastqRecord,
        mut processed: ProcessedRead,
    ) -> Option<ProcessedRead> {
        // Try to find a primer at the end of the read
        if let Some(primer_match) = self.primer_aligner.find_primer_at_end(
            &original_read.sequence,
            &self.reverse_primers,
            self.search_window,
        ) {
            tracing::debug!(
                "Found primer '{}' at end of read '{}' (positions {}-{})",
                primer_match.primer_name,
                original_read.name,
                primer_match.read_start,
                primer_match.read_end,
            );
            
            // Calculate how much to cut from the processed read
            // This is complex because the processed read may already be shorter
            let original_len = original_read.sequence.len();
            let processed_len = processed.sequence.len();
            
            if primer_match.read_start < original_len {
                // Calculate relative position in processed read
                let cut_from_original_end = original_len - primer_match.read_start;
                let new_len = processed_len.saturating_sub(cut_from_original_end);
                
                if new_len > 0 {
                    processed.sequence.truncate(new_len);
                    processed.quality.truncate(new_len);
                    
                    // Record removed coordinates
                    for i in new_len..processed_len {
                        processed.removed_coords.push(i as u64);
                    }
                }
            }
        }
        
        Some(processed)
    }

    fn handle_suspicious_both(
        &self,
        original_read: &FastqRecord,
        processed: ProcessedRead,
    ) -> Option<ProcessedRead> {
        // Handle start first, then end
        let after_start = self.handle_suspicious_start(original_read, processed)?;
        self.handle_suspicious_end(original_read, after_start)
    }

    fn check_alignment_status(&self, _processed: &ProcessedRead) -> ReadAlignmentStatus {
        // This would be implemented with actual alignment coordinates
        // For now, return Normal as placeholder
        ReadAlignmentStatus::Normal
    }
}
```

## Integration with Main Pipeline

```rust
// filepath: src/cutting/pipeline.rs

/// Processing mode for reads
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum ProcessingMode {
    /// Standard processing (current behavior)
    Standard,
    /// Enhanced processing with primer mismatch detection
    Enhanced,
}

/// Main processing pipeline
pub struct ProcessingPipeline {
    mode: ProcessingMode,
    standard_cutter: ParallelCutter,
    enhanced_processor: Option<EnhancedReadProcessor>,
}

impl ProcessingPipeline {
    pub fn new(
        mode: ProcessingMode,
        cutter: ParallelCutter,
        primers: Option<(Vec<(String, Vec<u8>)>, Vec<(String, Vec<u8>)>)>,
    ) -> Self {
        let enhanced_processor = match mode {
            ProcessingMode::Standard => None,
            ProcessingMode::Enhanced => {
                let (fw, rv) = primers.expect("Enhanced mode requires primer sequences");
                Some(EnhancedReadProcessor::new(cutter.clone(), fw, rv))
            }
        };
        
        Self {
            mode,
            standard_cutter: cutter,
            enhanced_processor,
        }
    }

    pub fn process_reads(&self, reads: Vec<FastqRecord>) -> Vec<ProcessedRead> {
        match self.mode {
            ProcessingMode::Standard => {
                self.standard_cutter.process_reads(reads)
            }
            ProcessingMode::Enhanced => {
                let processor = self.enhanced_processor.as_ref().unwrap();
                reads
                    .into_par_iter()
                    .filter_map(|read| processor.process_read(&read))
                    .collect()
            }
        }
    }
}
```

## CLI Arguments for Enhanced Mode

```rust
// filepath: src/cli/args.rs (additions)

#[derive(Parser, Debug)]
pub struct Args {
    // ... existing args ...

    /// Processing mode
    #[arg(long, value_enum, default_value = "standard")]
    pub processing_mode: ProcessingMode,

    /// Position tolerance for suspicious read detection (bp)
    #[arg(long, default_value = "5")]
    pub position_tolerance: u64,

    /// Search window size for primer alignment (bp)
    #[arg(long, default_value = "50")]
    pub primer_search_window: usize,

    /// Minimum primer alignment score fraction (0.0-1.0)
    #[arg(long, default_value = "0.7")]
    pub min_primer_score: f64,

    /// Maximum primer mismatches allowed
    #[arg(long, default_value = "3")]
    pub max_primer_mismatches: usize,
}
```

## Testing Strategy

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_primer_alignment_perfect_match() {
        let aligner = PrimerReadAligner::new(PrimerAlignmentConfig::default());
        
        let primer = b"ACGTACGT";
        let read = b"NNNACGTACGTNNNNNNNNN";
        
        let result = aligner.align_primer_to_read("test_primer", primer, read);
        
        assert!(result.is_some());
        let result = result.unwrap();
        assert!(result.is_valid_match);
        assert_eq!(result.read_start, 3);
        assert_eq!(result.read_end, 11);
        assert_eq!(result.mismatches, 0);
    }

    #[test]
    fn test_primer_alignment_with_mismatches() {
        let aligner = PrimerReadAligner::new(PrimerAlignmentConfig::default());
        
        let primer = b"ACGTACGT";
        let read = b"NNNACGAACGTNNNNNNNNN"; // One mismatch (T->A)
        
        let result = aligner.align_primer_to_read("test_primer", primer, read);
        
        assert!(result.is_some());
        let result = result.unwrap();
        assert!(result.is_valid_match); // Should still be valid with 1 mismatch
        assert_eq!(result.mismatches, 1);
    }

    #[test]
    fn test_primer_alignment_too_many_mismatches() {
        let config = PrimerAlignmentConfig {
            max_mismatches: 1,
            ..Default::default()
        };
        let aligner = PrimerReadAligner::new(config);
        
        let primer = b"ACGTACGT";
        let read = b"NNNAAAAAAGANNNNNNNNN"; // Many mismatches
        
        let result = aligner.align_primer_to_read("test_primer", primer, read);
        
        // May find an alignment but it shouldn't be valid
        if let Some(result) = result {
            assert!(!result.is_valid_match);
        }
    }

    #[test]
    fn test_suspicious_read_detection() {
        let mut forward_primers = HashSet::new();
        forward_primers.insert(100u64);
        forward_primers.insert(200u64);
        
        let mut reverse_primers = HashSet::new();
        reverse_primers.insert(300u64);
        reverse_primers.insert(400u64);
        
        // Read starts at expected position
        let status = check_read_alignment(102, 298, &forward_primers, &reverse_primers, 5);
        assert!(matches!(status, ReadAlignmentStatus::Normal));
        
        // Read starts too far from expected
        let status = check_read_alignment(150, 298, &forward_primers, &reverse_primers, 5);
        assert!(matches!(status, ReadAlignmentStatus::SuspiciousStart { .. }));
    }
}
```

## Performance Considerations

1. **Lazy Evaluation**: Only perform enhanced alignment on suspicious reads
2. **Caching**: Cache primer reverse complements
3. **SIMD**: rust-bio uses SIMD acceleration for alignments
4. **Parallelization**: Use rayon for parallel processing of suspicious reads
5. **Early Exit**: Skip enhanced processing if standard processing succeeds cleanly

## Future Enhancements

1. **Machine Learning**: Train classifier for suspicious read detection
2. **Adapter Detection**: Extend to detect and remove adapter contamination
3. **Quality-Aware Alignment**: Weight alignment scores by base quality
4. **Report Generation**: Generate detailed reports on primer mismatches