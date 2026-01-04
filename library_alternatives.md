# Updated Library Alternatives and Dependencies

## Complete Dependency List (Updated)

### Cargo.toml

```toml
[package]
name = "ampligone"
version = "3.0.0"
edition = "2021"
authors = ["RIVM Bioinformatics"]
description = "Accurate primer removal from NGS reads in amplicon experiments"
license = "AGPL-3.0"
repository = "https://github.com/RIVM-bioinformatics/AmpliGone"
keywords = ["bioinformatics", "ngs", "primers", "amplicon"]
categories = ["science", "command-line-utilities"]
rust-version = "1.75"

[dependencies]
# CLI and argument parsing
clap = { version = "4.5", features = ["derive", "cargo", "wrap_help", "color"] }
clap_complete = "4.5"

# I/O operations
flate2 = "1.0"
gzp = { version = "0.11", features = ["deflate_rust"] }  # Parallel gzip
memmap2 = "0.9"

# Bioinformatics - UPDATED: Using rust-bio instead of parasailors
rust-htslib = { version = "0.47", features = ["bam", "bgzip"] }
bio = "2.0"  # Comprehensive bioinformatics library with alignment support
needletail = "0.6"  # Fast FASTQ/FASTA parsing
noodles = { version = "0.72", features = ["fastq", "fasta", "bam"] }  # Alternative I/O

# Alignment - minimap2 bindings
minimap2 = "0.1"  # Rust bindings (evaluate first)
# If custom FFI needed:
# bindgen = "0.69" (in build-dependencies)

# Parallelism
rayon = "1.10"
crossbeam = "0.8"
crossbeam-channel = "0.5"
num_cpus = "1.16"

# Data structures
hashbrown = { version = "0.14", features = ["rayon"] }
ahash = "0.8"
phf = { version = "0.11", features = ["macros"] }
indexmap = "2.2"

# Serialization (for potential JSON/config support)
serde = { version = "1.0", features = ["derive"] }
serde_json = "1.0"

# Error handling
thiserror = "1.0"
anyhow = "1.0"

# Logging and progress
tracing = "0.1"
tracing-subscriber = { version = "0.3", features = ["env-filter", "fmt"] }
indicatif = { version = "0.17", features = ["rayon"] }
console = "0.15"

# Utilities
once_cell = "1.19"
regex = "1.10"
itertools = "0.13"

[dev-dependencies]
criterion = { version = "0.5", features = ["html_reports"] }
proptest = "1.4"
tempfile = "3.10"
assert_cmd = "2.0"
predicates = "3.1"
insta = "1.36"  # Snapshot testing

[build-dependencies]
# Only if using custom minimap2 FFI:
# bindgen = "0.69"
# cc = "1.0"

[profile.release]
lto = "thin"
codegen-units = 1
panic = "abort"
strip = true
opt-level = 3

[profile.release-with-debug]
inherits = "release"
debug = true
strip = false

[profile.bench]
inherits = "release"
debug = true

[[bench]]
name = "benchmarks"
harness = false
```

---

## Key Library Changes from Original Plan

### Semi-Global Alignment: parasailors → rust-bio

**Reason for change**: `parasailors` is no longer maintained and may have compatibility issues with newer Rust versions.

**rust-bio advantages**:
- Actively maintained (part of rust-bio ecosystem)
- Pure Rust implementation (no FFI complexity)
- SIMD acceleration via `simdeez`
- Comprehensive alignment algorithms
- Well-tested and documented

#### Usage Comparison

**Old approach (parasailors)**:
```rust
// NOT RECOMMENDED - parasailors is unmaintained
use parasailors::{Matrix, Profile, semi_global_alignment};

fn align(primer: &[u8], reference: &[u8]) -> AlignmentResult {
    let matrix = Matrix::create("ACGT", 2, -1);
    let profile = Profile::new(primer, &matrix);
    semi_global_alignment(&profile, reference, 8, 30)
}
```

**New approach (rust-bio)**:
```rust
// filepath: src/primer/alignment.rs
use bio::alignment::pairwise::{Aligner, Scoring};
use bio::alignment::AlignmentOperation;

/// Semi-global alignment of primer to reference using rust-bio
pub fn align_primer_to_reference(
    primer: &[u8],
    reference: &[u8],
    match_score: i32,
    mismatch_penalty: i32,
    gap_open: i32,
    gap_extend: i32,
) -> PrimerAlignment {
    // Create scoring scheme
    let scoring = Scoring::new(gap_open, gap_extend, |a: u8, b: u8| {
        if a == b { match_score } else { mismatch_penalty }
    });
    
    let mut aligner = Aligner::with_scoring(scoring);
    
    // Semi-global alignment: gaps at pattern (primer) ends are free
    // This finds the best match of the primer within the reference
    let alignment = aligner.semiglobal(primer, reference);
    
    // Parse alignment results
    let mut matches = 0;
    let mut mismatches = 0;
    let mut insertions = 0;
    let mut deletions = 0;
    
    for op in &alignment.operations {
        match op {
            AlignmentOperation::Match => matches += 1,
            AlignmentOperation::Subst => mismatches += 1,
            AlignmentOperation::Ins => insertions += 1,
            AlignmentOperation::Del => deletions += 1,
            AlignmentOperation::Xclip(_) | AlignmentOperation::Yclip(_) => {}
        }
    }
    
    PrimerAlignment {
        start: alignment.ystart,
        end: alignment.yend,
        score: alignment.score,
        matches,
        mismatches,
        insertions,
        deletions,
        cigar: alignment_to_cigar(&alignment.operations),
    }
}

#[derive(Debug, Clone)]
pub struct PrimerAlignment {
    pub start: usize,
    pub end: usize,
    pub score: i32,
    pub matches: usize,
    pub mismatches: usize,
    pub insertions: usize,
    pub deletions: usize,
    pub cigar: String,
}

fn alignment_to_cigar(ops: &[AlignmentOperation]) -> String {
    use std::fmt::Write;
    
    let mut cigar = String::new();
    let mut current_op: Option<char> = None;
    let mut count = 0;
    
    for op in ops {
        let op_char = match op {
            AlignmentOperation::Match => '=',
            AlignmentOperation::Subst => 'X',
            AlignmentOperation::Ins => 'I',
            AlignmentOperation::Del => 'D',
            AlignmentOperation::Xclip(n) => {
                if *n > 0 {
                    write!(cigar, "{}S", n).unwrap();
                }
                continue;
            }
            AlignmentOperation::Yclip(n) => {
                if *n > 0 {
                    write!(cigar, "{}S", n).unwrap();
                }
                continue;
            }
        };
        
        if current_op == Some(op_char) {
            count += 1;
        } else {
            if let Some(prev_op) = current_op {
                write!(cigar, "{}{}", count, prev_op).unwrap();
            }
            current_op = Some(op_char);
            count = 1;
        }
    }
    
    if let Some(op) = current_op {
        write!(cigar, "{}{}", count, op).unwrap();
    }
    
    cigar
}
```

---

## Updated Module-Specific Libraries

### 1. FASTQ/FASTA Parsing

| Library | Use Case | Notes |
|---------|----------|-------|
| **needletail** | Primary FASTQ/FASTA parsing | Fastest option |
| **noodles** | Alternative with more features | Better error handling |
| **bio::io** | Fallback option | Part of rust-bio |

```rust
// Primary: needletail for speed
use needletail::{parse_fastx_file, Sequence};

pub fn read_fastq_needletail(path: &Path) -> impl Iterator<Item = FastqRecord> {
    parse_fastx_file(path)
        .expect("Failed to open file")
        .filter_map(|r| r.ok())
        .map(|record| FastqRecord {
            name: String::from_utf8_lossy(record.id()).into_owned(),
            sequence: record.seq().to_vec(),
            quality: record.qual().map(|q| q.to_vec()).unwrap_or_default(),
        })
}

// Alternative: noodles for paired-end with better validation
use noodles::fastq;

pub fn read_fastq_noodles(path: &Path) -> io::Result<impl Iterator<Item = io::Result<FastqRecord>>> {
    let file = File::open(path)?;
    let reader = if path.extension() == Some("gz".as_ref()) {
        fastq::Reader::new(BufReader::new(GzDecoder::new(file)))
    } else {
        fastq::Reader::new(BufReader::new(file))
    };
    
    Ok(reader.records().map(|r| {
        r.map(|record| FastqRecord {
            name: record.name().to_string(),
            sequence: record.sequence().to_vec(),
            quality: record.quality_scores().to_vec(),
        })
    }))
}
```

### 2. Sequence Alignment (rust-bio)

```rust
// filepath: src/alignment/bio_wrapper.rs
use bio::alignment::pairwise::{Aligner, Scoring};
use bio::alignment::sparse::HashMapBanded;

/// Wrapper for rust-bio alignment functionality
pub struct BioAligner {
    scoring: Scoring<fn(u8, u8) -> i32>,
    band_size: Option<usize>,
}

impl BioAligner {
    pub fn new(
        match_score: i32,
        mismatch_penalty: i32,
        gap_open: i32,
        gap_extend: i32,
    ) -> Self {
        let scoring = Scoring::new(
            gap_open,
            gap_extend,
            move |a: u8, b: u8| if a == b { match_score } else { mismatch_penalty },
        );
        
        Self {
            scoring,
            band_size: None,
        }
    }

    /// Enable banded alignment for better performance on long sequences
    pub fn with_band(mut self, band_size: usize) -> Self {
        self.band_size = Some(band_size);
        self
    }

    /// Perform semi-global alignment (pattern vs text)
    /// Free gaps at ends of pattern, not text
    pub fn semiglobal(&self, pattern: &[u8], text: &[u8]) -> bio::alignment::Alignment {
        let mut aligner = Aligner::with_scoring(self.scoring.clone());
        aligner.semiglobal(pattern, text)
    }

    /// Perform global alignment
    pub fn global(&self, seq1: &[u8], seq2: &[u8]) -> bio::alignment::Alignment {
        let mut aligner = Aligner::with_scoring(self.scoring.clone());
        aligner.global(seq1, seq2)
    }

    /// Perform local alignment
    pub fn local(&self, seq1: &[u8], seq2: &[u8]) -> bio::alignment::Alignment {
        let mut aligner = Aligner::with_scoring(self.scoring.clone());
        aligner.local(seq1, seq2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_semiglobal_alignment() {
        let aligner = BioAligner::new(2, -4, -6, -2);
        
        let pattern = b"ACGT";
        let text = b"NNNACGTNNN";
        
        let alignment = aligner.semiglobal(pattern, text);
        
        assert!(alignment.score > 0);
        assert_eq!(alignment.ystart, 3); // Start in text
        assert_eq!(alignment.yend, 7);   // End in text
    }
}
```

### 3. IUPAC Ambiguity Handling

```rust
// filepath: src/utils/iupac.rs
use phf::phf_map;

/// IUPAC ambiguity code mappings
static IUPAC_CODES: phf::Map<u8, &'static [u8]> = phf_map! {
    b'A' => b"A",
    b'C' => b"C",
    b'G' => b"G",
    b'T' => b"T",
    b'U' => b"T",  // RNA
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

/// Expand a sequence with IUPAC ambiguity codes into all possible unambiguous sequences
pub fn expand_ambiguous(seq: &[u8]) -> Vec<Vec<u8>> {
    let mut results = vec![Vec::with_capacity(seq.len())];
    
    for &base in seq {
        let base_upper = base.to_ascii_uppercase();
        let options = IUPAC_CODES
            .get(&base_upper)
            .copied()
            .unwrap_or(&[base_upper]);
        
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

/// Check if a sequence contains any IUPAC ambiguity codes
pub fn has_ambiguous_bases(seq: &[u8]) -> bool {
    seq.iter().any(|&b| {
        let upper = b.to_ascii_uppercase();
        !matches!(upper, b'A' | b'C' | b'G' | b'T')
    })
}

/// Count the number of possible sequences from an ambiguous sequence
pub fn ambiguity_count(seq: &[u8]) -> usize {
    seq.iter()
        .map(|&b| {
            IUPAC_CODES
                .get(&b.to_ascii_uppercase())
                .map(|opts| opts.len())
                .unwrap_or(1)
        })
        .product()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_expand_ambiguous() {
        let seq = b"ACR";  // R = A or G
        let expanded = expand_ambiguous(seq);
        
        assert_eq!(expanded.len(), 2);
        assert!(expanded.contains(&b"ACA".to_vec()));
        assert!(expanded.contains(&b"ACG".to_vec()));
    }

    #[test]
    fn test_no_ambiguity() {
        let seq = b"ACGT";
        let expanded = expand_ambiguous(seq);
        
        assert_eq!(expanded.len(), 1);
        assert_eq!(expanded[0], b"ACGT");
    }

    #[test]
    fn test_multiple_ambiguities() {
        let seq = b"RY";  // R=AG, Y=CT -> 4 combinations
        let expanded = expand_ambiguous(seq);
        
        assert_eq!(expanded.len(), 4);
    }
}
```

### 4. DNA Sequence Operations

```rust
// filepath: src/utils/dna.rs
use std::arch::x86_64::*;

/// Complement table for DNA bases
static COMPLEMENT: [u8; 256] = {
    let mut table = [0u8; 256];
    table[b'A' as usize] = b'T';
    table[b'T' as usize] = b'A';
    table[b'G' as usize] = b'C';
    table[b'C' as usize] = b'G';
    table[b'a' as usize] = b't';
    table[b't' as usize] = b'a';
    table[b'g' as usize] = b'c';
    table[b'c' as usize] = b'g';
    table[b'N' as usize] = b'N';
    table[b'n' as usize] = b'n';
    table
};

/// Calculate reverse complement of a DNA sequence
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| COMPLEMENT[b as usize])
        .collect()
}

/// Calculate reverse complement in place (modifies sequence)
pub fn reverse_complement_inplace(seq: &mut [u8]) {
    let len = seq.len();
    for i in 0..len / 2 {
        let j = len - 1 - i;
        let a = COMPLEMENT[seq[i] as usize];
        let b = COMPLEMENT[seq[j] as usize];
        seq[i] = b;
        seq[j] = a;
    }
    if len % 2 == 1 {
        let mid = len / 2;
        seq[mid] = COMPLEMENT[seq[mid] as usize];
    }
}

/// SIMD-accelerated reverse complement for longer sequences
#[cfg(target_arch = "x86_64")]
pub fn reverse_complement_simd(seq: &[u8]) -> Vec<u8> {
    if seq.len() < 32 {
        return reverse_complement(seq);
    }
    
    // For longer sequences, use SIMD
    // This is a placeholder - actual implementation would use AVX2/SSE4
    reverse_complement(seq)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT"); // Palindrome
        assert_eq!(reverse_complement(b"AACG"), b"CGTT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
    }

    #[test]
    fn test_reverse_complement_inplace() {
        let mut seq = b"AACG".to_vec();
        reverse_complement_inplace(&mut seq);
        assert_eq!(seq, b"CGTT");
    }
}
```

---

## Dependency Comparison Summary

| Feature | Python Library | Rust Library | Notes |
|---------|---------------|--------------|-------|
| CLI parsing | argparse + Rich | clap | More features, derive macros |
| FASTQ I/O | custom + pysam | needletail / noodles | 4-5x faster |
| BAM I/O | pysam | rust-htslib | Equivalent functionality |
| Semi-global alignment | parasail (via parasailors) | **rust-bio** | Pure Rust, maintained |
| Read alignment | mappy (minimap2) | minimap2-rs / FFI | Requires evaluation |
| Parallel processing | parmap | rayon | Much better performance |
| Progress bars | Rich | indicatif | Similar functionality |
| Gzip compression | pgzip | gzp | Parallel compression |

---

## Migration Checklist Update

### Phase 1: Core Infrastructure
- [ ] Set up Cargo project with updated dependencies
- [ ] Implement CLI with clap 4.5
- [ ] Implement FastqRecord and basic I/O
- [ ] Add logging with tracing

### Phase 2: Sequence Operations
- [ ] Implement IUPAC expansion with phf maps
- [ ] Implement reverse complement (with SIMD option)
- [ ] Port primer search using rust-bio semi-global alignment
- [ ] Create PrimerIndex structure

### Phase 3: Alignment Integration
- [ ] Evaluate minimap2-rs crate
- [ ] If needed, implement custom FFI bindings
- [ ] Create thread-safe aligner wrapper
- [ ] Implement CIGAR parsing

### Phase 4: Read Processing
- [ ] Port CuttingParameters and types
- [ ] Implement cut_read function
- [ ] Add parallel processing with rayon
- [ ] Implement streaming/chunked processing

### Phase 5: Paired-End Support
- [ ] Implement PairedFastqReader
- [ ] Implement PairedFastqWriter
- [ ] Add paired-end processing logic
- [ ] Add proper R1/R2 handling

### Phase 6: Enhanced Features
- [ ] Implement suspicious read detection
- [ ] Add primer-to-read alignment
- [ ] Integrate enhanced processing mode
- [ ] Add comprehensive testing

### Phase 7: Testing & Documentation
- [ ] Port all unit tests
- [ ] Add integration tests
- [ ] Benchmark against Python version
- [ ] Write user documentation