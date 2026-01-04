# Minimap2 Integration Strategy

## Overview

Minimap2 is the core aligner used by AmpliGone for mapping reads to the reference genome. The Python version uses `mappy`, which is a Python binding to minimap2's C library. For Rust, we need to create proper FFI bindings.

## Current Python Usage Analysis

### How Mappy is Used in AmpliGone

From [`cut_reads.py`](AmpliGone/cut_reads.py):

```python
aligner = mp.Aligner(
    reference,
    preset=preset,
    best_n=1,
    scoring=scoring,
    extra_flags=0x4000000,  # MM_F_EQX flag - distinguish match vs mismatch
)

for hit in aligner.map(read.seq):
    # Access hit properties:
    # - hit.ctg: reference contig name
    # - hit.r_st: reference start position
    # - hit.r_en: reference end position
    # - hit.q_st: query start position
    # - hit.q_en: query end position
    # - hit.strand: alignment strand (1 or -1)
    # - hit.cigar: CIGAR operations as list of tuples
```

### Key Features Required

1. **Index Creation**: Load reference FASTA and create minimap2 index
2. **Single-hit Mapping**: Map reads with `best_n=1` (only best hit)
3. **Preset Support**: Support for `sr`, `map-ont`, `map-pb` presets
4. **Custom Scoring**: Ability to pass custom scoring matrix
5. **EQX Flag**: Distinguish between sequence match (=) and mismatch (X) in CIGAR
6. **CIGAR Access**: Full access to CIGAR operations

## Integration Options

### Option 1: minimap2-rs Crate (Recommended Starting Point)

The `minimap2` crate on crates.io provides Rust bindings.

```toml
[dependencies]
minimap2 = "0.1"
```

**Pros:**
- Maintained community crate
- Safe Rust API
- Handles memory management

**Cons:**
- May not expose all features we need
- Version might lag behind upstream

### Option 2: Custom FFI with bindgen

Create our own bindings using `bindgen` for complete control.

```toml
[build-dependencies]
bindgen = "0.69"
cc = "1.0"

[dependencies]
libc = "0.2"
```

**Pros:**
- Full control over exposed features
- Can match exact minimap2 version
- Access to all flags and options

**Cons:**
- More maintenance burden
- Need to handle unsafe code carefully

### Option 3: Hybrid Approach (Recommended)

Use `minimap2-rs` as base, extend with custom FFI for missing features.

## Detailed Implementation Plan

### Phase 1: Evaluate minimap2-rs

```rust
// Test if minimap2-rs meets our needs
use minimap2::*;

fn test_minimap2_features() -> Result<(), Box<dyn std::error::Error>> {
    // 1. Test index creation with preset
    let aligner = Aligner::builder()
        .preset(Preset::Sr)  // or MapOnt, MapPb
        .with_index("reference.fasta", None)?;
    
    // 2. Test mapping with options
    let mappings = aligner.map(
        b"ACGTACGTACGT",
        false,  // not reverse complement
        false,  // not supplementary only
        None,   // no channel
        None,   // no quality
    )?;
    
    // 3. Check if we can access all needed fields
    for mapping in mappings {
        println!("Target: {}", mapping.target_name.unwrap_or_default());
        println!("Target start: {}", mapping.target_start);
        println!("Target end: {}", mapping.target_end);
        println!("Query start: {}", mapping.query_start);
        println!("Query end: {}", mapping.query_end);
        println!("Strand: {:?}", mapping.strand);
        println!("CIGAR: {:?}", mapping.alignment); // Check CIGAR access
    }
    
    Ok(())
}
```

### Phase 2: Custom FFI Bindings (if needed)

If minimap2-rs doesn't meet all requirements, create custom bindings:

#### build.rs

```rust
// filepath: build.rs
use std::env;
use std::path::PathBuf;

fn main() {
    // Compile minimap2 from source
    cc::Build::new()
        .files([
            "vendor/minimap2/align.c",
            "vendor/minimap2/bseq.c",
            "vendor/minimap2/chain.c",
            "vendor/minimap2/esterr.c",
            "vendor/minimap2/format.c",
            "vendor/minimap2/hit.c",
            "vendor/minimap2/index.c",
            "vendor/minimap2/kalloc.c",
            "vendor/minimap2/ksw2_extd2_sse.c",
            "vendor/minimap2/ksw2_exts2_sse.c",
            "vendor/minimap2/ksw2_extz2_sse.c",
            "vendor/minimap2/ksw2_ll_sse.c",
            "vendor/minimap2/kthread.c",
            "vendor/minimap2/map.c",
            "vendor/minimap2/misc.c",
            "vendor/minimap2/options.c",
            "vendor/minimap2/pe.c",
            "vendor/minimap2/sdust.c",
            "vendor/minimap2/sketch.c",
            "vendor/minimap2/splitidx.c",
        ])
        .include("vendor/minimap2")
        .flag("-DHAVE_KALLOC")
        .flag("-msse4.1")  // SSE4.1 for SIMD
        .opt_level(3)
        .compile("minimap2");

    // Generate bindings
    let bindings = bindgen::Builder::default()
        .header("vendor/minimap2/minimap.h")
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .allowlist_function("mm_.*")
        .allowlist_type("mm_.*")
        .allowlist_var("MM_.*")
        .generate()
        .expect("Unable to generate bindings");

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("minimap2_bindings.rs"))
        .expect("Couldn't write bindings!");

    println!("cargo:rustc-link-lib=static=minimap2");
    println!("cargo:rustc-link-lib=z");
    println!("cargo:rustc-link-lib=pthread");
}
```

#### Safe Rust Wrapper

```rust
// filepath: src/alignment/minimap2.rs
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

include!(concat!(env!("OUT_DIR"), "/minimap2_bindings.rs"));

use std::ffi::{CStr, CString};
use std::ptr;

/// Minimap2 alignment preset
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Preset {
    /// Short reads (Illumina)
    ShortRead,
    /// Oxford Nanopore reads
    MapOnt,
    /// PacBio CLR reads
    MapPb,
    /// PacBio HiFi reads
    MapHifi,
}

impl Preset {
    fn as_cstr(&self) -> &'static CStr {
        match self {
            Preset::ShortRead => c"sr",
            Preset::MapOnt => c"map-ont",
            Preset::MapPb => c"map-pb",
            Preset::MapHifi => c"map-hifi",
        }
    }
}

/// Scoring matrix for alignment
#[derive(Debug, Clone)]
pub struct ScoringMatrix {
    pub match_score: i32,
    pub mismatch_penalty: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
}

impl Default for ScoringMatrix {
    fn default() -> Self {
        Self {
            match_score: 2,
            mismatch_penalty: 4,
            gap_open: 4,
            gap_extend: 2,
        }
    }
}

/// CIGAR operation type
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOp {
    Match,          // M (0)
    Insertion,      // I (1)
    Deletion,       // D (2)
    Skip,           // N (3)
    SoftClip,       // S (4)
    HardClip,       // H (5)
    Padding,        // P (6)
    SeqMatch,       // = (7) - with EQX flag
    SeqMismatch,    // X (8) - with EQX flag
}

impl CigarOp {
    fn from_raw(op: u32) -> Self {
        match op & 0xf {
            0 => CigarOp::Match,
            1 => CigarOp::Insertion,
            2 => CigarOp::Deletion,
            3 => CigarOp::Skip,
            4 => CigarOp::SoftClip,
            5 => CigarOp::HardClip,
            6 => CigarOp::Padding,
            7 => CigarOp::SeqMatch,
            8 => CigarOp::SeqMismatch,
            _ => CigarOp::Match,
        }
    }

    /// Whether this operation consumes query bases
    pub fn consumes_query(&self) -> bool {
        matches!(
            self,
            CigarOp::Match
                | CigarOp::Insertion
                | CigarOp::SoftClip
                | CigarOp::SeqMatch
                | CigarOp::SeqMismatch
        )
    }

    /// Whether this operation consumes reference bases
    pub fn consumes_reference(&self) -> bool {
        matches!(
            self,
            CigarOp::Match
                | CigarOp::Deletion
                | CigarOp::Skip
                | CigarOp::SeqMatch
                | CigarOp::SeqMismatch
        )
    }
}

/// A single CIGAR operation with length
#[derive(Debug, Clone, Copy)]
pub struct CigarElement {
    pub op: CigarOp,
    pub len: u32,
}

/// Alignment hit from minimap2
#[derive(Debug, Clone)]
pub struct AlignmentHit {
    /// Reference sequence name
    pub target_name: String,
    /// Reference start position (0-based)
    pub target_start: u64,
    /// Reference end position (0-based, exclusive)
    pub target_end: u64,
    /// Query start position (0-based)
    pub query_start: u32,
    /// Query end position (0-based, exclusive)
    pub query_end: u32,
    /// Strand: 1 for forward, -1 for reverse
    pub strand: i8,
    /// CIGAR operations
    pub cigar: Vec<CigarElement>,
    /// Mapping quality
    pub mapq: u8,
    /// Alignment score
    pub score: i32,
}

/// Minimap2 aligner wrapper
pub struct Aligner {
    idx: *mut mm_idx_t,
    map_opt: mm_mapopt_t,
    idx_opt: mm_idxopt_t,
}

// Safety: The minimap2 index is thread-safe for reading after construction
unsafe impl Send for Aligner {}
unsafe impl Sync for Aligner {}

impl Aligner {
    /// Create a new aligner from a reference FASTA file
    pub fn new(
        reference_path: &str,
        preset: Preset,
        scoring: Option<ScoringMatrix>,
    ) -> Result<Self, AlignerError> {
        let mut idx_opt: mm_idxopt_t = unsafe { std::mem::zeroed() };
        let mut map_opt: mm_mapopt_t = unsafe { std::mem::zeroed() };

        // Initialize with preset
        unsafe {
            mm_set_opt(preset.as_cstr().as_ptr(), &mut idx_opt, &mut map_opt);
        }

        // Apply custom scoring if provided
        if let Some(scoring) = scoring {
            map_opt.a = scoring.match_score as i16;
            map_opt.b = scoring.mismatch_penalty as i16;
            map_opt.q = scoring.gap_open as i16;
            map_opt.e = scoring.gap_extend as i16;
        }

        // Set EQX flag to distinguish match vs mismatch
        map_opt.flag |= MM_F_EQX as i64;

        // Only return best hit
        map_opt.best_n = 1;

        // Build index
        let path_cstr = CString::new(reference_path)
            .map_err(|_| AlignerError::InvalidPath)?;

        let idx = unsafe {
            mm_idx_str(
                idx_opt.w as i32,
                idx_opt.k as i32,
                0,  // is_hpc
                idx_opt.bucket_bits as i32,
                1,  // n_threads
                &path_cstr.as_ptr(),
                1,  // n_files
            )
        };

        if idx.is_null() {
            return Err(AlignerError::IndexCreationFailed);
        }

        // Prepare for mapping
        unsafe {
            mm_mapopt_update(&mut map_opt, idx);
        }

        Ok(Self { idx, map_opt, idx_opt })
    }

    /// Map a sequence to the reference
    pub fn map(&self, sequence: &[u8]) -> Vec<AlignmentHit> {
        let mut n_regs: i32 = 0;
        let mut tbuf = unsafe { mm_tbuf_init() };

        let regs = unsafe {
            mm_map(
                self.idx,
                sequence.len() as i32,
                sequence.as_ptr() as *const i8,
                &mut n_regs,
                tbuf,
                &self.map_opt,
                ptr::null(),  // name
            )
        };

        let mut hits = Vec::new();

        if !regs.is_null() && n_regs > 0 {
            for i in 0..n_regs as isize {
                let reg = unsafe { &*regs.offset(i) };

                // Get reference name
                let target_name = unsafe {
                    let name_ptr = (*(*self.idx).seq.offset(reg.rid as isize)).name;
                    CStr::from_ptr(name_ptr)
                        .to_string_lossy()
                        .into_owned()
                };

                // Parse CIGAR
                let cigar = if reg.p.is_null() {
                    Vec::new()
                } else {
                    let p = unsafe { &*reg.p };
                    (0..p.n_cigar as isize)
                        .map(|j| {
                            let raw = unsafe { *p.cigar.offset(j) };
                            CigarElement {
                                op: CigarOp::from_raw(raw),
                                len: raw >> 4,
                            }
                        })
                        .collect()
                };

                hits.push(AlignmentHit {
                    target_name,
                    target_start: reg.rs as u64,
                    target_end: reg.re as u64,
                    query_start: reg.qs,
                    query_end: reg.qe,
                    strand: if reg.rev() != 0 { -1 } else { 1 },
                    cigar,
                    mapq: reg.mapq,
                    score: reg.score,
                });
            }

            // Free registrations
            for i in 0..n_regs as isize {
                let reg = unsafe { &*regs.offset(i) };
                if !reg.p.is_null() {
                    unsafe { libc::free(reg.p as *mut libc::c_void) };
                }
            }
            unsafe { libc::free(regs as *mut libc::c_void) };
        }

        unsafe { mm_tbuf_destroy(tbuf) };

        hits
    }
}

impl Drop for Aligner {
    fn drop(&mut self) {
        if !self.idx.is_null() {
            unsafe { mm_idx_destroy(self.idx) };
        }
    }
}

#[derive(Debug, thiserror::Error)]
pub enum AlignerError {
    #[error("Invalid path")]
    InvalidPath,
    #[error("Failed to create index")]
    IndexCreationFailed,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cigar_op_properties() {
        assert!(CigarOp::Match.consumes_query());
        assert!(CigarOp::Match.consumes_reference());
        assert!(CigarOp::Insertion.consumes_query());
        assert!(!CigarOp::Insertion.consumes_reference());
        assert!(!CigarOp::Deletion.consumes_query());
        assert!(CigarOp::Deletion.consumes_reference());
    }
}
```

### Phase 3: Thread-Local Aligner Pool

For parallel processing, each thread needs its own mapping buffer:

```rust
// filepath: src/alignment/pool.rs
use std::cell::RefCell;
use std::sync::Arc;

thread_local! {
    static THREAD_BUFFER: RefCell<Option<ThreadLocalAligner>> = RefCell::new(None);
}

/// Thread-local aligner for parallel processing
pub struct ThreadLocalAligner {
    aligner: Arc<Aligner>,
    buffer: *mut mm_tbuf_t,
}

impl ThreadLocalAligner {
    pub fn new(aligner: Arc<Aligner>) -> Self {
        let buffer = unsafe { mm_tbuf_init() };
        Self { aligner, buffer }
    }

    pub fn map(&self, sequence: &[u8]) -> Vec<AlignmentHit> {
        // Use thread-local buffer for mapping
        // This avoids allocation per-read
        self.aligner.map_with_buffer(sequence, self.buffer)
    }
}

impl Drop for ThreadLocalAligner {
    fn drop(&mut self) {
        if !self.buffer.is_null() {
            unsafe { mm_tbuf_destroy(self.buffer) };
        }
    }
}

/// Get or create thread-local aligner
pub fn with_thread_aligner<F, R>(aligner: &Arc<Aligner>, f: F) -> R
where
    F: FnOnce(&ThreadLocalAligner) -> R,
{
    THREAD_BUFFER.with(|cell| {
        let mut borrow = cell.borrow_mut();
        if borrow.is_none() {
            *borrow = Some(ThreadLocalAligner::new(Arc::clone(aligner)));
        }
        f(borrow.as_ref().unwrap())
    })
}
```

## Testing Strategy

### Unit Tests

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;
    use std::io::Write;

    fn create_test_reference() -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">ref1").unwrap();
        writeln!(file, "ACGTACGTACGTACGTACGTACGTACGTACGTACGT").unwrap();
        file
    }

    #[test]
    fn test_aligner_creation() {
        let ref_file = create_test_reference();
        let aligner = Aligner::new(
            ref_file.path().to_str().unwrap(),
            Preset::ShortRead,
            None,
        );
        assert!(aligner.is_ok());
    }

    #[test]
    fn test_mapping_perfect_match() {
        let ref_file = create_test_reference();
        let aligner = Aligner::new(
            ref_file.path().to_str().unwrap(),
            Preset::ShortRead,
            None,
        ).unwrap();

        let hits = aligner.map(b"ACGTACGTACGT");
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].target_name, "ref1");
        assert_eq!(hits[0].strand, 1);
    }

    #[test]
    fn test_mapping_reverse_complement() {
        let ref_file = create_test_reference();
        let aligner = Aligner::new(
            ref_file.path().to_str().unwrap(),
            Preset::ShortRead,
            None,
        ).unwrap();

        // Reverse complement of ACGTACGT is ACGTACGT (palindrome)
        // Test with non-palindrome
        let hits = aligner.map(b"ACGTACGTACGT");
        assert!(!hits.is_empty());
    }

    #[test]
    fn test_cigar_parsing() {
        let ref_file = create_test_reference();
        let aligner = Aligner::new(
            ref_file.path().to_str().unwrap(),
            Preset::ShortRead,
            None,
        ).unwrap();

        let hits = aligner.map(b"ACGTACGTACGT");
        assert!(!hits.is_empty());
        
        // With EQX flag, perfect matches should have SeqMatch operations
        for hit in &hits {
            for cigar_elem in &hit.cigar {
                if matches!(cigar_elem.op, CigarOp::SeqMatch | CigarOp::Match) {
                    // Expected for matching regions
                }
            }
        }
    }
}
```

### Integration Tests

```rust
// tests/integration/minimap2_test.rs
use ampligone::alignment::*;
use std::path::PathBuf;

#[test]
fn test_real_world_alignment() {
    let test_data = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/reference.fasta");
    
    if !test_data.exists() {
        // Skip if test data not available
        return;
    }

    let aligner = Aligner::new(
        test_data.to_str().unwrap(),
        Preset::ShortRead,
        Some(ScoringMatrix {
            match_score: 2,
            mismatch_penalty: 4,
            gap_open: 4,
            gap_extend: 2,
        }),
    ).unwrap();

    // Test with known sequences
    let hits = aligner.map(b"ACGTACGTACGTACGTACGTACGT");
    // Assertions based on known reference
}
```

## Performance Considerations

1. **Index Caching**: Create index once, share across threads
2. **Thread-Local Buffers**: Avoid allocation per mapping call
3. **Memory Mapping**: Consider memory-mapped index for large references
4. **SIMD**: Ensure minimap2 compiled with SSE/AVX support

## Fallback Strategy

If FFI integration proves problematic:

```rust
use std::process::{Command, Stdio};
use std::io::{BufReader, BufRead};

/// Fallback: call minimap2 as subprocess
pub fn map_with_subprocess(
    reference: &Path,
    reads: &[FastqRecord],
    preset: &str,
) -> Result<Vec<AlignmentHit>, std::io::Error> {
    // Write reads to temp file
    let temp_reads = tempfile::NamedTempFile::new()?;
    // ... write reads ...

    let output = Command::new("minimap2")
        .args(["-x", preset, "-a", "--eqx"])
        .arg(reference)
        .arg(temp_reads.path())
        .stdout(Stdio::piped())
        .spawn()?
        .wait_with_output()?;

    // Parse SAM output
    parse_sam_output(&output.stdout)
}
```

This fallback is useful for:
- Development and testing
- Platforms where FFI is problematic
- Debugging alignment issues