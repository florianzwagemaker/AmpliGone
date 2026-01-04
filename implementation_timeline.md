# AmpliGone Rust Implementation Timeline

## Overview

This document provides a detailed week-by-week implementation plan for rewriting AmpliGone in Rust, including milestones, deliverables, and acceptance criteria.

## Phase 1: Foundation (Weeks 1-3)

### Week 1: Project Setup and CLI

**Objectives:**
- Initialize Rust project with proper structure
- Implement complete CLI argument parsing
- Set up development infrastructure

**Tasks:**

1. **Project Initialization**
   ```bash
   cargo new ampligone
   cd ampligone
   # Set up directory structure as outlined in rust_rewrite_overview.md
   ```

2. **Cargo.toml Setup**
   - Add all dependencies from library_alternatives_updated.md
   - Configure release profiles
   - Set up benchmarks

3. **CLI Implementation** (`src/cli/args.rs`)
   - Port all arguments from Python's `args.py`
   - Add new paired-end arguments (`--input2`, `--output2`)
   - Add new processing mode arguments
   - Implement file validation
   - Add shell completion generation

4. **Logging Setup** (`src/utils/logging.rs`)
   - Configure tracing with colored output
   - Implement verbosity levels
   - Add progress bar utilities

**Deliverables:**
- [ ] Compiling project with all dependencies
- [ ] `ampligone --help` shows all arguments
- [ ] `ampligone --version` works
- [ ] Shell completions generate correctly

**Acceptance Criteria:**
```bash
# All these should work:
ampligone --help
ampligone --version
ampligone --input test.fq --output out.fq --reference ref.fa --primers primers.fa
# Should fail with proper error:
ampligone --input nonexistent.fq  # File not found error
```

---

### Week 2: I/O Operations

**Objectives:**
- Implement FASTQ reading/writing
- Implement BAM reading
- Implement BED file handling
- Implement FASTA parsing

**Tasks:**

1. **FASTQ I/O** (`src/io/fastq.rs`)
   ```rust
   // Core types
   pub struct FastqRecord { ... }
   pub struct FastqReader<R: BufRead> { ... }
   pub struct FastqWriter<W: Write> { ... }
   
   // Features:
   // - Gzip detection and handling
   // - Streaming iteration
   // - Buffered writing
   // - Parallel gzip output
   ```

2. **BAM Reading** (`src/io/bam.rs`)
   ```rust
   pub struct BamReader { ... }
   // - Use rust-htslib
   // - Handle unmapped reads
   // - Reverse complement for reverse strand
   ```

3. **BED File Handling** (`src/io/bed.rs`)
   ```rust
   pub struct BedRecord { ... }
   pub fn read_bed(path: &Path) -> Result<Vec<BedRecord>>
   pub fn write_bed(path: &Path, records: &[BedRecord]) -> Result<()>
   ```

4. **FASTA Parsing** (`src/io/fasta.rs`)
   ```rust
   pub struct FastaRecord { ... }
   pub fn read_fasta(path: &Path) -> Result<Vec<FastaRecord>>
   ```

**Deliverables:**
- [ ] FastqReader correctly parses test files
- [ ] FastqWriter produces valid FASTQ
- [ ] Gzip handling works for both read/write
- [ ] BAM reading extracts sequences correctly
- [ ] BED parsing matches Python output

**Tests:**
```rust
#[test]
fn test_fastq_roundtrip() {
    let records = vec![...];
    let temp = tempfile::NamedTempFile::new().unwrap();
    write_fastq(&temp.path(), &records).unwrap();
    let read_back = read_fastq(&temp.path()).collect::<Vec<_>>();
    assert_eq!(records, read_back);
}
```

---

### Week 3: DNA Utilities and Primer Structures

**Objectives:**
- Implement DNA sequence operations
- Implement IUPAC handling
- Create primer data structures
- Implement primer index

**Tasks:**

1. **DNA Operations** (`src/utils/dna.rs`)
   - Reverse complement
   - Complement only
   - DNA validation
   - Base quality conversion

2. **IUPAC Handling** (`src/utils/iupac.rs`)
   - Ambiguity expansion
   - Ambiguity detection
   - phf compile-time maps

3. **Primer Types** (`src/primer/types.rs`)
   ```rust
   pub struct PrimerCoordinates { ... }
   pub enum Strand { Forward, Reverse }
   ```

4. **Primer Index** (`src/primer/index.rs`)
   ```rust
   pub struct PrimerIndex {
       forward: HashMap<String, HashSet<u64>>,
       reverse: HashMap<String, HashSet<u64>>,
   }
   ```

**Deliverables:**
- [ ] reverse_complement matches Python behavior
- [ ] IUPAC expansion produces correct combinations
- [ ] PrimerIndex can be built from coordinates
- [ ] All unit tests pass

---

## Phase 2: Core Alignment (Weeks 4-6)

### Week 4: Primer Search with rust-bio

**Objectives:**
- Implement semi-global alignment using rust-bio
- Port primer coordinate finding
- Handle ambiguous primers

**Tasks:**

1. **Bio Aligner Wrapper** (`src/primer/alignment.rs`)
   ```rust
   pub struct PrimerAligner {
       scoring: Scoring,
       error_rate: f64,
   }
   
   impl PrimerAligner {
       pub fn find_primer_coords(
           &self,
           primer_seq: &[u8],
           reference: &[u8],
       ) -> Option<PrimerCoordinates>
   }
   ```

2. **Primer Search** (`src/primer/search.rs`)
   - Port `coord_list_gen` from Python
   - Handle IUPAC ambiguity
   - Choose best fitting coordinates

3. **CIGAR Parsing** (`src/utils/cigar.rs`)
   - Parse CIGAR strings
   - Count matches/mismatches/gaps
   - Calculate alignment statistics

**Deliverables:**
- [ ] Semi-global alignment produces same results as parasail
- [ ] Primer coordinates match Python output
- [ ] IUPAC primers correctly expanded and searched

**Validation:**
```python
# Run both Python and Rust on same inputs
python -m AmpliGone.fasta2bed --primers test.fa --reference ref.fa --output python.bed
ampligone-fasta2bed --primers test.fa --reference ref.fa --output rust.bed
diff python.bed rust.bed  # Should be identical
```

---

### Week 5: Minimap2 Integration

**Objectives:**
- Evaluate and integrate minimap2 bindings
- Implement thread-safe aligner
- Port alignment preset detection

**Tasks:**

1. **Evaluate minimap2-rs**
   - Test basic functionality
   - Check CIGAR access
   - Verify EQX flag support
   - Test performance

2. **If minimap2-rs insufficient, implement FFI**
   - See minimap2_integration_strategy.md
   - Create safe Rust wrapper
   - Handle memory management

3. **Alignment Preset Detection** (`src/alignment/preset.rs`)
   - Port from `alignmentpreset.py`
   - Calculate read statistics
   - Determine optimal preset

4. **Thread-Safe Aligner** (`src/alignment/pool.rs`)
   - Thread-local buffers
   - Shared index
   - Proper cleanup

**Deliverables:**
- [ ] Minimap2 alignments work correctly
- [ ] Preset detection matches Python
- [ ] Multi-threaded mapping works
- [ ] Memory usage is reasonable

---

### Week 6: Read Cutting Logic

**Objectives:**
- Port core cutting algorithm
- Implement position checking functions
- Handle all amplicon types

**Tasks:**

1. **Cutting Types** (`src/cutting/types.rs`)
   ```rust
   pub struct Read { ... }
   pub struct CuttingParameters { ... }
   pub struct CuttingResult { ... }
   pub enum AmpliconType { EndToEnd, EndToMid, Fragmented }
   ```

2. **Position Checking** (`src/cutting/position.rs`)
   ```rust
   pub fn position_in_or_before_primer(...) -> bool
   pub fn position_in_or_after_primer(...) -> bool
   ```

3. **Cut Read Function** (`src/cutting/read.rs`)
   - Port from Python's `cut_read`
   - Handle CIGAR operations
   - Track removed coordinates

4. **Amplicon Type Handling**
   - End-to-end: cut both ends
   - End-to-mid: cut based on strand
   - Fragmented: use lookaround

**Deliverables:**
- [ ] cut_read produces identical output to Python
- [ ] All amplicon types handled correctly
- [ ] Edge cases (short reads, no primers) handled

---

## Phase 3: Parallel Processing (Weeks 7-8)

### Week 7: Single-End Parallel Processing

**Objectives:**
- Implement parallel read processing
- Add streaming/chunked processing
- Optimize memory usage

**Tasks:**

1. **Parallel Cutter** (`src/cutting/parallel.rs`)
   ```rust
   pub struct ParallelCutter {
       threads: usize,
       aligner: Arc<Aligner>,
       primer_index: Arc<PrimerIndex>,
       // ...
   }
   
   impl ParallelCutter {
       pub fn process_reads(&self, reads: Vec<FastqRecord>) -> Vec<ProcessedRead>
       pub fn process_streaming<R, W>(&self, reader: R, writer: W) -> Stats
   }
   ```

2. **Chunked Processing**
   - Read in chunks of N reads
   - Process chunk in parallel
   - Write results immediately
   - Continue with next chunk

3. **Memory Optimization**
   - Reuse buffers where possible
   - Avoid unnecessary copies
   - Monitor memory usage

**Deliverables:**
- [ ] Parallel processing works correctly
- [ ] Speedup near-linear with thread count
- [ ] Memory usage bounded
- [ ] Streaming mode works

**Benchmarks:**
```rust
#[bench]
fn bench_parallel_processing(b: &mut Bencher) {
    // Test with varying thread counts
}
```

---

### Week 8: Paired-End Support

**Objectives:**
- Implement paired-end reading
- Implement paired-end processing
- Implement paired-end writing

**Tasks:**

1. **Paired Reader** (`src/paired/reader.rs`)
   - Synchronize R1 and R2 files
   - Validate read names match
   - Handle errors gracefully

2. **Paired Processor** (`src/paired/processor.rs`)
   - Process pairs together
   - Handle orphan reads
   - Implement processing modes (paired, relaxed)

3. **Paired Writer** (`src/paired/writer.rs`)
   - Synchronized output
   - Track statistics
   - Handle orphans appropriately

4. **Integration**
   - CLI arguments for paired mode
   - Main function branching
   - Statistics reporting

**Deliverables:**
- [ ] Paired-end reads processed correctly
- [ ] Read pairing validated
- [ ] Orphan handling works as expected
- [ ] Statistics accurate

---

## Phase 4: Enhanced Features (Weeks 9-10)

### Week 9: Enhanced Primer Detection

**Objectives:**
- Implement suspicious read detection
- Implement primer-to-read alignment
- Integrate enhanced processing mode

**Tasks:**

1. **Suspicious Read Detection** (`src/cutting/detection.rs`)
   - Compare alignment positions to primer sites
   - Flag reads starting/ending too far from primers
   - Calculate distance metrics

2. **Primer-to-Read Alignment** (`src/cutting/primer_alignment.rs`)
   - Use rust-bio semi-global alignment
   - Search for primers in read sequence
   - Handle reverse complement

3. **Enhanced Processor** (`src/cutting/enhanced.rs`)
   - Integrate detection and alignment
   - Additional cutting based on findings
   - Comprehensive logging

**Deliverables:**
- [ ] Suspicious reads correctly identified
- [ ] Primer alignment finds matches
- [ ] Enhanced mode properly cleans reads

---

### Week 10: Testing and Validation

**Objectives:**
- Comprehensive testing against Python
- Performance benchmarking
- Edge case handling

**Tasks:**

1. **Comparison Testing**
   - Run both versions on identical inputs
   - Compare output files byte-by-byte
   - Document any intentional differences

2. **Edge Case Testing**
   - Very short reads
   - Very long reads
   - No primers found
   - Empty files
   - Malformed input

3. **Performance Testing**
   - Benchmark various file sizes
   - Compare memory usage
   - Test different thread counts

4. **Fuzzing**
   - Use proptest for property-based testing
   - Fuzz with random inputs
   - Ensure no panics

**Deliverables:**
- [ ] All comparison tests pass
- [ ] Performance meets expectations
- [ ] No crashes on edge cases
- [ ] Fuzzing finds no issues

---

## Phase 5: Documentation and Release (Weeks 11-12)

### Week 11: Documentation

**Objectives:**
- User documentation
- API documentation
- Migration guide

**Tasks:**

1. **User Documentation**
   - Installation instructions
   - Usage examples
   - Changelog from Python version

2. **API Documentation**
   - Rustdoc for all public items
   - Examples in documentation
   - README updates

3. **Migration Guide**
   - Differences from Python version
   - New features (paired-end)
   - Known issues

---

### Week 12: Release Preparation

**Objectives:**
- Binary distribution
- Conda recipe
- CI/CD pipeline

**Tasks:**

1. **Binary Distribution**
   - Build for Linux (x86_64, aarch64)
   - Build for macOS (x86_64, aarch64)
   - Build for Windows (x86_64)

2. **Conda Recipe**
   - Update for Rust build
   - Test installation
   - Submit to bioconda

3. **CI/CD**
   - GitHub Actions workflow
   - Automated testing
   - Release automation

**Deliverables:**
- [ ] Binaries available for download
- [ ] Conda package works
- [ ] CI/CD fully functional
- [ ] Version 3.0.0 released

---

## Summary Timeline

| Week | Phase | Focus |
|------|-------|-------|
| 1 | Foundation | Project setup, CLI |
| 2 | Foundation | I/O operations |
| 3 | Foundation | DNA utilities, primer structures |
| 4 | Core Alignment | Primer search with rust-bio |
| 5 | Core Alignment | Minimap2 integration |
| 6 | Core Alignment | Read cutting logic |
| 7 | Parallel Processing | Single-end parallel |
| 8 | Parallel Processing | Paired-end support |
| 9 | Enhanced Features | Enhanced primer detection |
| 10 | Enhanced Features | Testing and validation |
| 11 | Release | Documentation |
| 12 | Release | Release preparation |

## Risk Mitigation

| Risk | Mitigation |
|------|------------|
| minimap2 FFI complexity | Start with minimap2-rs, have FFI fallback ready |
| Performance not meeting expectations | Continuous benchmarking, profile-guided optimization |
| Paired-end edge cases | Extensive testing with real-world data |
| rust-bio alignment differences | Validate against parasail results |
| Memory issues with large files | Implement streaming from the start |