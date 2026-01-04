# AmpliGone Rust Rewrite: Overview and Architecture

## Executive Summary

This document outlines the comprehensive plan to rewrite AmpliGone from Python to Rust. The primary goals are:

1. **Performance**: Achieve significant speedup through Rust's zero-cost abstractions and memory safety
2. **Memory Efficiency**: Reduce memory footprint for large datasets
3. **Parallelism**: Leverage Rust's fearless concurrency for better multi-threading
4. **Maintainability**: Strong type system and compile-time guarantees

## Current Python Architecture

### Module Dependency Graph

```
__main__.py
    ├── args.py (CLI argument parsing)
    ├── io_ops.py (FASTQ/BAM I/O, SequenceReads class)
    ├── fasta2bed.py (Primer coordinate finding)
    ├── alignmentpreset.py (Automatic preset detection)
    ├── alignmentmatrix.py (Scoring matrix handling)
    ├── cut_reads.py (Core primer removal logic)
    ├── cutlery.py (Helper functions)
    └── log.py (Logging utilities)
```

### Core Data Flow

```
Input FASTQ/BAM → Index Reads → Find Primers → Determine Preset → 
    → Parallel Read Processing → Cut Primers → Write Output FASTQ
```

## Proposed Rust Architecture

### Crate Structure

```
ampligone/
├── Cargo.toml
├── src/
│   ├── main.rs              # Entry point, CLI handling
│   ├── lib.rs               # Library root, public API
│   ├── cli/
│   │   ├── mod.rs
│   │   └── args.rs          # Argument parsing with clap
│   ├── io/
│   │   ├── mod.rs
│   │   ├── fastq.rs         # FASTQ reading/writing
│   │   ├── bam.rs           # BAM reading
│   │   ├── bed.rs           # BED file handling
│   │   └── fasta.rs         # FASTA parsing
│   ├── primer/
│   │   ├── mod.rs
│   │   ├── search.rs        # Primer coordinate finding
│   │   ├── index.rs         # Primer indexing structures
│   │   └── alignment.rs     # Semi-global alignment for primers
│   ├── alignment/
│   │   ├── mod.rs
│   │   ├── minimap2.rs      # Minimap2 FFI wrapper
│   │   ├── preset.rs        # Automatic preset detection
│   │   └── scoring.rs       # Scoring matrix handling
│   ├── cutting/
│   │   ├── mod.rs
│   │   ├── read.rs          # Single read cutting logic
│   │   ├── parallel.rs      # Parallel processing orchestration
│   │   └── types.rs         # AmpliconType enum, CuttingParameters
│   ├── paired/
│   │   ├── mod.rs
│   │   ├── reader.rs        # Paired-end FASTQ reading
│   │   └── processor.rs     # R1/R2 coordinated processing
│   └── utils/
│       ├── mod.rs
│       ├── complement.rs    # DNA complement operations
│       └── cigar.rs         # CIGAR string parsing
├── benches/
│   └── benchmarks.rs        # Performance benchmarks
└── tests/
    ├── integration/
    └── data/
```

### Key Design Decisions

#### 1. Zero-Copy Where Possible
Use `&str` and `&[u8]` slices instead of owned strings where possible. Leverage Rust's lifetime system to avoid unnecessary allocations.

#### 2. SIMD-Accelerated Operations
Use SIMD for:
- Reverse complement calculations
- Quality score processing
- Sequence matching

#### 3. Memory-Mapped I/O
Use memory-mapped files for large FASTQ inputs to reduce memory pressure and improve I/O performance.

#### 4. Streaming Processing
Process reads in chunks rather than loading everything into memory:
- Read chunks of N reads
- Process in parallel
- Write output immediately
- Continue with next chunk

#### 5. Type-Safe Enums
Replace string-based type checks with proper Rust enums:

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AmpliconType {
    EndToEnd,
    EndToMid,
    Fragmented,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentPreset {
    ShortRead,      // sr
    MapOnt,         // map-ont
    MapPb,          // map-pb (if needed)
}
```

## Migration Strategy

### Phase 1: Foundation (Weeks 1-2)
- Set up Rust project structure
- Implement CLI argument parsing
- Implement basic I/O (FASTQ reading/writing)
- Establish testing infrastructure

### Phase 2: Core Logic (Weeks 3-5)
- Implement minimap2 FFI bindings
- Port primer search functionality
- Implement read cutting logic
- Add BED file support

### Phase 3: Parallelism & Optimization (Weeks 6-7)
- Implement parallel read processing with Rayon
- Add SIMD optimizations
- Implement streaming/chunked processing
- Memory optimization

### Phase 4: New Features (Weeks 8-9)
- Paired-end read support (R1/R2)
- Enhanced primer mismatch detection
- Semi-global alignment for primer cleaning

### Phase 5: Testing & Documentation (Weeks 10-11)
- Comprehensive testing against Python version
- Performance benchmarking
- Documentation
- CI/CD setup

### Phase 6: Release Preparation (Week 12)
- Binary distribution setup
- Conda recipe creation
- Final testing and bug fixes

## Performance Expectations

Based on similar bioinformatics tool rewrites (e.g., fastp, minimap2, seqtk):

| Metric | Python (Current) | Rust (Expected) |
|--------|-----------------|-----------------|
| Single-threaded speed | 1x (baseline) | 5-10x |
| Multi-threaded scaling | Limited by GIL | Near-linear |
| Memory usage | High (pandas DataFrames) | 50-70% reduction |
| Startup time | ~1-2s (Python import) | <100ms |

## Compatibility Considerations

### Input/Output Compatibility
- Maintain identical output format for drop-in replacement
- Support same input file formats (FASTQ, FASTQ.GZ, BAM)
- Identical BED output format for `--export-primers`

### CLI Compatibility
- Keep same argument names and semantics
- Add new arguments for paired-end support
- Maintain backwards compatibility with existing workflows

### Behavioral Compatibility
- Implement comprehensive test suite comparing Python vs Rust outputs
- Ensure identical primer detection results
- Match read cutting behavior exactly