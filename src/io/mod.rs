pub mod bam;
pub mod bed;
pub mod fasta;
pub mod fastq;

// Re-export commonly used types
pub use fastq::{
    FastqReader, FastqRecord, FastqWriter, PairedFastqIterator, PairedFastqReader,
    PairedFastqWriter,
};