use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Represents a single FASTQ record containing read information
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastqRecord {
    /// Read identifier (without the '@' prefix)
    pub name: String,
    /// Nucleotide sequence
    pub sequence: String,
    /// Optional comment line (usually empty)
    pub comment: String,
    /// Quality scores as ASCII characters
    pub qualities: String,
}

impl FastqRecord {
    /// Create a new FastqRecord
    pub fn new(name: String, sequence: String, qualities: String) -> Self {
        Self {
            name,
            sequence,
            comment: String::new(),
            qualities,
        }
    }

    /// Create a FastqRecord with a comment
    pub fn with_comment(
        name: String,
        sequence: String,
        comment: String,
        qualities: String,
    ) -> Self {
        Self {
            name,
            sequence,
            comment,
            qualities,
        }
    }

    /// Validate that the record is well-formed
    pub fn validate(&self) -> Result<()> {
        if self.sequence.len() != self.qualities.len() {
            anyhow::bail!(
                "Sequence and quality lengths don't match for read '{}': {} vs {}",
                self.name,
                self.sequence.len(),
                self.qualities.len()
            );
        }

        if self.sequence.is_empty() {
            anyhow::bail!("Empty sequence for read '{}'", self.name);
        }

        Ok(())
    }

    /// Get the length of the read
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Check if the read is empty
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }
}

/// Reader for FASTQ files (supports both plain and gzipped files)
pub struct FastqReader {
    reader: Box<dyn BufRead>,
    path: String,
}

impl FastqReader {
    /// Open a FASTQ file for reading (automatically detects gzip compression)
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path_ref = path.as_ref();
        let path_str = path_ref.to_string_lossy().to_string();
        let file = File::open(path_ref)
            .with_context(|| format!("Failed to open FASTQ file: {}", path_str))?;

        let reader: Box<dyn BufRead> = if path_str.ends_with(".gz") {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        Ok(Self {
            reader,
            path: path_str,
        })
    }

    /// Read the next FASTQ record from the file
    pub fn read_record(&mut self) -> Result<Option<FastqRecord>> {
        let mut header = String::new();
        let bytes_read = self.reader.read_line(&mut header)?;

        // End of file
        if bytes_read == 0 {
            return Ok(None);
        }

        // Parse header line (starts with '@')
        if !header.starts_with('@') {
            anyhow::bail!(
                "Expected FASTQ header starting with '@', got: {}",
                header.trim()
            );
        }

        let name = header[1..].split_whitespace().next().unwrap_or("").to_string();

        // Read sequence line
        let mut sequence = String::new();
        self.reader.read_line(&mut sequence)?;
        sequence = sequence.trim().to_string();

        // Read '+' separator line
        let mut separator = String::new();
        self.reader.read_line(&mut separator)?;
        if !separator.starts_with('+') {
            anyhow::bail!(
                "Expected '+' separator in FASTQ record for read '{}', got: {}",
                name,
                separator.trim()
            );
        }

        // Read quality line
        let mut qualities = String::new();
        self.reader.read_line(&mut qualities)?;
        qualities = qualities.trim().to_string();

        let record = FastqRecord::new(name, sequence, qualities);
        record.validate()?;

        Ok(Some(record))
    }

    /// Read all records from the file into a Vec
    pub fn read_all(&mut self) -> Result<Vec<FastqRecord>> {
        let mut records = Vec::new();
        while let Some(record) = self.read_record()? {
            records.push(record);
        }
        Ok(records)
    }

    /// Get an iterator over all records in the file
    pub fn records(self) -> FastqRecordIterator {
        FastqRecordIterator { reader: self }
    }
}

/// Iterator over FASTQ records
pub struct FastqRecordIterator {
    reader: FastqReader,
}

impl Iterator for FastqRecordIterator {
    type Item = Result<FastqRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

/// Writer for FASTQ files (supports both plain and gzipped files)
pub struct FastqWriter {
    writer: Box<dyn Write>,
    records_written: usize,
}

impl FastqWriter {
    /// Create a new FASTQ writer (automatically detects gzip compression from filename)
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path_ref = path.as_ref();
        let path_str = path_ref.to_string_lossy().to_string();
        let file = File::create(path_ref)
            .with_context(|| format!("Failed to create FASTQ file: {}", path_str))?;

        let writer: Box<dyn Write> = if path_str.ends_with(".gz") {
            Box::new(BufWriter::new(GzEncoder::new(
                file,
                Compression::default(),
            )))
        } else {
            Box::new(BufWriter::new(file))
        };

        Ok(Self {
            writer,
            records_written: 0,
        })
    }

    /// Write a single FASTQ record to the file
    pub fn write_record(&mut self, record: &FastqRecord) -> Result<()> {
        record.validate()?;

        writeln!(self.writer, "@{}", record.name)?;
        writeln!(self.writer, "{}", record.sequence)?;
        writeln!(self.writer, "+")?;
        writeln!(self.writer, "{}", record.qualities)?;

        self.records_written += 1;
        Ok(())
    }

    /// Write multiple FASTQ records to the file
    pub fn write_records(&mut self, records: &[FastqRecord]) -> Result<()> {
        for record in records {
            self.write_record(record)?;
        }
        Ok(())
    }

    /// Get the number of records written so far
    pub fn records_written(&self) -> usize {
        self.records_written
    }

    /// Flush the writer to ensure all data is written
    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush()?;
        Ok(())
    }
}

/// Paired-end FASTQ reader for reading R1 and R2 files simultaneously
pub struct PairedFastqReader {
    r1_reader: FastqReader,
    r2_reader: FastqReader,
}

impl PairedFastqReader {
    /// Create a new paired-end FASTQ reader from two file paths
    pub fn from_paths<P: AsRef<Path>>(r1_path: P, r2_path: P) -> Result<Self> {
        let r1_reader = FastqReader::from_path(r1_path)?;
        let r2_reader = FastqReader::from_path(r2_path)?;

        Ok(Self {
            r1_reader,
            r2_reader,
        })
    }

    /// Read the next pair of records (R1 and R2)
    pub fn read_pair(&mut self) -> Result<Option<(FastqRecord, FastqRecord)>> {
        let r1 = self.r1_reader.read_record()?;
        let r2 = self.r2_reader.read_record()?;

        match (r1, r2) {
            (Some(r1_record), Some(r2_record)) => Ok(Some((r1_record, r2_record))),
            (None, None) => Ok(None),
            (Some(_), None) => {
                anyhow::bail!("R1 file has more records than R2 file")
            }
            (None, Some(_)) => {
                anyhow::bail!("R2 file has more records than R1 file")
            }
        }
    }

    /// Read all paired records from both files
    pub fn read_all_pairs(&mut self) -> Result<Vec<(FastqRecord, FastqRecord)>> {
        let mut pairs = Vec::new();
        while let Some(pair) = self.read_pair()? {
            pairs.push(pair);
        }
        Ok(pairs)
    }

    /// Get an iterator over all paired records
    pub fn pairs(self) -> PairedFastqIterator {
        PairedFastqIterator { reader: self }
    }
}

/// Iterator over paired FASTQ records
pub struct PairedFastqIterator {
    reader: PairedFastqReader,
}

impl Iterator for PairedFastqIterator {
    type Item = Result<(FastqRecord, FastqRecord)>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_pair() {
            Ok(Some(pair)) => Some(Ok(pair)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

/// Paired-end FASTQ writer for writing R1 and R2 files simultaneously
pub struct PairedFastqWriter {
    r1_writer: FastqWriter,
    r2_writer: FastqWriter,
    pairs_written: usize,
}

impl PairedFastqWriter {
    /// Create a new paired-end FASTQ writer from two file paths
    pub fn from_paths<P: AsRef<Path>>(
        r1_path: P,
        r2_path: P,
        _threads: usize,
    ) -> Result<Self> {
        let r1_writer = FastqWriter::from_path(r1_path)?;
        let r2_writer = FastqWriter::from_path(r2_path)?;

        Ok(Self {
            r1_writer,
            r2_writer,
            pairs_written: 0,
        })
    }

    /// Write a pair of FASTQ records (R1 and R2)
    pub fn write_pair(&mut self, r1: &FastqRecord, r2: &FastqRecord) -> Result<()> {
        self.r1_writer.write_record(r1)?;
        self.r2_writer.write_record(r2)?;
        self.pairs_written += 1;
        Ok(())
    }

    /// Write multiple pairs of FASTQ records
    pub fn write_pairs(&mut self, pairs: &[(FastqRecord, FastqRecord)]) -> Result<()> {
        for (r1, r2) in pairs {
            self.write_pair(r1, r2)?;
        }
        Ok(())
    }

    /// Get the number of pairs written so far
    pub fn pairs_written(&self) -> usize {
        self.pairs_written
    }

    /// Flush both writers to ensure all data is written
    pub fn flush(&mut self) -> Result<()> {
        self.r1_writer.flush()?;
        self.r2_writer.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_fastq_record_creation() {
        let record = FastqRecord::new(
            "read1".to_string(),
            "ATCGATCG".to_string(),
            "IIIIIIII".to_string(),
        );

        assert_eq!(record.name, "read1");
        assert_eq!(record.sequence, "ATCGATCG");
        assert_eq!(record.qualities, "IIIIIIII");
        assert_eq!(record.len(), 8);
    }

    #[test]
    fn test_fastq_record_validation() {
        let valid = FastqRecord::new(
            "read1".to_string(),
            "ATCG".to_string(),
            "IIII".to_string(),
        );
        assert!(valid.validate().is_ok());

        let invalid = FastqRecord::new(
            "read2".to_string(),
            "ATCG".to_string(),
            "III".to_string(),
        );
        assert!(invalid.validate().is_err());

        let empty = FastqRecord::new("read3".to_string(), "".to_string(), "".to_string());
        assert!(empty.validate().is_err());
    }

    #[test]
    fn test_fastq_read_write() -> Result<()> {
        // Create a temporary FASTQ file
        let mut temp_file = NamedTempFile::new()?;
        writeln!(temp_file, "@read1")?;
        writeln!(temp_file, "ATCGATCG")?;
        writeln!(temp_file, "+")?;
        writeln!(temp_file, "IIIIIIII")?;
        writeln!(temp_file, "@read2")?;
        writeln!(temp_file, "GCTAGCTA")?;
        writeln!(temp_file, "+")?;
        writeln!(temp_file, "HHHHHHHH")?;
        temp_file.flush()?;

        // Read records
        let mut reader = FastqReader::from_path(temp_file.path())?;
        let records = reader.read_all()?;

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, "read1");
        assert_eq!(records[0].sequence, "ATCGATCG");
        assert_eq!(records[1].name, "read2");
        assert_eq!(records[1].sequence, "GCTAGCTA");

        // Write records to a new file
        let temp_out = NamedTempFile::new()?;
        let mut writer = FastqWriter::from_path(temp_out.path())?;
        writer.write_records(&records)?;
        writer.flush()?;

        // Read back and verify
        let mut reader2 = FastqReader::from_path(temp_out.path())?;
        let records2 = reader2.read_all()?;

        assert_eq!(records, records2);

        Ok(())
    }

    #[test]
    fn test_paired_fastq_operations() -> Result<()> {
        // Create temporary R1 and R2 files
        let mut r1_file = NamedTempFile::new()?;
        writeln!(r1_file, "@read1/1")?;
        writeln!(r1_file, "ATCGATCG")?;
        writeln!(r1_file, "+")?;
        writeln!(r1_file, "IIIIIIII")?;
        r1_file.flush()?;

        let mut r2_file = NamedTempFile::new()?;
        writeln!(r2_file, "@read1/2")?;
        writeln!(r2_file, "GCTAGCTA")?;
        writeln!(r2_file, "+")?;
        writeln!(r2_file, "HHHHHHHH")?;
        r2_file.flush()?;

        // Read paired records
        let mut paired_reader =
            PairedFastqReader::from_paths(r1_file.path(), r2_file.path())?;
        let pairs = paired_reader.read_all_pairs()?;

        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0].0.name, "read1/1");
        assert_eq!(pairs[0].1.name, "read1/2");

        // Write paired records
        let r1_out = NamedTempFile::new()?;
        let r2_out = NamedTempFile::new()?;
        let mut paired_writer =
            PairedFastqWriter::from_paths(r1_out.path(), r2_out.path(), 1)?;
        paired_writer.write_pairs(&pairs)?;
        paired_writer.flush()?;

        assert_eq!(paired_writer.pairs_written(), 1);

        Ok(())
    }
}
