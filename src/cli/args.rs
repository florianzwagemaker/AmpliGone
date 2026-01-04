use clap::Parser;
use clap::ValueEnum;
use std::path::PathBuf;

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum AmpliconType {
    #[value(name = "end-to-end")]
    EndToEnd,
    #[value(name = "end-to-mid")]
    EndToMid,
    #[value(name = "fragmented")]
    Fragmented,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum AlignmentPreset {
    #[value(name = "sr")]
    Sr,
    #[value(name = "map-ont")]
    MapOnt,
    #[value(name = "map-pb")]
    MapPb,
    #[value(name = "splice")]
    Splice,
}

#[derive(Debug, Clone)]
pub struct AlignmentScoring {
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open1: i32,
    pub gap_extend1: i32,
    pub gap_open2: Option<i32>,
    pub gap_extend2: Option<i32>,
    pub mma: Option<i32>,
}

impl Default for AlignmentScoring {
    fn default() -> Self {
        Self {
            match_score: 2,
            mismatch_score: 4,
            gap_open1: 4,
            gap_extend1: 2,
            gap_open2: None,
            gap_extend2: None,
            mma: None,
        }
    }
}

#[derive(Parser, Debug)]
#[command(name = "ampligone")]
#[command(author, version, about = "A tool for removing amplicon sequences from sequencing data")]
#[command(long_about = "AmpliGone: An accurate and efficient tool to remove primers from NGS reads in reference-based experiments")]
pub struct Args {
    /// Input file with reads in either FastQ or BAM format
    #[arg(short, long, value_name = "FILE")]
    pub input: PathBuf,

    /// Second input file for paired-end reads (R2)
    #[arg(long, value_name = "FILE")]
    pub input2: Option<PathBuf>,

    /// Output (FastQ) file with cleaned reads
    #[arg(short, long, value_name = "FILE")]
    pub output: PathBuf,

    /// Second output file for paired-end reads (R2)
    #[arg(long, value_name = "FILE")]
    pub output2: Option<PathBuf>,

    /// Input Reference genome in FASTA format
    #[arg(short = 'r', long = "reference", value_name = "FILE")]
    pub reference: PathBuf,

    /// Used primer sequences in FASTA format or primer coordinates in BED format.
    /// Note that using bed-files overrides error-rate and ambiguity functionality
    #[arg(short = 'p', long = "primers", value_name = "FILE")]
    pub primers: PathBuf,

    /// Define the amplicon-type, either being 'end-to-end', 'end-to-mid', or 'fragmented'
    #[arg(short = 'a', long = "amplicon-type", value_enum, default_value = "end-to-end")]
    pub amplicon_type: AmpliconType,

    /// If set, primers closely positioned to each other in the same orientation will be
    /// virtually combined into a single primer
    #[arg(long = "virtual-primers", default_value = "false")]
    pub virtual_primers: bool,

    /// The number of bases to look around a primer-site to consider it part of a fragment.
    /// Only used if amplicon-type is 'fragmented'. Default is 10
    #[arg(long = "fragment-lookaround-size", value_name = "N")]
    pub fragment_lookaround_size: Option<usize>,

    /// The maximum allowed error rate (as a percentage) for the primer search.
    /// Use 0 for exact primer matches. (0.1 = 10% error rate)
    /// Note that this is only used if the primer-file is in FASTA format
    #[arg(short = 'e', long = "error-rate", default_value = "0.1", value_name = "N")]
    pub error_rate: f64,

    /// The preset to use for alignment of reads against the reference.
    /// This can be either 'sr', 'map-ont', 'map-pb', or 'splice'
    #[arg(long = "alignment-preset", value_enum, value_name = "PRESET")]
    pub alignment_preset: Option<AlignmentPreset>,

    /// The scoring matrix to use for alignment of reads.
    /// Format: match=4 mismatch=3 gap_o1=2 gap_e1=1 [gap_o2=X gap_e2=Y mma=Z]
    #[arg(long = "alignment-scoring", value_name = "KEY=VALUE", num_args = 1..)]
    pub alignment_scoring: Option<Vec<String>>,

    /// Output BED file with found primer coordinates if they are actually cut from the reads
    #[arg(long = "export-primers", value_name = "FILE")]
    pub export_primers: Option<PathBuf>,

    /// Number of threads you wish to use
    #[arg(short = 't', long = "threads", default_value_t = num_cpus::get().max(2), value_name = "N")]
    pub threads: usize,

    /// If set, AmpliGone will always create the output files even if there is nothing to output
    /// (for example when an empty input-file is given).
    /// This is useful in (automated) pipelines
    #[arg(long = "to")]
    pub to: bool,

    /// Prints more information, like DEBUG statements, to the terminal
    #[arg(short = 'v', long = "verbose")]
    pub verbose: bool,

    /// Prints less information, like only WARNING and ERROR statements, to the terminal
    #[arg(short = 'q', long = "quiet")]
    pub quiet: bool,
}

impl Args {
    /// Validate all command-line arguments
    pub fn validate(&self) -> Result<(), String> {
        let mut errors: Vec<String> = Vec::new();

        if let Err(e) = self.validate_input_extensions() {
            errors.push(e);
        }
        if let Err(e) = self.validate_output_extensions() {
            errors.push(e);
        }
        if let Err(e) = self.validate_files_exist() {
            errors.push(e);
        }
        if let Err(e) = self.validate_paired_end_consistency() {
            errors.push(e);
        }
        if let Err(e) = self.validate_verbosity_flags() {
            errors.push(e);
        }
        if let Err(e) = self.validate_error_rate() {
            errors.push(e);
        }
        if let Err(e) = self.validate_threads() {
            errors.push(e);
        }

        if errors.is_empty() {
            Ok(())
        } else {
            return Err(errors.join("\n"))
        }
    }

    fn validate_paired_end_consistency(&self) -> Result<(), String> {
        let input2_provided = self.input2.is_some();
        let output2_provided = self.output2.is_some();

        if input2_provided != output2_provided {
            return Err(
                "Both --input2 and --output2 must be provided for paired-end mode.".to_string(),
            );
        }
        Ok(())
    }

    fn validate_input_extensions(&self) -> Result<(), String> {
        let valid_extensions = ["fastq", "fq", "fastq.gz", "fq.gz", "bam"];
        Self::validate_extension(&self.input, &valid_extensions, "input")?;
        if let Some(input2) = &self.input2 {
            Self::validate_extension(input2, &valid_extensions, "input2")?;
        }
        Ok(())
    }

    fn validate_output_extensions(&self) -> Result<(), String> {
        let valid_extensions = ["fastq", "fq", "fastq.gz", "fq.gz"];
        Self::validate_extension(&self.output, &valid_extensions, "output")?;
        if let Some(output2) = &self.output2 {
            Self::validate_extension(output2, &valid_extensions, "output2")?;
        }
        Ok(())
    }

    fn validate_extension(
        path: &PathBuf,
        valid_extensions: &[&str],
        file_type: &str,
    ) -> Result<(), String> {
        let file_str = path.to_string_lossy();
        let is_valid = valid_extensions
            .iter()
            .any(|ext| file_str.ends_with(ext));

        if !is_valid {
            return Err(format!(
                "File '{}' doesn't end with one of {:?}",
                file_str, valid_extensions
            ));
        }
        Ok(())
    }

    fn validate_files_exist(&self) -> Result<(), String> {
        Self::check_file_exists(&self.input, "Input")?;
        if let Some(input2) = &self.input2 {
            Self::check_file_exists(input2, "Input2")?;
        }
        Self::check_file_exists(&self.reference, "Reference")?;
        Self::check_file_exists(&self.primers, "Primers")?;
        Ok(())
    }

    fn check_file_exists(path: &PathBuf, file_type: &str) -> Result<(), String> {
        if !path.exists() {
            let msg = format!("{} file does not exist: {}", file_type, path.display());
            return Err(msg);
        }
        Ok(())
    }

    fn validate_verbosity_flags(&self) -> Result<(), String> {
        if self.verbose && self.quiet {
            return Err(
                "AmpliGone was given both the '--verbose' and '--quiet' flags. \
                Please only use one of these flags at a time."
                    .to_string(),
            );
        }
        Ok(())
    }

    fn validate_error_rate(&self) -> Result<(), String> {
        if !(0.0..=1.0).contains(&self.error_rate) {
            return Err(format!(
                "Error rate must be between 0.0 and 1.0, got: {}",
                self.error_rate
            ));
        }
        Ok(())
    }

    fn validate_threads(&self) -> Result<(), String> {
        if self.threads < 2 {
            return Err(
                "AmpliGone requires a minimum of 2 threads for execution.".to_string(),
            );
        }
        Ok(())
    }

    /// Parse alignment scoring parameters from command-line strings
    pub fn parse_alignment_scoring(&self) -> Result<Option<AlignmentScoring>, String> {
        let Some(scoring_args) = &self.alignment_scoring else {
            return Ok(None);
        };

        let mut scoring = AlignmentScoring::default();

        for arg in scoring_args {
            let parts: Vec<&str> = arg.split('=').collect();
            if parts.len() != 2 {
                return Err(format!("Invalid scoring parameter format: {}", arg));
            }

            let key = parts[0];
            let value: i32 = parts[1]
                .parse()
                .map_err(|_| format!("Invalid numeric value for {}: {}", key, parts[1]))?;

            match key {
                "match" => scoring.match_score = value,
                "mismatch" => scoring.mismatch_score = value,
                "gap_o1" => scoring.gap_open1 = value,
                "gap_e1" => scoring.gap_extend1 = value,
                "gap_o2" => scoring.gap_open2 = Some(value),
                "gap_e2" => scoring.gap_extend2 = Some(value),
                "mma" => scoring.mma = Some(value),
                _ => return Err(format!("Unknown scoring parameter: {}", key)),
            }
        }

        Ok(Some(scoring))
    }

    /// Get the effective fragment lookaround size based on amplicon type
    pub fn get_fragment_lookaround_size(&self) -> usize {
        match self.amplicon_type {
            AmpliconType::Fragmented => self.fragment_lookaround_size.unwrap_or(10),
            _ => 10000, // Large value for non-fragmented types
        }
    }

    /// Check if we're operating in paired-end mode
    pub fn is_paired_end(&self) -> bool {
        self.input2.is_some() && self.output2.is_some()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_validate_error_rate() {
        let mut args = create_test_args();
        args.error_rate = 0.5;
        assert!(args.validate_error_rate().is_ok());

        args.error_rate = 1.5;
        assert!(args.validate_error_rate().is_err());

        args.error_rate = -0.1;
        assert!(args.validate_error_rate().is_err());
    }

    #[test]
    fn test_validate_threads() {
        let mut args = create_test_args();
        args.threads = 2;
        assert!(args.validate_threads().is_ok());

        args.threads = 1;
        assert!(args.validate_threads().is_err());
    }

    #[test]
    fn test_validate_verbosity_flags() {
        let mut args = create_test_args();
        args.verbose = true;
        args.quiet = true;
        assert!(args.validate_verbosity_flags().is_err());

        args.quiet = false;
        assert!(args.validate_verbosity_flags().is_ok());
    }

    fn create_test_args() -> Args {
        Args {
            input: PathBuf::from("test.fastq"),
            input2: None,
            output: PathBuf::from("out.fastq"),
            output2: None,
            reference: PathBuf::from("ref.fasta"),
            primers: PathBuf::from("primers.fasta"),
            amplicon_type: AmpliconType::EndToEnd,
            virtual_primers: false,
            fragment_lookaround_size: None,
            error_rate: 0.1,
            alignment_preset: None,
            alignment_scoring: None,
            export_primers: None,
            threads: 2,
            to: false,
            verbose: false,
            quiet: false,
        }
    }
}