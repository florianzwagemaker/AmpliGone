mod cli;
mod io;
mod cutting;
mod alignment;
mod utils;
mod primer;

use clap::Parser;
use cli::args::Args;
use io::{PairedFastqReader, PairedFastqWriter};
use tracing::{error, info};
use utils::logging::init_logging;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Parse CLI arguments using the Args struct defined in src/cli/args.rs
    let args = Args::parse();

    // Initialize logging based on verbosity flags
    init_logging(args.verbose, args.quiet);

    // Validate all arguments
    if let Err(e) = args.validate() {
        error!("Argument validation error(s): {}", e);
        std::process::exit(1);
    }

    // Dispatch to appropriate runner based on input files
    if args.is_paired_end() {
        run_paired_end(&args)?;
    } else {
        run_single_end(&args)?;
    }

    Ok(())
}

/// Handles processing for paired-end sequencing data (R1 and R2)
fn run_paired_end(args: &Args) -> Result<(), Box<dyn std::error::Error>> {
    let input2 = args.input2.as_ref().unwrap();
    let output2 = args.output2.as_ref().unwrap();

    info!("Running in paired-end mode");
    info!("R1 input: {}", args.input.display());
    info!("R2 input: {}", input2.display());

    // Create paired reader from src/paired/reader.rs
    let reader = PairedFastqReader::from_paths(&args.input, input2)?;

    // Create paired writer
    let mut writer = PairedFastqWriter::from_paths(&args.output, output2, args.threads)?;

    // TODO: Initialize the core cutting logic
    // let cutter = create_cutter(args)?;
    
    // TODO: Create the paired-end processor
    // let processor = PairedEndProcessor::new(cutter);

    // TODO: Execute processing loop with a chunk size of 10,000 for parallel efficiency
    // let stats = processor.process_and_write(reader, &mut writer, 10000)?;

    info!("Paired-end processing not yet fully implemented");

    Ok(())
}

/// Handles processing for single-end sequencing data
fn run_single_end(args: &Args) -> Result<(), Box<dyn std::error::Error>> {
    info!("Running in single-end mode");
    info!("Input: {}", args.input.display());
    
    // TODO: Single-end implementation would follow a similar pattern:
    // 1. Create FastqReader
    // 2. Create FastqWriter
    // 3. Create Cutter/Processor
    // 4. Process and write
    
    info!("Single-end processing not yet fully implemented");
    
    Ok(())
}