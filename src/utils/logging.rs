use tracing::Level;
use tracing_subscriber::{fmt, EnvFilter};

/// Initialize logging based on verbosity flags
///
/// # Arguments
///
/// * `verbose` - If true, sets log level to DEBUG
/// * `quiet` - If true, sets log level to WARN
///
/// If neither flag is set, the default level is INFO.
pub fn init_logging(verbose: bool, quiet: bool) {
    let level = if verbose {
        Level::DEBUG
    } else if quiet {
        Level::WARN
    } else {
        Level::INFO
    };

    let filter = EnvFilter::from_default_env()
        .add_directive(level.into());

    fmt()
        .with_env_filter(filter)
        .with_target(false)
        .init();
}
