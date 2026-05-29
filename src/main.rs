mod cli;

use std::process::ExitCode;

// Delegates process startup to the CLI module.
fn main() -> ExitCode {
    cli::launch()
}
