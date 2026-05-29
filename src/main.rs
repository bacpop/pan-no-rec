use anyhow::Error;

// Delegates process startup to the CLI module.
fn main() -> Result<(), Error> {
    pan_no_rec::main()
}
