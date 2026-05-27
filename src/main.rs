use clap::Parser;
use pangenome_recombination::{CompareError, compare_alignments};
use std::error::Error;
use std::fmt;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

#[derive(Debug, Parser)]
#[command(version, about = "Compare pangenome multiple-sequence alignments")]
struct Args {
    #[arg(long, value_name = "PATH")]
    msa_list: PathBuf,

    #[arg(long, default_value_t = 1, value_parser = parse_threads)]
    threads: usize,
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("{error}");
            ExitCode::FAILURE
        }
    }
}

fn run() -> Result<(), CliError> {
    let args = Args::parse();
    let aln_paths = read_msa_list(&args.msa_list)?;
    let _hits = compare_alignments(&aln_paths, args.threads)?;
    Ok(())
}

fn read_msa_list(path: &Path) -> Result<Vec<PathBuf>, CliError> {
    let contents = fs::read_to_string(path).map_err(|source| CliError::ReadMsaList {
        path: path.to_path_buf(),
        source,
    })?;
    Ok(parse_msa_list(path, &contents))
}

fn parse_msa_list(list_path: &Path, contents: &str) -> Vec<PathBuf> {
    let base_dir = list_path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));

    contents
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(|line| {
            let path = PathBuf::from(line);
            if path.is_absolute() {
                path
            } else {
                base_dir.join(path)
            }
        })
        .collect()
}

fn parse_threads(value: &str) -> Result<usize, String> {
    let threads = value
        .parse::<usize>()
        .map_err(|source| format!("invalid thread count: {source}"))?;

    if threads == 0 {
        Err("thread count must be at least 1".to_string())
    } else {
        Ok(threads)
    }
}

#[derive(Debug)]
enum CliError {
    ReadMsaList {
        path: PathBuf,
        source: std::io::Error,
    },
    Compare(CompareError),
}

impl fmt::Display for CliError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CliError::ReadMsaList { path, source } => {
                write!(f, "failed to read MSA list '{}': {source}", path.display())
            }
            CliError::Compare(error) => write!(f, "{error}"),
        }
    }
}

impl Error for CliError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            CliError::ReadMsaList { source, .. } => Some(source),
            CliError::Compare(error) => Some(error),
        }
    }
}

impl From<CompareError> for CliError {
    fn from(error: CompareError) -> Self {
        CliError::Compare(error)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_msa_list_ignores_blank_lines_and_comments() {
        let observed = parse_msa_list(
            Path::new("fixtures/list.txt"),
            "\n  \n# comment\n  # indented comment\nalpha.aln\n\nbeta.aln\n",
        );

        assert_eq!(
            observed,
            vec![
                PathBuf::from("fixtures/alpha.aln"),
                PathBuf::from("fixtures/beta.aln")
            ]
        );
    }

    #[test]
    fn parse_msa_list_resolves_relative_paths_against_list_directory() {
        let observed = parse_msa_list(
            Path::new("data/lists/msa.txt"),
            "../gene.aln\nnested/gene.aln",
        );

        assert_eq!(
            observed,
            vec![
                PathBuf::from("data/lists/../gene.aln"),
                PathBuf::from("data/lists/nested/gene.aln")
            ]
        );
    }

    #[test]
    fn parse_msa_list_preserves_absolute_paths() {
        let absolute = std::env::current_dir().unwrap().join("gene.aln");
        let observed = parse_msa_list(Path::new("data/list.txt"), &absolute.to_string_lossy());

        assert_eq!(observed, vec![absolute]);
    }
}
