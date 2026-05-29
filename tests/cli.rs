use std::fs;
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

// Writes a FASTA alignment fixture into the temporary directory.
fn write_alignment(dir: &TempDir, name: &str, contents: &str) {
    fs::write(dir.path().join(name), contents).unwrap();
}

#[test]
// Verifies CLI TSV output is stable across thread counts.
fn cli_prints_presence_absence_tsv_stably_across_thread_counts() {
    let dir = TempDir::new().unwrap();
    write_alignment(
        &dir,
        "gene_low.aln",
        ">beta\nAAAAAAAAAA\n>alpha\nAAAAAAAAAA\n",
    );
    write_alignment(
        &dir,
        "gene_middle.aln",
        ">beta\nCAAAAAAAAA\n>alpha\nAAAAAAAAAA\n",
    );
    write_alignment(
        &dir,
        "gene_high.aln",
        ">beta\nCCCCCCCCAA\n>alpha\nAAAAAAAAAA\n",
    );
    let list_path = dir.path().join("msa-list.txt");
    fs::write(
        &list_path,
        "\n# input alignments\ngene_low.aln\ngene_middle.aln\ngene_high.aln\n",
    )
    .unwrap();

    let one_thread = run_cli(&list_path, "1");
    let two_threads = run_cli(&list_path, "2");

    let expected = concat!(
        "gene\talpha\tbeta\n",
        "gene_low\t0\t0\n",
        "gene_middle\t0\t0\n",
        "gene_high\t0\t0\n"
    );
    assert_eq!(one_thread, expected);
    assert_eq!(two_threads, expected);
}

#[test]
// Verifies Panaroo directory input discovers top-level .aln.fas genes.
fn cli_accepts_panaroo_dir_top_level_alignments() {
    let dir = TempDir::new().unwrap();
    write_alignment(
        &dir,
        "gene_low.aln.fas",
        ">beta\nAAAAAAAAAA\n>alpha\nAAAAAAAAAA\n",
    );
    write_alignment(
        &dir,
        "gene_high.aln.fas",
        ">beta\nCCCCCCCCAA\n>alpha\nAAAAAAAAAA\n",
    );
    write_alignment(&dir, "ignored.fas", ">alpha\nAAAA\n>beta\nAAAA\n");

    let observed = run_cli_panaroo_dir(dir.path());

    let expected = concat!(
        "gene\talpha\tbeta\n",
        "gene_high\t0\t0\n",
        "gene_low\t0\t0\n"
    );
    assert_eq!(observed, expected);
}

#[test]
// Verifies Panaroo directory input falls back to aligned_gene_sequences.
fn cli_accepts_panaroo_dir_fallback_alignments() {
    let dir = TempDir::new().unwrap();
    fs::write(
        dir.path().join("ignored.fas"),
        ">alpha\nAAAA\n>beta\nAAAA\n",
    )
    .unwrap();
    let fallback = dir.path().join("aligned_gene_sequences");
    fs::create_dir(&fallback).unwrap();
    fs::write(
        fallback.join("gene_low.aln.fas"),
        ">beta\nAAAAAAAAAA\n>alpha\nAAAAAAAAAA\n",
    )
    .unwrap();

    let observed = run_cli_panaroo_dir(dir.path());

    let expected = concat!("gene\talpha\tbeta\n", "gene_low\t0\t0\n");
    assert_eq!(observed, expected);
}

#[test]
// Verifies clap requires one input source.
fn cli_rejects_missing_input_source() {
    let output = Command::new(env!("CARGO_BIN_EXE_pan-no-rec"))
        .arg("--threads")
        .arg("1")
        .output()
        .unwrap();

    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("--msa-list"),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("--panaroo-dir"),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
// Verifies clap rejects simultaneous input sources.
fn cli_rejects_multiple_input_sources() {
    let dir = TempDir::new().unwrap();
    let list_path = dir.path().join("msa-list.txt");
    fs::write(&list_path, "").unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_pan-no-rec"))
        .arg("--msa-list")
        .arg(&list_path)
        .arg("--panaroo-dir")
        .arg(dir.path())
        .output()
        .unwrap();

    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("cannot be used with"),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
// Verifies clap rejects invalid thread counts before inference starts.
fn cli_rejects_zero_threads() {
    let output = Command::new(env!("CARGO_BIN_EXE_pan-no-rec"))
        .arg("--msa-list")
        .arg("unused.txt")
        .arg("--threads")
        .arg("0")
        .output()
        .unwrap();

    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("thread count must be at least 1"),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

// Runs the binary with one MSA list and thread count.
fn run_cli(list_path: &Path, threads: &str) -> String {
    let output = Command::new(env!("CARGO_BIN_EXE_pan-no-rec"))
        .arg("--msa-list")
        .arg(list_path)
        .arg("--threads")
        .arg(threads)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    String::from_utf8(output.stdout).unwrap()
}

// Runs the binary with one Panaroo directory.
fn run_cli_panaroo_dir(path: &Path) -> String {
    let output = Command::new(env!("CARGO_BIN_EXE_pan-no-rec"))
        .arg("--panaroo-dir")
        .arg(path)
        .arg("--threads")
        .arg("1")
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    String::from_utf8(output.stdout).unwrap()
}
