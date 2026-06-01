use std::fs;
use std::path::Path;
use std::process::{Command, Output};
use tempfile::TempDir;

// Writes a FASTA alignment fixture into the temporary directory.
fn write_alignment(dir: &TempDir, name: &str, contents: &str) {
    let alignment_dir = dir.path().join("aligned_gene_sequences");
    fs::create_dir_all(&alignment_dir).unwrap();
    fs::write(alignment_dir.join(name), contents).unwrap();
}

// Writes the required Panaroo Rtab header fixture.
fn write_rtab(dir: &TempDir, sample_names: &[&str]) {
    fs::write(
        dir.path().join("gene_presence_absence.Rtab"),
        format!("Gene\t{}\n", sample_names.join("\t")),
    )
    .unwrap();
}

#[test]
// Verifies CLI TSV output is stable across thread counts.
fn cli_prints_presence_absence_tsv_stably_across_thread_counts() {
    let dir = TempDir::new().unwrap();
    write_rtab(&dir, &["alpha", "beta"]);
    write_alignment(
        &dir,
        "gene_low.aln.fas",
        ">beta\nAAAAAAAAAA\n>alpha\nAAAAAAAAAA\n",
    );
    write_alignment(
        &dir,
        "gene_middle.aln.fas",
        ">beta\nCAAAAAAAAA\n>alpha\nAAAAAAAAAA\n",
    );
    write_alignment(
        &dir,
        "gene_high.aln.fas",
        ">beta\nCCCCCCCCAA\n>alpha\nAAAAAAAAAA\n",
    );

    let one_thread = run_cli(dir.path(), "1");
    let two_threads = run_cli(dir.path(), "2");

    let expected = concat!(
        "gene\talpha\tbeta\n",
        "gene_high\t0\t0\n",
        "gene_low\t0\t0\n",
        "gene_middle\t0\t0\n"
    );
    assert_eq!(one_thread, expected);
    assert_eq!(two_threads, expected);
}

#[test]
// Verifies Panaroo directory input discovers standard aligned_gene_sequences genes.
fn cli_accepts_panaroo_dir_aligned_gene_sequences() {
    let dir = TempDir::new().unwrap();
    write_rtab(&dir, &["alpha", "beta"]);
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
    fs::write(
        dir.path().join("ignored.aln.fas"),
        ">alpha\nAAAA\n>beta\nAAAA\n",
    )
    .unwrap();

    let observed = run_cli_panaroo_dir(dir.path());

    let expected = concat!(
        "gene\talpha\tbeta\n",
        "gene_high\t0\t0\n",
        "gene_low\t0\t0\n"
    );
    assert_eq!(observed, expected);
}

#[test]
// Verifies Panaroo directory input requires aligned_gene_sequences.
fn cli_rejects_panaroo_dir_without_aligned_gene_sequences() {
    let dir = TempDir::new().unwrap();
    write_rtab(&dir, &["alpha", "beta"]);
    fs::write(
        dir.path().join("gene_low.aln.fas"),
        ">beta\nAAAAAAAAAA\n>alpha\nAAAAAAAAAA\n",
    )
    .unwrap();

    let output = panaroo_command(dir.path())
        .arg("--threads")
        .arg("1")
        .output()
        .unwrap();

    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("aligned_gene_sequences"),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
// Verifies first and longest paralog modes can change selected duplicate records.
fn cli_paralog_mode_first_and_longest_affect_output() {
    let dir = TempDir::new().unwrap();
    write_rtab(&dir, &["alpha", "beta", "delta", "gamma"]);
    write_alignment(
        &dir,
        "gene_low.aln.fas",
        concat!(
            ">alpha\nAAAAAAAA--\n",
            ">beta\nAAAAAAAA--\n",
            ">gamma\nAAAAAAAA--\n",
            ">delta\nAAAAAAAA--\n",
        ),
    );
    write_alignment(
        &dir,
        "gene_dup.aln.fas",
        concat!(
            ">alpha;first\nAAAAAAAA--\n",
            ">beta\nAAAAAAAA--\n",
            ">gamma\nAAAAAAAA--\n",
            ">delta\nAAAAAAAA--\n",
            ">alpha;second\nCCCCCCCCAA\n",
        ),
    );

    let first = stdout_from_success(run_cli_with_paralog_mode(dir.path(), "first"));
    let longest = stdout_from_success(run_cli_with_paralog_mode(dir.path(), "longest"));

    let expected_first = concat!(
        "gene\talpha\tbeta\tdelta\tgamma\n",
        "gene_dup\t0\t0\t0\t0\n",
        "gene_low\t0\t0\t0\t0\n"
    );
    let expected_longest = concat!(
        "gene\talpha\tbeta\tdelta\tgamma\n",
        "gene_dup\t1\t1\t1\t1\n",
        "gene_low\t0\t0\t0\t0\n"
    );
    assert_eq!(first, expected_first);
    assert_eq!(longest, expected_longest);
}

#[test]
// Verifies skip mode retains paralogous genes after removing duplicated samples.
fn cli_paralog_mode_skip_retains_paralogous_alignments() {
    let dir = TempDir::new().unwrap();
    write_rtab(&dir, &["alpha", "beta"]);
    write_alignment(
        &dir,
        "gene_dup.aln.fas",
        ">alpha;first\nAAAA\n>beta\nAAAA\n>alpha;second\nCCCC\n",
    );
    write_alignment(&dir, "gene_clean.aln.fas", ">alpha\nAAAA\n>beta\nAAAA\n");

    let observed = stdout_from_success(run_cli_with_paralog_mode(dir.path(), "skip"));

    assert_eq!(
        observed,
        "gene\talpha\tbeta\ngene_clean\t0\t0\ngene_dup\t0\t0\n"
    );
}

#[test]
// Verifies clap rejects unsupported paralog mode values.
fn cli_rejects_invalid_paralog_mode() {
    let output = Command::new(env!("CARGO_BIN_EXE_pan-no-rec"))
        .arg("--panaroo-dir")
        .arg("unused")
        .arg("--paralog-mode")
        .arg("invalid")
        .output()
        .unwrap();

    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("invalid value"),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
// Verifies duplicate-containing alignments emit one summary warning and report.
fn cli_summarizes_paralog_warnings_and_writes_default_report() {
    let dir = TempDir::new().unwrap();
    write_rtab(&dir, &["alpha", "beta"]);
    write_alignment(
        &dir,
        "gene_dup_a.aln.fas",
        ">alpha;first\nAAAA\n>beta\nAAAA\n>alpha;second\nCCCC\n",
    );
    write_alignment(&dir, "gene_clean.aln.fas", ">alpha\nAAAA\n>beta\nAAAA\n");
    write_alignment(
        &dir,
        "gene_dup_b.aln.fas",
        concat!(
            ">alpha;first\nAAAA\n",
            ">beta;first\nAAAA\n",
            ">alpha;second\nCCCC\n",
            ">beta;second\nCCCC\n",
        ),
    );

    let output = run_cli_with_paralog_mode(dir.path(), "first");

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(stderr.matches("contained paralogs").count(), 1);
    assert!(
        stderr.contains("2 alignments contained paralogs; wrote paralog report to 'paralogs.txt'; using paralog mode 'first'"),
        "stderr: {stderr}"
    );
    assert!(
        !stderr.contains("alignment 'gene_dup_a' contains paralogs"),
        "stderr: {stderr}"
    );
    assert_eq!(
        fs::read_to_string(dir.path().join("paralogs.txt")).unwrap(),
        "gene\tparalogs\ngene_dup_a\t1\ngene_dup_b\t2\n"
    );
}

#[test]
// Verifies --paralog-report selects the report path.
fn cli_writes_custom_paralog_report_path() {
    let dir = TempDir::new().unwrap();
    write_rtab(&dir, &["alpha", "beta"]);
    write_alignment(
        &dir,
        "gene_dup.aln.fas",
        ">alpha;first\nAAAA\n>beta\nAAAA\n>alpha;second\nCCCC\n",
    );
    let report_path = dir.path().join("custom-paralogs.tsv");

    let output = panaroo_command(dir.path())
        .arg("--threads")
        .arg("1")
        .arg("--paralog-mode")
        .arg("first")
        .arg("--paralog-report")
        .arg(&report_path)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        fs::read_to_string(report_path).unwrap(),
        "gene\tparalogs\ngene_dup\t1\n"
    );
}

#[test]
// Verifies no report is created when no paralogs are present.
fn cli_no_paralogs_does_not_create_report() {
    let dir = TempDir::new().unwrap();
    write_rtab(&dir, &["alpha", "beta"]);
    write_alignment(&dir, "gene_clean.aln.fas", ">alpha\nAAAA\n>beta\nAAAA\n");
    let report_path = dir.path().join("no-paralogs.tsv");

    let output = panaroo_command(dir.path())
        .arg("--threads")
        .arg("1")
        .arg("--paralog-report")
        .arg(&report_path)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(!report_path.exists());
    assert!(
        !String::from_utf8_lossy(&output.stderr).contains("contained paralogs"),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
// Verifies --quiet suppresses the summary warning but still writes the report.
fn cli_quiet_suppresses_paralog_warning_but_writes_report() {
    let dir = TempDir::new().unwrap();
    write_rtab(&dir, &["alpha", "beta"]);
    write_alignment(
        &dir,
        "gene_dup.aln.fas",
        ">alpha;first\nAAAA\n>beta\nAAAA\n>alpha;second\nCCCC\n",
    );
    let report_path = dir.path().join("quiet-paralogs.tsv");

    let output = panaroo_command(dir.path())
        .arg("--threads")
        .arg("1")
        .arg("--quiet")
        .arg("--paralog-report")
        .arg(&report_path)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        !String::from_utf8_lossy(&output.stderr).contains("contained paralogs"),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        fs::read_to_string(report_path).unwrap(),
        "gene\tparalogs\ngene_dup\t1\n"
    );
}

#[test]
// Verifies clap requires Panaroo input.
fn cli_rejects_missing_input_source() {
    let output = Command::new(env!("CARGO_BIN_EXE_pan-no-rec"))
        .arg("--threads")
        .arg("1")
        .output()
        .unwrap();

    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("--panaroo-dir"),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
// Verifies gene_presence_absence.Rtab controls sample column order.
fn cli_uses_rtab_header_sample_order() {
    let dir = TempDir::new().unwrap();
    write_rtab(&dir, &["beta", "alpha"]);
    write_alignment(&dir, "gene_clean.aln.fas", ">alpha\nAAAA\n>beta\nAAAA\n");

    let observed = run_cli_panaroo_dir(dir.path());

    assert_eq!(observed, "gene\tbeta\talpha\ngene_clean\t0\t0\n");
}

#[test]
// Verifies --max-entropy removes only alignments above the threshold.
fn cli_max_entropy_filters_high_entropy_genes_and_logs_summary() {
    let dir = TempDir::new().unwrap();
    write_rtab(&dir, &["alpha", "beta"]);
    write_alignment(&dir, "gene_low.aln.fas", ">alpha\nAAAA\n>beta\nAAAA\n");
    write_alignment(&dir, "gene_equal.aln.fas", ">alpha\nAAAA\n>beta\nAAAA\n");
    write_alignment(&dir, "gene_high.aln.fas", ">alpha\nAAAA\n>beta\nCCCC\n");
    fs::write(
        dir.path().join("alignment_entropy.csv"),
        "gene,entropy\ngene_low,0.1\ngene_equal,0.5\ngene_high,0.5001\n",
    )
    .unwrap();

    let output = panaroo_command(dir.path())
        .arg("--threads")
        .arg("1")
        .arg("--max-entropy")
        .arg("0.5")
        .arg("--verbose")
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        "gene\talpha\tbeta\ngene_equal\t0\t0\ngene_low\t0\t0\n"
    );
    assert!(
        String::from_utf8_lossy(&output.stderr)
            .contains("Filtered 1 alignments with entropy > 0.5"),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
// Verifies --max-entropy requires Panaroo's entropy CSV.
fn cli_max_entropy_errors_when_entropy_csv_is_missing() {
    let dir = TempDir::new().unwrap();
    write_rtab(&dir, &["alpha", "beta"]);
    write_alignment(&dir, "gene_low.aln.fas", ">alpha\nAAAA\n>beta\nAAAA\n");

    let output = panaroo_command(dir.path())
        .arg("--threads")
        .arg("1")
        .arg("--max-entropy")
        .arg("0.5")
        .output()
        .unwrap();

    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("alignment_entropy.csv"),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
// Verifies genes absent from alignment_entropy.csv are kept with one warning.
fn cli_max_entropy_keeps_missing_entropy_rows_and_warns() {
    let dir = TempDir::new().unwrap();
    write_rtab(&dir, &["alpha", "beta"]);
    write_alignment(&dir, "gene_known.aln.fas", ">alpha\nAAAA\n>beta\nAAAA\n");
    write_alignment(&dir, "gene_missing.aln.fas", ">alpha\nAAAA\n>beta\nAAAA\n");
    fs::write(
        dir.path().join("alignment_entropy.csv"),
        "gene,entropy\ngene_known,0.1\n",
    )
    .unwrap();

    let output = panaroo_command(dir.path())
        .arg("--threads")
        .arg("1")
        .arg("--max-entropy")
        .arg("0.5")
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        "gene\talpha\tbeta\ngene_known\t0\t0\ngene_missing\t0\t0\n"
    );
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("1 alignments lacked entropy metadata"),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
// Verifies clap rejects invalid thread counts before inference starts.
fn cli_rejects_zero_threads() {
    let output = Command::new(env!("CARGO_BIN_EXE_pan-no-rec"))
        .arg("--panaroo-dir")
        .arg("unused")
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

// Runs the binary with one Panaroo directory and thread count.
fn run_cli(path: &Path, threads: &str) -> String {
    let output = panaroo_command(path)
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

// Runs the binary with one Panaroo directory and selected paralog mode.
fn run_cli_with_paralog_mode(path: &Path, mode: &str) -> Output {
    panaroo_command(path)
        .arg("--threads")
        .arg("1")
        .arg("--paralog-mode")
        .arg(mode)
        .output()
        .unwrap()
}

// Creates a command for a Panaroo directory run in that directory.
fn panaroo_command(path: &Path) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_pan-no-rec"));
    command.current_dir(path);
    command.arg("--panaroo-dir").arg(path);
    command
}

// Returns stdout for successful command output.
fn stdout_from_success(output: Output) -> String {
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    String::from_utf8(output.stdout).unwrap()
}

// Runs the binary with one Panaroo directory.
fn run_cli_panaroo_dir(path: &Path) -> String {
    let output = panaroo_command(path)
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
