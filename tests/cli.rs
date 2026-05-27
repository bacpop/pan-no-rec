use std::fs;
use std::process::Command;
use tempfile::TempDir;

fn write_alignment(dir: &TempDir, name: &str, contents: &str) {
    fs::write(dir.path().join(name), contents).unwrap();
}

#[test]
fn cli_runs_comparison_from_msa_list_without_stdout() {
    let dir = TempDir::new().unwrap();
    write_alignment(&dir, "gene_low.aln", ">s1\nAAAAAAAAAA\n>s2\nAAAAAAAAAA\n");
    write_alignment(
        &dir,
        "gene_middle.aln",
        ">s1\nAAAAAAAAAA\n>s2\nCAAAAAAAAA\n",
    );
    write_alignment(&dir, "gene_high.aln", ">s1\nAAAAAAAAAA\n>s2\nCCCCCCCCAA\n");
    let list_path = dir.path().join("msa-list.txt");
    fs::write(
        &list_path,
        "\n# input alignments\ngene_low.aln\ngene_middle.aln\ngene_high.aln\n",
    )
    .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_pangenome_recombination"))
        .arg("--msa-list")
        .arg(&list_path)
        .arg("--threads")
        .arg("2")
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stdout.is_empty());
}
