from pathlib import Path
import shutil
import tarfile
import unittest

import pangenome_recombination


TESTS_DIR = Path(__file__).resolve().parent
ARCHIVE = TESTS_DIR / "gene_alignments.tar.bz2"
EXTRACTED_DIR = TESTS_DIR / "gene_alignments"


def alignment_paths(directory):
    return sorted(
        path
        for path in directory.rglob("*.aln.fas")
        if not path.name.startswith("._")
        and has_unique_normalized_sample_ids(path)
    )


def has_unique_normalized_sample_ids(path):
    seen = set()
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            record_id = line[1:].split(None, 1)[0]
            sample_id = record_id.split(";", 1)[0]
            if sample_id in seen:
                return False
            seen.add(sample_id)
    return bool(seen)


def ensure_gene_alignments():
    paths = alignment_paths(EXTRACTED_DIR)
    if paths:
        return paths

    if not ARCHIVE.is_file():
        raise FileNotFoundError(f"missing fixture archive: {ARCHIVE}")

    EXTRACTED_DIR.mkdir(parents=True, exist_ok=True)
    with tarfile.open(ARCHIVE, "r:bz2") as archive:
        for member in archive.getmembers():
            member_path = Path(member.name)
            if (
                member.isdir()
                or member_path.is_absolute()
                or ".." in member_path.parts
                or member_path.name.startswith("._")
                or not member_path.name.endswith(".aln.fas")
            ):
                continue

            target = EXTRACTED_DIR / member_path
            target.parent.mkdir(parents=True, exist_ok=True)
            source = archive.extractfile(member)
            if source is None:
                continue
            with source, target.open("wb") as output:
                shutil.copyfileobj(source, output)

    paths = alignment_paths(EXTRACTED_DIR)
    if not paths:
        raise AssertionError(f"no usable *.aln.fas files extracted from {ARCHIVE}")
    return paths


class GeneAlignmentSmokeTest(unittest.TestCase):
    def test_compare_alignments_fixture(self):
        paths = [str(path) for path in ensure_gene_alignments()]

        result = pangenome_recombination.compare_alignments(paths, threads=2)

        self.assertIsInstance(result, dict)
        for gene, pairs in result.items():
            self.assertIsInstance(gene, str)
            self.assertIsInstance(pairs, list)
            for pair in pairs:
                self.assertIsInstance(pair, tuple)
                self.assertEqual(len(pair), 2)
                self.assertIsInstance(pair[0], str)
                self.assertIsInstance(pair[1], str)


if __name__ == "__main__":
    unittest.main()
