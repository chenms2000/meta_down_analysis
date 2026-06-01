import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "download_databases.py"
SPEC = importlib.util.spec_from_file_location("download_databases", SCRIPT)
download_databases = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = download_databases
SPEC.loader.exec_module(download_databases)


class DownloadCatalogTests(unittest.TestCase):
    def test_source_selection_blocks_licensed_without_acceptance(self):
        source = {
            "id": "hmdb",
            "layer": "licensed_enhancement",
            "enabled_by_default": False,
        }
        selected, reason = download_databases.source_selected(
            source,
            requested_sources={"hmdb"},
            requested_layers={"licensed_enhancement"},
            include_disabled=True,
            accept_licensed=False,
        )
        self.assertFalse(selected)
        self.assertIn("requires --accept-licensed", reason)

    def test_html_link_parser_resolves_relative_links(self):
        parser = download_databases.LinkParser("https://example.org/data/")
        parser.feed('<a href="file.tsv.gz"> File TSV </a><a href="../">Parent Directory</a>')
        self.assertEqual(parser.links[0].url, "https://example.org/data/file.tsv.gz")
        self.assertEqual(parser.links[0].text, "File TSV")

    def test_relative_path_uses_anchor_text_for_download_links(self):
        rel = download_databases.relative_from_base(
            "https://example.org/databases/lmsd/download?file=abc",
            "https://example.org/databases/lmsd/download",
            text="LMSD 2026-05-11 (ZIP)",
        )
        self.assertEqual(rel, "LMSD_2026-05-11_(ZIP)")

    def test_relative_path_can_prefer_anchor_text(self):
        rel = download_databases.relative_from_base(
            "https://ndownloader.figshare.com/files/12345",
            "https://www.bridgedb.org/mapping-databases/metabolite-mappings.html",
            text="metabolites_20200809.bridge",
            prefer_text=True,
        )
        self.assertEqual(rel, "metabolites_20200809.bridge")

    def test_manifest_write_round_trip(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = {"release_id": "test", "files": []}
            path = download_databases.write_manifest(manifest, Path(tmp), "test")
            self.assertTrue(path.exists())
            self.assertIn('"release_id": "test"', path.read_text(encoding="utf-8"))

    def test_ensure_within_accepts_child_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            child = root / "source" / "release" / "file.txt"
            download_databases.ensure_within(root, child)

    def test_index_path_prefix_keeps_dataset_partitions_apart(self):
        spec = {
            "url": "https://example.org/output/target/",
            "path_prefix": "target",
            "include_regex": r"\.parquet$",
        }
        original = download_databases.parse_index_links
        try:
            download_databases.parse_index_links = lambda *args, **kwargs: [
                download_databases.Link("https://example.org/output/target/part-000.parquet", "part-000.parquet")
            ]
            files = download_databases.resolve_index_file_spec("opentargets_core", spec, 1, 1, True)
        finally:
            download_databases.parse_index_links = original
        self.assertEqual(files[0].relative_path, "target/part-000.parquet")


if __name__ == "__main__":
    unittest.main()
