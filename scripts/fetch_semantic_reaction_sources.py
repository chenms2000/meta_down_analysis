#!/usr/bin/env python
"""Download optional semantic reaction sources for database_accuracy_store.v2.

The files cached by this script are used to increase real substrate/product
coverage. They are optional at service startup: the v2 builder parses them when
present and falls back to Rhea TSV/Reactome PE mappings when absent.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import urllib.request
from pathlib import Path
from typing import Any


DEFAULT_SOURCES = {
    "rhea_biopax": {
        "source_name": "Rhea BioPAX level 3",
        "source_url": "https://ftp.expasy.org/databases/rhea/biopax/rhea-biopax.owl.gz",
        "relative_path": "raw_lake/rhea/{release_id}/rhea-biopax.owl.gz",
        "license_id": "CC-BY-4.0",
    },
    "rhea_rdf": {
        "source_name": "Rhea RDF",
        "source_url": "https://ftp.expasy.org/databases/rhea/rdf/rhea.rdf.gz",
        "relative_path": "raw_lake/rhea/{release_id}/rhea.rdf.gz",
        "license_id": "CC-BY-4.0",
    },
    "reactome_sbml": {
        "source_name": "Reactome human reactions SBML",
        "source_url": "https://reactome.org/download/current/homo_sapiens.3.1.sbml.tgz",
        "relative_path": "raw_lake/reactome/{release_id}/homo_sapiens.3.1.sbml.tgz",
        "license_id": "Reactome Terms/CC",
    },
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parser_hash() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]


def download_file(url: str, output: Path, *, force: bool = False) -> dict[str, Any]:
    output.parent.mkdir(parents=True, exist_ok=True)
    status = "cached"
    if force or not output.exists():
        request = urllib.request.Request(url, headers={"User-Agent": "metabo-data-end/semantic-source-fetcher"})
        tmp = output.with_suffix(output.suffix + ".tmp")
        with urllib.request.urlopen(request, timeout=180) as response, tmp.open("wb") as handle:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                handle.write(chunk)
        tmp.replace(output)
        status = "downloaded"
    return {
        "path": str(output),
        "status": status,
        "bytes": output.stat().st_size,
        "checksum": sha256_file(output),
    }


def fetch_sources(workspace: Path, release_id: str, *, force: bool = False, include_reactome_sbml: bool = False) -> Path:
    records = []
    selected = ["rhea_biopax", "rhea_rdf"]
    if include_reactome_sbml:
        selected.append("reactome_sbml")
    for key in selected:
        spec = DEFAULT_SOURCES[key]
        output = workspace / spec["relative_path"].format(release_id=release_id)
        result = download_file(spec["source_url"], output, force=force)
        records.append(
            {
                "source_key": key,
                "source_name": spec["source_name"],
                "source_url": spec["source_url"],
                "download_date": dt.datetime.now(dt.UTC).replace(microsecond=0).isoformat(),
                "checksum": result["checksum"],
                "license_id": spec["license_id"],
                "parser_hash": parser_hash(),
                "file_path": result["path"],
                "bytes": result["bytes"],
                "status": result["status"],
            }
        )
    manifest = {
        "manifest_version": "semantic_reaction_sources.v1",
        "release_id": release_id,
        "workspace": str(workspace),
        "sources": records,
    }
    manifest_path = workspace / "raw_lake" / "semantic_reaction_sources_manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch optional RDF/BioPAX/SBML semantic reaction sources.")
    parser.add_argument("--workspace", type=Path, default=Path("."))
    parser.add_argument("--release-id", required=True)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--include-reactome-sbml", action="store_true", help="Also download the large Reactome human SBML tgz.")
    args = parser.parse_args()
    manifest = fetch_sources(
        args.workspace.resolve(),
        args.release_id,
        force=args.force,
        include_reactome_sbml=args.include_reactome_sbml,
    )
    print(json.dumps({"semantic_reaction_sources_manifest": str(manifest)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
