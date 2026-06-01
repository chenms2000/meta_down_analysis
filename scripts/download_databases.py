#!/usr/bin/env python3
"""Versioned downloader for metabolism-oncology source databases.

The script intentionally uses only the Python standard library so it can run on
fresh Windows/Linux machines before the rest of the data platform exists.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import datetime as dt
import hashlib
import html.parser
import json
import os
import re
import shutil
import sys
import time
import tomllib
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


USER_AGENT = "metabo-data-downloader/0.1 (+https://local.audit/raw-lake)"
LICENSED_LAYERS = {"licensed_enhancement", "controlled_access"}


@dataclass(frozen=True)
class Link:
    url: str
    text: str


@dataclass(frozen=True)
class PlannedFile:
    source_id: str
    url: str
    relative_path: str
    note: str = ""


class LinkParser(html.parser.HTMLParser):
    def __init__(self, base_url: str) -> None:
        super().__init__()
        self.base_url = base_url
        self.links: list[Link] = []
        self._href: str | None = None
        self._text_parts: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() != "a":
            return
        attrs_dict = dict(attrs)
        href = attrs_dict.get("href")
        if href:
            self._href = urllib.parse.urljoin(self.base_url, href)
            self._text_parts = []

    def handle_data(self, data: str) -> None:
        if self._href:
            self._text_parts.append(data)

    def handle_endtag(self, tag: str) -> None:
        if tag.lower() == "a" and self._href:
            text = " ".join("".join(self._text_parts).split())
            self.links.append(Link(url=self._href, text=text))
            self._href = None
            self._text_parts = []


def utc_release_id() -> str:
    return dt.datetime.now(dt.UTC).strftime("%Y%m%dT%H%M%SZ")


def read_catalog(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        data = tomllib.load(handle)
    if "sources" not in data:
        raise ValueError(f"Catalog {path} does not contain [[sources]] entries")
    return data


def source_selected(
    source: dict[str, Any],
    requested_sources: set[str] | None,
    requested_layers: set[str],
    include_disabled: bool,
    accept_licensed: bool,
) -> tuple[bool, str]:
    source_id = source["id"]
    layer = source.get("layer", "open_core")

    if requested_sources and source_id not in requested_sources:
        return False, "not requested"
    if "all" not in requested_layers and layer not in requested_layers:
        return False, f"layer {layer!r} not requested"
    if not include_disabled and not bool(source.get("enabled_by_default", True)):
        return False, "disabled by default"
    if layer in LICENSED_LAYERS and not accept_licensed:
        return False, f"licensed layer {layer!r} requires --accept-licensed"
    return True, "selected"


def sanitize_filename(value: str, fallback: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._() -]+", "_", value).strip(" ._")
    cleaned = re.sub(r"\s+", "_", cleaned)
    return cleaned or fallback


def basename_from_url(url: str) -> str:
    parsed = urllib.parse.urlparse(url)
    name = Path(urllib.parse.unquote(parsed.path)).name
    if name:
        return name
    query_name = urllib.parse.parse_qs(parsed.query).get("file", [""])[0]
    return sanitize_filename(query_name, "download.bin")


def request_url(url: str, timeout: int, headers: dict[str, str] | None = None) -> urllib.request.Request:
    merged_headers = {"User-Agent": USER_AGENT}
    if headers:
        merged_headers.update(headers)
    return urllib.request.Request(url, headers=merged_headers)


def fetch_text(url: str, timeout: int, retries: int) -> str:
    last_error: BaseException | None = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(request_url(url, timeout), timeout=timeout) as response:
                raw = response.read()
            return raw.decode("utf-8", errors="replace")
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last_error = exc
            if attempt == retries:
                break
            time.sleep(min(2 ** attempt, 30))
    raise RuntimeError(f"Failed to fetch index {url}: {last_error}") from last_error


def parse_index_links(url: str, timeout: int, retries: int) -> list[Link]:
    parser = LinkParser(url)
    parser.feed(fetch_text(url, timeout=timeout, retries=retries))
    return parser.links


def link_matches(link: Link, include_regex: str | None, exclude_regex: str | None) -> bool:
    parsed = urllib.parse.urlparse(link.url)
    name = Path(urllib.parse.unquote(parsed.path.rstrip("/"))).name
    candidate = f"{name} {link.text} {link.url}"
    if include_regex and not re.search(include_regex, candidate):
        return False
    if exclude_regex and re.search(exclude_regex, candidate):
        return False
    return True


def is_directory_link(link: Link) -> bool:
    parsed = urllib.parse.urlparse(link.url)
    name = Path(urllib.parse.unquote(parsed.path.rstrip("/"))).name
    if name in {"", ".", ".."} or link.text.lower().startswith("parent"):
        return False
    return parsed.path.endswith("/")


def relative_from_base(url: str, base_url: str, text: str = "", prefer_text: bool = False) -> str:
    if prefer_text and text:
        return sanitize_filename(text, basename_from_url(url))
    if url.startswith(base_url):
        rel = urllib.parse.unquote(url[len(base_url) :]).lstrip("/")
        if rel and not rel.startswith("?"):
            return rel
    fallback = basename_from_url(url)
    if fallback in {"download", "download.html"} and text:
        fallback = sanitize_filename(text, fallback)
    return fallback


def resolve_index_file_spec(
    source_id: str,
    spec: dict[str, Any],
    timeout: int,
    retries: int,
    resolve_indexes: bool,
) -> list[PlannedFile]:
    base_url = spec["url"]
    include_regex = spec.get("include_regex")
    exclude_regex = spec.get("exclude_regex")
    recursive = bool(spec.get("recursive", False))
    filename_from_text = bool(spec.get("filename_from_text", False))
    path_prefix = str(spec.get("path_prefix", "")).strip("/\\")
    max_depth = int(spec.get("max_depth", 2 if recursive else 0))

    if not resolve_indexes:
        return [PlannedFile(source_id, base_url, "__INDEX_NOT_RESOLVED__", "index; use --resolve-indexes")]

    planned: list[PlannedFile] = []
    seen_dirs: set[str] = set()
    seen_files: set[str] = set()

    def visit(index_url: str, depth: int) -> None:
        if index_url in seen_dirs:
            return
        seen_dirs.add(index_url)
        for link in parse_index_links(index_url, timeout=timeout, retries=retries):
            if is_directory_link(link):
                if recursive and depth < max_depth:
                    visit(link.url, depth + 1)
                continue
            if link_matches(link, include_regex, exclude_regex):
                if link.url in seen_files:
                    continue
                seen_files.add(link.url)
                rel = relative_from_base(link.url, base_url, text=link.text, prefer_text=filename_from_text)
                if path_prefix:
                    rel = f"{path_prefix}/{rel}"
                planned.append(PlannedFile(source_id, link.url, rel))

    try:
        visit(base_url, 0)
    except RuntimeError as exc:
        return [PlannedFile(source_id, base_url, "__INDEX_FAILED__", str(exc))]
    return planned


def resolve_api_cursor_spec(source_id: str, spec: dict[str, Any]) -> list[PlannedFile]:
    filename = spec.get("filename", "api-export.jsonl")
    return [PlannedFile(source_id, spec["url"], filename, "api_cursor")]


def resolve_source_files(
    source: dict[str, Any],
    timeout: int,
    retries: int,
    resolve_indexes: bool,
) -> list[PlannedFile]:
    source_id = source["id"]
    planned: list[PlannedFile] = []
    for spec in source.get("files", []):
        if not spec.get("enabled", True):
            continue
        kind = spec.get("kind", "url")
        if kind == "url":
            rel = spec.get("path") or basename_from_url(spec["url"])
            planned.append(PlannedFile(source_id, spec["url"], rel, spec.get("note", "")))
        elif kind == "index":
            planned.extend(resolve_index_file_spec(source_id, spec, timeout, retries, resolve_indexes))
        elif kind == "api_cursor":
            planned.extend(resolve_api_cursor_spec(source_id, spec))
        elif kind == "manual":
            continue
        else:
            raise ValueError(f"Unknown file spec kind {kind!r} in source {source_id}")
    return planned


def ensure_within(root: Path, target: Path) -> None:
    root_resolved = os.path.normcase(os.path.abspath(os.fspath(root)))
    target_resolved = os.path.normcase(os.path.abspath(os.fspath(target)))
    if os.path.commonpath([root_resolved, target_resolved]) != root_resolved:
        raise ValueError(f"Refusing to write outside {root_resolved}: {target_resolved}")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def md5_file(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_md5_sidecar(path: Path) -> str | None:
    text = path.read_text(encoding="utf-8", errors="replace")
    match = re.search(r"\b([0-9a-fA-F]{32})\b", text)
    return match.group(1).lower() if match else None


def download_url(
    planned: PlannedFile,
    destination: Path,
    timeout: int,
    retries: int,
    skip_existing: bool,
) -> dict[str, Any]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    part_path = destination.with_suffix(destination.suffix + ".part")

    if skip_existing and destination.exists():
        return {
            "url": planned.url,
            "path": str(destination),
            "status": "skipped_existing",
            "bytes": destination.stat().st_size,
            "sha256": sha256_file(destination),
        }

    last_error: BaseException | None = None
    for attempt in range(1, retries + 1):
        try:
            resume_from = part_path.stat().st_size if part_path.exists() else 0
            headers = {"Range": f"bytes={resume_from}-"} if resume_from else None
            request = request_url(planned.url, timeout=timeout, headers=headers)
            with urllib.request.urlopen(request, timeout=timeout) as response:
                status = getattr(response, "status", 200)
                mode = "ab" if resume_from and status == 206 else "wb"
                if resume_from and status != 206:
                    resume_from = 0
                with part_path.open(mode) as handle:
                    shutil.copyfileobj(response, handle, length=1024 * 1024)
            part_path.replace(destination)
            return {
                "url": planned.url,
                "path": str(destination),
                "status": "downloaded",
                "bytes": destination.stat().st_size,
                "sha256": sha256_file(destination),
                "resumed_from": resume_from,
            }
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last_error = exc
            if attempt == retries:
                break
            time.sleep(min(2 ** attempt, 30))

    return {
        "url": planned.url,
        "path": str(destination),
        "status": "failed",
        "error": str(last_error),
    }


def download_icite_cursor(
    planned: PlannedFile,
    destination: Path,
    spec: dict[str, Any],
    timeout: int,
    retries: int,
    max_pages: int | None,
) -> dict[str, Any]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    limit = int(spec.get("limit", 200))
    offset = int(spec.get("offset", 0))
    pages = 0
    records = 0
    last_error: BaseException | None = None

    with destination.open("w", encoding="utf-8") as out:
        while True:
            if max_pages is not None and pages >= max_pages:
                break
            params = {"limit": str(limit), "offset": str(offset), "legacy": "false"}
            url = planned.url + "?" + urllib.parse.urlencode(params)
            for attempt in range(1, retries + 1):
                try:
                    with urllib.request.urlopen(request_url(url, timeout), timeout=timeout) as response:
                        payload = json.loads(response.read().decode("utf-8"))
                    break
                except (urllib.error.URLError, TimeoutError, OSError, json.JSONDecodeError) as exc:
                    last_error = exc
                    if attempt == retries:
                        return {
                            "url": planned.url,
                            "path": str(destination),
                            "status": "failed",
                            "error": str(last_error),
                        }
                    time.sleep(min(2 ** attempt, 30))
            data = payload.get("data", payload if isinstance(payload, list) else [])
            if not data:
                break
            for row in data:
                out.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n")
            records += len(data)
            pages += 1
            max_pmid = max(int(row.get("pmid", offset)) for row in data if row.get("pmid") is not None)
            if max_pmid <= offset:
                break
            offset = max_pmid

    return {
        "url": planned.url,
        "path": str(destination),
        "status": "downloaded",
        "bytes": destination.stat().st_size,
        "sha256": sha256_file(destination),
        "pages": pages,
        "records": records,
        "final_offset": offset,
    }


def verify_md5_sidecars(entries: list[dict[str, Any]]) -> None:
    by_path = {Path(entry["path"]): entry for entry in entries if entry.get("path") and entry.get("status") in {"downloaded", "skipped_existing"}}
    for sidecar_path, sidecar_entry in list(by_path.items()):
        if sidecar_path.suffix.lower() != ".md5":
            continue
        expected = parse_md5_sidecar(sidecar_path)
        if not expected:
            sidecar_entry["md5_sidecar_status"] = "unreadable"
            continue
        target_name = sidecar_path.name[:-4]
        target_path = sidecar_path.with_name(target_name)
        target_entry = by_path.get(target_path)
        if not target_entry or not target_path.exists():
            sidecar_entry["md5_sidecar_status"] = "target_missing"
            continue
        actual = md5_file(target_path)
        target_entry["md5"] = actual
        target_entry["md5_expected"] = expected
        target_entry["md5_verified"] = actual == expected


def build_manifest(
    catalog_path: Path,
    release_id: str,
    output_root: Path,
    selected_sources: list[dict[str, Any]],
    entries: list[dict[str, Any]],
    dry_run: bool,
) -> dict[str, Any]:
    return {
        "release_id": release_id,
        "created_at_utc": dt.datetime.now(dt.UTC).isoformat(),
        "catalog": str(catalog_path),
        "output_root": str(output_root),
        "dry_run": dry_run,
        "sources": [
            {
                "id": source["id"],
                "name": source.get("name", source["id"]),
                "layer": source.get("layer", "open_core"),
                "official_url": source.get("official_url"),
                "license_policy": source.get("license_policy"),
            }
            for source in selected_sources
        ],
        "files": entries,
    }


def write_manifest(manifest: dict[str, Any], manifest_dir: Path, release_id: str) -> Path:
    manifest_dir.mkdir(parents=True, exist_ok=True)
    path = manifest_dir / f"{release_id}.manifest.json"
    path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")
    return path


def print_source_table(sources: Iterable[dict[str, Any]]) -> None:
    for source in sources:
        default = "yes" if source.get("enabled_by_default", True) else "no"
        print(f"{source['id']:<28} {source.get('layer', 'open_core'):<22} default={default:<3} {source.get('name', '')}")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download metabolism-oncology source databases into a versioned raw lake.")
    parser.add_argument("--catalog", default="config/source_catalog.toml", help="TOML source catalog.")
    parser.add_argument("--output-root", default="raw_lake", help="Root folder for downloaded files.")
    parser.add_argument("--manifest-dir", default="manifests", help="Folder for release manifests.")
    parser.add_argument("--release-id", default=utc_release_id(), help="Release/run id used in output paths.")
    parser.add_argument("--sources", help="Comma-separated source ids. Defaults to selected layers.")
    parser.add_argument("--layers", default="open_core", help="Comma-separated layers, or 'all'.")
    parser.add_argument("--include-disabled", action="store_true", help="Include sources marked enabled_by_default=false.")
    parser.add_argument("--accept-licensed", action="store_true", help="Allow licensed/controlled sources after manual review.")
    parser.add_argument("--list-sources", action="store_true", help="Print catalog source ids and exit.")
    parser.add_argument("--dry-run", action="store_true", help="Plan only; do not download.")
    parser.add_argument("--resolve-indexes", action="store_true", help="Fetch remote directory pages to enumerate dynamic files.")
    parser.add_argument("--yes", action="store_true", help="Required for real downloads.")
    parser.add_argument("--jobs", type=int, default=2, help="Parallel file downloads.")
    parser.add_argument("--timeout", type=int, default=90, help="Network timeout in seconds.")
    parser.add_argument("--retries", type=int, default=3, help="Network retry count.")
    parser.add_argument("--skip-existing", action="store_true", default=True, help="Skip files already present in the release folder.")
    parser.add_argument("--limit-per-source", type=int, help="Testing guard: cap planned files per source.")
    parser.add_argument("--max-api-pages", type=int, help="Testing guard for api_cursor sources.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    catalog_path = Path(args.catalog)
    output_root = Path(args.output_root)
    manifest_dir = Path(args.manifest_dir)
    catalog = read_catalog(catalog_path)
    all_sources = catalog["sources"]

    if args.list_sources:
        print_source_table(all_sources)
        return 0

    requested_sources = set(filter(None, (args.sources or "").split(","))) or None
    requested_layers = set(filter(None, args.layers.split(","))) or {"open_core"}

    selected_sources: list[dict[str, Any]] = []
    skipped: list[tuple[str, str]] = []
    for source in all_sources:
        selected, reason = source_selected(
            source,
            requested_sources=requested_sources,
            requested_layers=requested_layers,
            include_disabled=args.include_disabled,
            accept_licensed=args.accept_licensed,
        )
        if selected:
            selected_sources.append(source)
        else:
            skipped.append((source["id"], reason))

    if not selected_sources:
        print("No sources selected. Use --list-sources to inspect the catalog.", file=sys.stderr)
        return 2

    if not args.dry_run and not args.yes:
        print("Refusing to download without --yes. Run --dry-run first, then add --yes.", file=sys.stderr)
        return 2

    output_root.mkdir(parents=True, exist_ok=True)
    ensure_within(Path.cwd(), output_root)
    ensure_within(Path.cwd(), manifest_dir)

    planned_by_source: dict[str, list[PlannedFile]] = {}
    api_specs: dict[tuple[str, str], dict[str, Any]] = {}
    for source in selected_sources:
        planned = resolve_source_files(
            source,
            timeout=args.timeout,
            retries=args.retries,
            resolve_indexes=args.resolve_indexes,
        )
        if args.limit_per_source is not None:
            planned = planned[: args.limit_per_source]
        planned_by_source[source["id"]] = planned
        for spec in source.get("files", []):
            if spec.get("kind") == "api_cursor":
                api_specs[(source["id"], spec.get("filename", "api-export.jsonl"))] = spec

    entries: list[dict[str, Any]] = []
    if args.dry_run:
        for source in selected_sources:
            for planned in planned_by_source[source["id"]]:
                entries.append(
                    {
                        "source_id": planned.source_id,
                        "url": planned.url,
                        "relative_path": planned.relative_path,
                        "note": planned.note,
                        "status": "planned",
                    }
                )
        manifest = build_manifest(catalog_path, args.release_id, output_root, selected_sources, entries, dry_run=True)
        manifest_path = write_manifest(manifest, manifest_dir, args.release_id)
        print(f"Dry run planned {len(entries)} files across {len(selected_sources)} sources.")
        print(f"Manifest: {manifest_path}")
        unresolved = [entry for entry in entries if entry["relative_path"] == "__INDEX_NOT_RESOLVED__"]
        if unresolved:
            print(f"{len(unresolved)} index specs were not resolved. Add --resolve-indexes to enumerate remote files.")
        return 0

    futures: list[concurrent.futures.Future[dict[str, Any]]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        for source in selected_sources:
            source_root = output_root / source["id"] / args.release_id
            for planned in planned_by_source[source["id"]]:
                if planned.relative_path in {"__INDEX_NOT_RESOLVED__", "__INDEX_FAILED__"}:
                    entries.append(
                        {
                            "source_id": source["id"],
                            "url": planned.url,
                            "status": "skipped_index",
                            "note": planned.note or "run with --resolve-indexes",
                        }
                    )
                    continue
                destination = source_root / planned.relative_path
                ensure_within(output_root, destination)
                spec = api_specs.get((source["id"], planned.relative_path))
                if spec:
                    futures.append(
                        pool.submit(
                            download_icite_cursor,
                            planned,
                            destination,
                            spec,
                            args.timeout,
                            args.retries,
                            args.max_api_pages,
                        )
                    )
                else:
                    futures.append(
                        pool.submit(
                            download_url,
                            planned,
                            destination,
                            args.timeout,
                            args.retries,
                            args.skip_existing,
                        )
                    )

        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            result["source_id"] = next((planned.source_id for plans in planned_by_source.values() for planned in plans if planned.url == result.get("url")), None)
            entries.append(result)
            status = result.get("status")
            path = result.get("path", result.get("url"))
            print(f"[{status}] {path}")

    verify_md5_sidecars(entries)
    manifest = build_manifest(catalog_path, args.release_id, output_root, selected_sources, entries, dry_run=False)
    manifest_path = write_manifest(manifest, manifest_dir, args.release_id)
    failed = sum(1 for entry in entries if entry.get("status") == "failed")
    print(f"Finished {len(entries)} file jobs with {failed} failures.")
    print(f"Manifest: {manifest_path}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
