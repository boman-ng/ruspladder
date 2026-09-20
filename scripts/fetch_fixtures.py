#!/usr/bin/env python3
"""Fetch and verify the committed airway fixture without replacing existing data."""
import argparse
import json
import os
import shutil
import tarfile
import tempfile
import urllib.parse
import urllib.request
from pathlib import Path, PurePosixPath

from fixture_manifest import sha256, verify_files

MANIFEST = Path(__file__).resolve().parent.parent / "fixtures/airway-manifest.json"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-root", type=Path, default=Path(
        os.environ.get("RUSPLADDER_WORK_ROOT", Path.home() / "data/ruspladder")))
    return parser.parse_args(argv)


def fetch(work_root, manifest):
    fixtures = work_root / "fixtures"
    fixtures.mkdir(parents=True, exist_ok=True)
    archive = fixtures / PurePosixPath(urllib.parse.urlparse(manifest["source"]).path).name
    destination = fixtures / "airway"
    with tempfile.TemporaryDirectory(prefix=".airway-", dir=fixtures) as temporary:
        temporary = Path(temporary)
        download = archive if archive.exists() else temporary / "download"
        if not archive.exists():
            with urllib.request.urlopen(manifest["source"]) as source, download.open("wb") as output:
                shutil.copyfileobj(source, output)
        if sha256(download) != manifest["archive_sha256"]:
            raise ValueError(f"Archive checksum mismatch: {archive.name}")
        staged = temporary / "airway"
        staged.mkdir()
        expected = {item["path"] for item in manifest["files"]}
        if len(expected) != len(manifest["files"]) or any(
            len(PurePosixPath(name).parts) != 1 or name in (".", "..") for name in expected
        ):
            raise ValueError("Expected unique fixture basenames")
        found = set()
        with tarfile.open(download, "r:gz") as source:
            for member in source.getmembers():
                parts = PurePosixPath(member.name).parts
                if len(parts) != 4 or parts[:3] != ("airway", "inst", "extdata"):
                    continue
                name = parts[-1]
                if name not in expected:
                    continue
                if name in found or not member.isfile():
                    raise ValueError(f"Invalid or duplicate archive fixture: {name}")
                with source.extractfile(member) as content, (staged / name).open("wb") as output:
                    shutil.copyfileobj(content, output)
                found.add(name)
        if found != expected:
            raise ValueError(f"Missing archive fixtures: {sorted(expected - found)}")
        verify_files(staged, manifest)
        (staged / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        if destination.exists():
            verify_files(destination, manifest)
        else:
            staged.rename(destination)
        if not archive.exists():
            download.rename(archive)
    print(f"Verified {len(expected)} fixtures in {destination}")
    return destination


def main():
    args = parse_args()
    try:
        fetch(args.work_root, json.loads(MANIFEST.read_text()))
    except (OSError, ValueError, tarfile.TarError) as error:
        raise SystemExit(str(error)) from error


if __name__ == "__main__":
    main()
