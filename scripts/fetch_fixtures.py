#!/usr/bin/env python3
"""Fetch only the pinned small airway archive; retain relevant extdata files."""
import argparse
import hashlib
import json
import tarfile
import urllib.request
from pathlib import Path, PurePosixPath

URL = "https://bioconductor.org/packages/3.23/data/experiment/src/contrib/airway_1.32.0.tar.gz"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-root", type=Path, default=Path("/home/wubw/data/ruspladder"))
    args = parser.parse_args()
    archive = args.work_root / "fixtures/airway_1.32.0.tar.gz"
    destination = args.work_root / "fixtures/airway"
    destination.mkdir(parents=True, exist_ok=True)
    if not archive.exists():
        temporary = archive.with_suffix(".download")
        with urllib.request.urlopen(URL) as response, temporary.open("wb") as output:
            while chunk := response.read(1024 * 1024):
                output.write(chunk)
        temporary.rename(archive)
    files = []
    with tarfile.open(archive, "r:gz") as source:
        for member in source.getmembers():
            parts = PurePosixPath(member.name).parts
            if len(parts) != 4 or parts[:3] != ("airway", "inst", "extdata") or not member.isfile():
                continue
            name = parts[-1]
            if not (name.endswith((".bam", ".gtf")) or name in
                    ("sample_table.csv", "SraRunInfo_SRP033351.csv", "GSE52778_series_matrix.txt")):
                continue
            data = source.extractfile(member).read()
            (destination / name).write_bytes(data)
            files.append({"path": name, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()})
    manifest = {"source": URL, "archive_sha256": hashlib.sha256(archive.read_bytes()).hexdigest(), "files": files}
    (destination / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Fetched {len(files)} files ({sum(x['bytes'] for x in files)} bytes) into {destination}")


if __name__ == "__main__":
    main()
