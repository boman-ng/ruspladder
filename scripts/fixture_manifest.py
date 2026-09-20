"""Verify the committed public-data manifests before scientific comparisons."""
import hashlib
from pathlib import Path, PurePosixPath


def sha256(path):
    with Path(path).open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def verify_files(root, manifest):
    for item in manifest["files"]:
        relative = PurePosixPath(item["path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"Unsafe fixture path: {relative}")
        path = Path(root) / relative
        if not path.is_file():
            raise ValueError(f"Missing fixture: {relative}")
        if path.stat().st_size != item["bytes"] or sha256(path) != item["sha256"]:
            raise ValueError(f"Fixture checksum mismatch: {relative}")
