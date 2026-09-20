"""Fetch corresponding GCC runtime sources for release redistribution."""
import json
from pathlib import Path
import shutil
import tarfile
import tempfile
import urllib.request

from fixture_manifest import sha256


def validate_provenance(manifest):
    required = {name for name in manifest["libraries"]
                if name.startswith(("libgfortran", "libquadmath"))}
    records = manifest["runtime_libraries"]
    if {item["wheel_library"] for item in records} != required or len(records) != len(required):
        raise ValueError("Corresponding-source records must cover each bundled GCC runtime")
    for item in records:
        if manifest["libraries"][item["wheel_library"]] != item["wheel_library_sha256"]:
            raise ValueError(f"Native source provenance does not match pinned library: {item['wheel_library']}")


def package_sources(repo, work, out, version):
    provenance = repo / "licenses/native-runtime-provenance.json"
    manifest = json.loads(provenance.read_text())
    validate_provenance(manifest)
    libraries = manifest["runtime_libraries"]
    cache = work / "cache/native-sources"
    cache.mkdir(parents=True, exist_ok=True)
    sources = []
    for library in libraries:
        item = library["source"]
        source = cache / item["filename"]
        if not source.exists():
            with tempfile.TemporaryDirectory(dir=cache) as temporary:
                download = Path(temporary) / item["filename"]
                with urllib.request.urlopen(item["url"]) as response, download.open("wb") as output:
                    shutil.copyfileobj(response, output)
                if sha256(download) != item["sha256"]:
                    raise ValueError(f"Native source checksum mismatch: {item['filename']}")
                download.rename(source)
        if sha256(source) != item["sha256"]:
            raise ValueError(f"Cached native source checksum mismatch: {source.name}")
        sources.append(source)
    basename = f"ruspladder-v{version}-native-sources"
    archive = out / f"{basename}.tar.gz"
    with tarfile.open(archive, "w:gz") as output:
        for path in sources:
            output.add(path, arcname=f"{basename}/{path.name}")
        for name in ["README.md", "native-runtime-provenance.json", "LGPL-2.1.txt",
                     "NumPy-wheel-native-notices.txt"]:
            output.add(repo / "licenses" / name, arcname=f"{basename}/{name}")
    return archive
