#!/usr/bin/env python3
"""Extract the pinned NumPy reference's native BLAS kernels; no Python runtime dependency.

The scientific parity baseline currently targets Linux x86_64. The wheel is a
verified source of OpenBLAS 0.3.29 ILP64 and its redistributable native runtimes;
none of NumPy's Python modules or extension modules are installed or loaded.
"""
import hashlib
import json
import os
from pathlib import Path
import platform
import urllib.request
import zipfile

ROOT = Path(os.environ.get("RUSPLADDER_WORK_ROOT", Path.home() / "data/ruspladder"))
PROVENANCE = json.loads((Path(__file__).resolve().parent.parent /
                         "licenses/native-runtime-provenance.json").read_text())
URL = PROVENANCE["wheel"]["url"]
SHA256 = PROVENANCE["wheel"]["sha256"]
LIBRARIES = PROVENANCE["libraries"]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    if (platform.system(), platform.machine()) != ("Linux", "x86_64"):
        raise SystemExit("Pinned native parity baseline requires Linux x86_64")
    destination = ROOT / "native/openblas/lib"
    destination.mkdir(parents=True, exist_ok=True)
    if all((destination / name).exists() and digest(destination / name) == sha for name, sha in LIBRARIES.items()):
        print("Pinned native libraries verified:", destination)
        return
    archive = ROOT / "cache/native" / URL.rsplit("/", 1)[1]
    archive.parent.mkdir(parents=True, exist_ok=True)
    if not archive.exists():
        temporary = archive.with_suffix(".download")
        urllib.request.urlretrieve(URL, temporary)
        if digest(temporary) != SHA256:
            raise SystemExit("Native dependency archive checksum mismatch")
        temporary.replace(archive)
    if digest(archive) != SHA256:
        raise SystemExit("Cached native dependency archive checksum mismatch")
    with zipfile.ZipFile(archive) as wheel:
        for name, sha in LIBRARIES.items():
            data = wheel.read("numpy.libs/" + name)
            if hashlib.sha256(data).hexdigest() != sha:
                raise SystemExit("Native library checksum mismatch: " + name)
            (destination / name).write_bytes(data)
        (destination.parent / "LICENSE.txt").write_bytes(wheel.read("numpy-2.2.6.dist-info/LICENSE.txt"))
    print("Installed native libraries:", destination)

if __name__ == "__main__":
    main()
