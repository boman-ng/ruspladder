#!/usr/bin/env python3
"""Extract the pinned NumPy reference's native BLAS kernels; no Python runtime dependency.

The scientific parity baseline currently targets Linux x86_64. The wheel is a
verified source of OpenBLAS 0.3.29 ILP64 and its redistributable native runtimes;
none of NumPy's Python modules or extension modules are installed or loaded.
"""
import hashlib
import os
from pathlib import Path
import platform
import urllib.request
import zipfile

ROOT = Path(os.environ.get("RUSPLADDER_WORK_ROOT", Path.home() / "data/ruspladder"))
URL = "https://files.pythonhosted.org/packages/8c/3d/1e1db36cfd41f895d266b103df00ca5b3cbe965184df824dec5c08c6b803/numpy-2.2.6-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl"
SHA256 = "fd83c01228a688733f1ded5201c678f0c53ecc1006ffbc404db9f7a899ac6249"
LIBRARIES = {
    "libscipy_openblas64_-56d6093b.so": "0bd815d04b6b54990e3cccc7528fbb696456d09569f533d0390c13f0cdc4dd4a",
    "libgfortran-040039e1-0352e75f.so.5.0.0": "c6090048eccc763522c12ef016f81da6b627cb3a044f55cf0479a839c41c0980",
    "libquadmath-96973f99-934c22de.so.0.0.0": "6ed5137f412781ad7863439fb543613f620b43c32b63292a0029246162f5bbc6",
}

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
