#!/usr/bin/env python3
"""Package the Linux amd64 executable, native libraries and dependency notices."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tarfile
import tomllib


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parent.parent
    version = tomllib.loads((repo / "Cargo.toml").read_text())["package"]["version"]
    work = Path(os.environ.get("RUSPLADDER_WORK_ROOT", Path.home() / "data/ruspladder"))
    target = Path(os.environ.get("CARGO_TARGET_DIR", work / "target"))
    package = args.out.resolve() / f"ruspladder-v{version}-linux-amd64"
    package.mkdir(parents=True)
    lib = package / "lib"
    lib.mkdir()
    binary = package / "ruspladder"
    shutil.copy2(target / "release/ruspladder", binary)
    subprocess.run(["strip", str(binary)], check=True)
    for path in (work / "native/openblas/lib").glob("*.so*"):
        shutil.copy2(path, lib / path.name)
    # These two libraries are used by HTSlib/HDF5. Keep glibc and the
    # compiler's system runtime supplied by the host operating system.
    needed = subprocess.check_output(["ldd", str(binary)], text=True)
    for line in needed.splitlines():
        fields = line.split()
        if fields and fields[0].startswith(("libbz2.so", "libz.so")):
            shutil.copy2(fields[2], lib / fields[0])
    subprocess.run(["patchelf", "--set-rpath", "$ORIGIN/lib", str(binary)], check=True)
    for path in lib.glob("*.so*"):
        subprocess.run(["patchelf", "--set-rpath", "$ORIGIN", str(path)], check=True)
    for path in [binary, *lib.glob("*.so*")]:
        symbols = subprocess.check_output(["readelf", "--version-info", str(path)], text=True)
        required = [tuple(map(int, v)) for v in re.findall(r"Name: GLIBC_(\d+)\.(\d+)", symbols)]
        if required and max(required) > (2, 28):
            raise RuntimeError(f"{path.name} exceeds the glibc 2.28 release baseline: {max(required)}")
    for name in ["LICENSE", "README.md", "COMPATIBILITY.md", "NUMERICS.md", "CITATION.cff"]:
        shutil.copy2(repo / name, package / name)
    notices = package / "licenses"
    shutil.copytree(repo / "licenses", notices)
    shutil.copy2(repo / "vendor/x86-simd-sort/LICENSE.md", notices / "x86-simd-sort-BSD.txt")
    sysroot = Path(subprocess.check_output(["rustc", "--print", "sysroot"], text=True).strip())
    rust_notices = sysroot / "share/doc/rust"
    shutil.copy2(rust_notices / "COPYRIGHT-library.html", notices / "Rust-standard-library.html")
    shutil.copytree(rust_notices / "licenses", notices / "rust")
    for name in ["bzip2-libs", "zlib"]:
        shutil.copytree(Path("/usr/share/licenses") / name, notices / name)
    metadata = json.loads(subprocess.check_output(
        [str(repo / "scripts/cargo.sh"), "metadata", "--locked", "--format-version", "1"],
        cwd=repo, text=True))
    dependencies = []
    for dependency in sorted(metadata["packages"], key=lambda p: (p["name"], p["version"])):
        if not dependency["source"]:
            continue
        name = f'{dependency["name"]}-{dependency["version"]}'
        source = Path(dependency["manifest_path"]).parent
        # Include nested native-library notices as well as crate licenses.
        for path in source.rglob("*"):
            if path.is_file() and path.name.upper().startswith(("LICENSE", "LICENCE", "COPYING", "COPYRIGHT", "NOTICE")):
                dest = notices / "crates" / name / path.relative_to(source)
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, dest)
        dependencies.append({k: dependency[k] for k in ["name", "version", "license", "repository"]})
    (notices / "dependencies.json").write_text(json.dumps(dependencies, indent=2) + "\n")
    archive = package.with_name(package.name + ".tar.gz")
    with tarfile.open(archive, "w:gz") as output:
        output.add(package, arcname=package.name)
    with archive.open("rb") as source:
        digest = hashlib.file_digest(source, "sha256").hexdigest()
    archive.with_name(archive.name + ".sha256").write_text(f"{digest}  {archive.name}\n")
    print(archive)


if __name__ == "__main__":
    main()
