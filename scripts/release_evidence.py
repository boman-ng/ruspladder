#!/usr/bin/env python3
"""Bind successful public comparisons to the exact runtime and source attachments."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile
import tomllib

from fixture_manifest import sha256
from native_sources import validate_provenance

CASES = {"basic-pos-merge", "basic-neg-merge", "basic-pos-single",
         "basic-neg-single", "events-bam", "events-cram"}


def command(repo, *args):
    return subprocess.check_output(args, cwd=repo, text=True).strip()


def write_checksum(path):
    path.with_name(path.name + ".sha256").write_text(f"{sha256(path)}  {path.name}\n")


def build_provenance(repo, binary, lib):
    names = command(repo, "git", "ls-files", "--cached", "--others", "--exclude-standard", "-z")
    inputs = {name: sha256(repo / name) for name in sorted(set(names.split("\0"))) if name}
    return {
        "source_commit": command(repo, "git", "rev-parse", "HEAD"),
        "source_dirty": bool(command(repo, "git", "status", "--porcelain")),
        "source_files_sha256": inputs,
        "source_tree_sha256": hashlib.sha256(json.dumps(inputs, sort_keys=True).encode()).hexdigest(),
        "binary_sha256": sha256(binary),
        "native_libraries_sha256": {p.name: sha256(p) for p in sorted(lib.iterdir())},
        "toolchain": {tool: command(repo, tool, "--version").splitlines()[0]
                      for tool in ["rustc", "cargo", "cc", "ld", "patchelf", "python3"]},
        "bundled_system_packages": {
            name: command(repo, "rpm", "-q", "--qf", "%{NAME}-%{VERSION}-%{RELEASE}.%{ARCH}\n%{SOURCERPM}", name).splitlines()
            for name in ["zlib", "bzip2-libs"]
        },
    }


def check_report(report, binary_hash):
    cases = report.get("cases", [])
    if (report.get("result") != "pass" or len(cases) != len(CASES)
            or {case.get("case") for case in cases} != CASES
            or any(case.get("result") != "pass" for case in cases)):
        raise ValueError("All six successful public comparisons are required")
    if report.get("binary_sha256") != binary_hash:
        raise ValueError("Public comparisons used a different runtime binary")


def read_member(archive, name):
    member = archive.getmember(name)
    if not member.isfile():
        raise ValueError(f"Expected archive file: {name}")
    with archive.extractfile(member) as source:
        return source.read()


def verify(evidence_path):
    evidence = json.loads(evidence_path.read_text())
    if evidence.get("result") != "pass":
        raise ValueError("Missing successful release verification")
    check_report(evidence["public_tests"], evidence["build"]["binary_sha256"])
    for name, expected in evidence["artifacts"].items():
        if Path(name).name != name:
            raise ValueError("Expected artifact basename")
        path = evidence_path.parent / name
        if sha256(path) != expected:
            raise ValueError(f"Release artifact checksum mismatch: {name}")
        if path.with_name(name + ".sha256").read_text() != f"{expected}  {name}\n":
            raise ValueError(f"Release checksum file mismatch: {name}")
    runtime = evidence_path.parent / evidence["runtime_archive"]
    base = runtime.name.removesuffix(".tar.gz")
    with tarfile.open(runtime, "r:gz") as archive:
        build = json.loads(read_member(archive, f"{base}/build-provenance.json"))
        if build != evidence["build"]:
            raise ValueError("Runtime provenance differs from verification attachment")
        if hashlib.sha256(read_member(archive, f"{base}/ruspladder")).hexdigest() != build["binary_sha256"]:
            raise ValueError("Runtime binary differs from tested binary")
        for name, expected in build["native_libraries_sha256"].items():
            if hashlib.sha256(read_member(archive, f"{base}/lib/{name}")).hexdigest() != expected:
                raise ValueError(f"Runtime library differs from provenance: {name}")
        for name in ["LICENSE", "licenses/LGPL-2.1.txt", "licenses/NumPy-wheel-native-notices.txt",
                     "licenses/README.md", "licenses/dependencies.json"]:
            if not read_member(archive, f"{base}/{name}"):
                raise ValueError(f"Empty redistribution material: {name}")
        native = json.loads(read_member(archive, f"{base}/licenses/native-runtime-provenance.json"))
        validate_provenance(native)
    sources = evidence_path.parent / evidence["source_archive"]
    base = sources.name.removesuffix(".tar.gz")
    with tarfile.open(sources, "r:gz") as archive:
        if json.loads(read_member(archive, f"{base}/native-runtime-provenance.json")) != native:
            raise ValueError("Source attachment provenance mismatch")
        for item in native["runtime_libraries"]:
            source = item["source"]
            with archive.extractfile(f"{base}/{source['filename']}") as data:
                if hashlib.file_digest(data, "sha256").hexdigest() != source["sha256"]:
                    raise ValueError(f"Corresponding source mismatch: {source['filename']}")
    return evidence


def finalize(repo, out, report_path):
    version = tomllib.loads((repo / "Cargo.toml").read_text())["package"]["version"]
    prefix = f"ruspladder-v{version}"
    runtime = out / f"{prefix}-linux-amd64.tar.gz"
    sources = out / f"{prefix}-native-sources.tar.gz"
    with tarfile.open(runtime, "r:gz") as archive:
        build = json.loads(read_member(archive, f"{prefix}-linux-amd64/build-provenance.json"))
    report = json.loads(report_path.read_text())
    check_report(report, build["binary_sha256"])
    upstream = json.loads((repo / "fixtures/upstream-manifest.json").read_text())
    if report["reference_commit"] != upstream["commit"]:
        raise ValueError("Public comparisons used a different upstream commit")
    for case in report["cases"]:
        case.pop("commands", None)  # Absolute local paths do not belong in public evidence.
    evidence = {
        "result": "pass", "runtime_archive": runtime.name, "source_archive": sources.name,
        "build": build, "public_tests": report,
        "input_manifests_sha256": {name: sha256(repo / "fixtures" / name)
                                   for name in ["airway-manifest.json", "upstream-manifest.json"]},
        "artifacts": {p.name: sha256(p) for p in [runtime, sources]},
        "checks": ["package members and hashes", "native source hashes", "six public comparisons"],
    }
    path = out / f"{prefix}-verification.json"
    path.write_text(json.dumps(evidence, indent=2) + "\n")
    verify(path)
    write_checksum(path)
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--verify", type=Path)
    parser.add_argument("--require-commit")
    args = parser.parse_args()
    if args.verify:
        evidence = verify(args.verify)
        if args.require_commit and (evidence["build"]["source_dirty"] or
                                    evidence["build"]["source_commit"] != args.require_commit):
            raise SystemExit("Release evidence must come from the clean tagged commit")
        print("Verified release evidence:", args.verify)
    elif args.out and args.report:
        print(finalize(Path(__file__).resolve().parent.parent, args.out, args.report))
    else:
        parser.error("use --verify, or both --out and --report")


if __name__ == "__main__":
    main()
