"""Small offline tests for fixture integrity and release evidence."""
import copy
import hashlib
import io
import json
import os
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch

from fetch_fixtures import fetch, parse_args
from fixture_manifest import sha256, verify_files
from native_sources import package_sources, validate_provenance
from release_evidence import CASES, check_report, verify, write_checksum


class FixtureTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.archive = self.root / "airway.tar.gz"
        self.work = self.root / "work"
        self.files = {"sample.bam": b"alignment", "annotation.gtf": b"annotation"}
        with tarfile.open(self.archive, "w:gz") as output:
            for name, content in self.files.items():
                info = tarfile.TarInfo("airway/inst/extdata/" + name)
                info.size = len(content)
                output.addfile(info, io.BytesIO(content))
        self.manifest = {
            "source": self.archive.as_uri(), "archive_sha256": sha256(self.archive),
            "files": [{"path": name, "bytes": len(content),
                       "sha256": hashlib.sha256(content).hexdigest()}
                      for name, content in self.files.items()],
        }

    def test_directory_precedence(self):
        with patch.dict(os.environ, {"RUSPLADDER_WORK_ROOT": str(self.work)}):
            self.assertEqual(parse_args([]).work_root, self.work)
            self.assertEqual(parse_args(["--work-root", "/explicit"]).work_root, Path("/explicit"))
        with patch.dict(os.environ, {}, clear=True), patch.object(Path, "home", return_value=self.root):
            self.assertEqual(parse_args([]).work_root, self.root / "data/ruspladder")

    def test_verified_download_and_offline_cache_reuse(self):
        destination = fetch(self.work, self.manifest)
        verify_files(destination, self.manifest)
        before = (destination / "sample.bam").stat().st_mtime_ns
        self.archive.unlink()
        fetch(self.work, self.manifest)
        self.assertEqual((destination / "sample.bam").stat().st_mtime_ns, before)
        self.assertEqual(json.loads((destination / "manifest.json").read_text()), self.manifest)

    def test_corrupt_download_is_not_installed(self):
        self.archive.write_bytes(b"corrupt")
        with self.assertRaisesRegex(ValueError, "Archive checksum"):
            fetch(self.work, self.manifest)
        self.assertFalse((self.work / "fixtures/airway").exists())
        self.assertFalse((self.work / "fixtures/airway.tar.gz").exists())
        self.assertEqual(list((self.work / "fixtures").iterdir()), [])

    def test_corrupt_cache_preserves_destination(self):
        destination = fetch(self.work, self.manifest)
        (self.work / "fixtures/airway.tar.gz").write_bytes(b"corrupt")
        with self.assertRaisesRegex(ValueError, "Archive checksum"):
            fetch(self.work, self.manifest)
        self.assertEqual((destination / "sample.bam").read_bytes(), self.files["sample.bam"])

    def test_missing_archive_member_is_not_installed(self):
        self.manifest["files"].append({"path": "missing.bam", "bytes": 0, "sha256": "0" * 64})
        with self.assertRaisesRegex(ValueError, "Missing archive fixtures"):
            fetch(self.work, self.manifest)
        self.assertFalse((self.work / "fixtures/airway").exists())

    def test_wrong_file_hash_is_not_installed(self):
        self.manifest["files"][0]["sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "Fixture checksum"):
            fetch(self.work, self.manifest)
        self.assertFalse((self.work / "fixtures/airway").exists())

    def test_existing_modified_fixture_is_not_overwritten(self):
        destination = fetch(self.work, self.manifest)
        (destination / "sample.bam").write_bytes(b"user data")
        with self.assertRaisesRegex(ValueError, "Fixture checksum"):
            fetch(self.work, self.manifest)
        self.assertEqual((destination / "sample.bam").read_bytes(), b"user data")

    def test_missing_upstream_file_and_changed_same_size_file(self):
        destination = fetch(self.work, self.manifest)
        (destination / "sample.bam").unlink()
        with self.assertRaisesRegex(ValueError, "Missing fixture"):
            verify_files(destination, self.manifest)
        (destination / "sample.bam").write_bytes(b"ALIGNMENT")
        with self.assertRaisesRegex(ValueError, "Fixture checksum"):
            verify_files(destination, self.manifest)

    def test_manifest_cannot_escape_destination(self):
        manifest = copy.deepcopy(self.manifest)
        manifest["files"][0]["path"] = "../outside"
        with self.assertRaisesRegex(ValueError, "basenames"):
            fetch(self.work, manifest)
        with self.assertRaisesRegex(ValueError, "Unsafe fixture"):
            verify_files(self.root, manifest)


class NativeSourceTests(unittest.TestCase):
    def test_sources_download_and_corrupt_cache_rejection(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            repo, work, out = root / "repo", root / "work", root / "out"
            (repo / "licenses").mkdir(parents=True)
            out.mkdir()
            source = root / "gcc.src.rpm"
            source.write_bytes(b"test source package")
            provenance = {"libraries": {"libquadmath-test.so": "test-hash"},
                          "runtime_libraries": [{"wheel_library": "libquadmath-test.so",
                                                 "wheel_library_sha256": "test-hash", "source": {
                "filename": source.name, "url": source.as_uri(), "sha256": sha256(source)}}]}
            (repo / "licenses/native-runtime-provenance.json").write_text(json.dumps(provenance))
            for name in ["README.md", "LGPL-2.1.txt", "NumPy-wheel-native-notices.txt"]:
                (repo / "licenses" / name).write_text("test notice")
            archive = package_sources(repo, work, out, "test")
            with tarfile.open(archive) as package:
                self.assertEqual(package.extractfile("ruspladder-vtest-native-sources/gcc.src.rpm").read(),
                                 source.read_bytes())
            source.unlink()
            package_sources(repo, work, out, "test")  # No network needed for valid cache.
            (work / "cache/native-sources/gcc.src.rpm").write_bytes(b"wrong source")
            with self.assertRaisesRegex(ValueError, "Cached native source checksum"):
                package_sources(repo, work, out, "test")

    def test_changed_runtime_requires_new_source_correspondence(self):
        manifest = {"libraries": {"libquadmath-test.so": "new"}, "runtime_libraries": []}
        with self.assertRaisesRegex(ValueError, "cover each"):
            validate_provenance(manifest)
        manifest["runtime_libraries"] = [{"wheel_library": "libquadmath-test.so",
                                          "wheel_library_sha256": "old"}]
        with self.assertRaisesRegex(ValueError, "does not match"):
            validate_provenance(manifest)


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.runtime = self.root / "ruspladder-vtest-linux-amd64.tar.gz"
        self.sources = self.root / "ruspladder-vtest-native-sources.tar.gz"
        self.path = self.root / "verification.json"
        digest = lambda b: hashlib.sha256(b).hexdigest()
        self.build = {"source_commit": "test", "source_dirty": False,
                      "binary_sha256": digest(b"binary"),
                      "native_libraries_sha256": {"libtest.so": digest(b"library")}}
        self.native = {"libraries": {"libquadmath-test.so": "test-hash"},
                       "runtime_libraries": [{"wheel_library": "libquadmath-test.so",
                                              "wheel_library_sha256": "test-hash", "source": {
            "filename": "gcc.src.rpm", "sha256": digest(b"source")}}]}
        self.runtime_files = {"ruspladder": b"binary", "lib/libtest.so": b"library",
                              "build-provenance.json": json.dumps(self.build).encode(),
                              "licenses/native-runtime-provenance.json": json.dumps(self.native).encode()}
        for name in ["LICENSE", "licenses/LGPL-2.1.txt", "licenses/NumPy-wheel-native-notices.txt",
                     "licenses/README.md", "licenses/dependencies.json"]:
            self.runtime_files[name] = b"test notice"
        self.source_files = {"gcc.src.rpm": b"source",
                             "native-runtime-provenance.json": json.dumps(self.native).encode()}
        self.report = {"result": "pass", "binary_sha256": self.build["binary_sha256"],
                       "cases": [{"case": name, "result": "pass"} for name in sorted(CASES)]}
        self.write_artifacts()

    def write_artifacts(self):
        for path, files in [(self.runtime, self.runtime_files), (self.sources, self.source_files)]:
            with tarfile.open(path, "w:gz") as archive:
                for name, value in files.items():
                    entry = tarfile.TarInfo(path.name.removesuffix(".tar.gz") + "/" + name)
                    entry.size = len(value)
                    archive.addfile(entry, io.BytesIO(value))
            write_checksum(path)
        self.evidence = {"result": "pass", "runtime_archive": self.runtime.name,
                         "source_archive": self.sources.name, "build": self.build,
                         "public_tests": self.report,
                         "artifacts": {p.name: sha256(p) for p in [self.runtime, self.sources]}}
        self.path.write_text(json.dumps(self.evidence))

    def test_complete_evidence_and_modified_archive(self):
        verify(self.path)
        with self.runtime.open("ab") as file:
            file.write(b"altered")
        with self.assertRaisesRegex(ValueError, "artifact checksum"):
            verify(self.path)

    def test_missing_case_and_wrong_binary_are_rejected(self):
        report = copy.deepcopy(self.report)
        report["cases"].pop()
        with self.assertRaisesRegex(ValueError, "six successful"):
            check_report(report, self.build["binary_sha256"])
        with self.assertRaisesRegex(ValueError, "different runtime"):
            check_report(self.report, "different")

    def test_wrong_library_inside_correctly_hashed_archive(self):
        self.runtime_files["lib/libtest.so"] = b"different library"
        self.write_artifacts()
        with self.assertRaisesRegex(ValueError, "library differs"):
            verify(self.path)

    def test_wrong_source_inside_correctly_hashed_archive(self):
        self.source_files["gcc.src.rpm"] = b"different source"
        self.write_artifacts()
        with self.assertRaisesRegex(ValueError, "Corresponding source mismatch"):
            verify(self.path)

    def test_missing_license_is_rejected(self):
        del self.runtime_files["licenses/LGPL-2.1.txt"]
        self.write_artifacts()
        with self.assertRaises(KeyError):
            verify(self.path)


if __name__ == "__main__":
    unittest.main()
