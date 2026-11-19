#!/usr/bin/env python3
"""Audit GetFEM 5.5 Savannah tarball against the immutable official Git tag.

This does not modify the SCRIMP Dockerfile, caches, sources or images.
It fetches into a disposable temporary directory, checks the archive
hash and gzip integrity, compares overlapping Git-tracked files by
*content* (not gzip bytes or tar timestamps), and applies all five
approved GetFEM patches to that disposable extracted source tree.

This is a SOURCE audit, not a cryptographic signature from upstream:
a valid gzip file, and even equal sources in selected paths, do not
by themselves authenticate a release artifact. Review differences
before changing the pinned Docker SHA256.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
import tarfile
import tempfile

REPO = Path(__file__).resolve().parents[1]
SOURCE_URL = (
    "https://download-mirror.savannah.gnu.org/releases/getfem/stable/"
    "getfem-5.5.tar.gz"
)
GIT_MIRROR = "https://github.com/getfem/getfem.git"
TAG = "v5.5"
TAG_COMMIT = "c242f7fe0989e9a90846c63eaccc6b495d774724"
KNOWN_SAVANNAH_SHA = "0e042770b86152b09be4bfd0ee160b8d10d49c3f8a0b9b06d55f1ffaa1d30fc4"
ORIGINAL_FOSSIES_SHA = "f6c09fb2edd112b6abe04e0485534f5e965d6f8f9ad853308b7120bfc6722ab5"

CORE_ROOTS = ("src/", "interface/", "tests/", "contrib/", "config/", "m4/")
CORE_TOP = {
    "configure.ac", "Makefile.am", "autogen.sh",
    "src/Makefile.am", "interface/Makefile.am",
}


def run(argv: list[str], **kwargs: object) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(argv, check=True, **kwargs)  # type: ignore[arg-type]


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def is_core(path: str) -> bool:
    return path.startswith(CORE_ROOTS) or path in CORE_TOP


def audit(archive: Path, checkout: Path, extraction: Path) -> int:
    print("\n=== Official reference checkout ===", flush=True)
    checkout.mkdir()
    run(["git", "-C", str(checkout), "init", "-q"])
    run(["git", "-C", str(checkout), "remote", "add", "origin", GIT_MIRROR])
    run(["git", "-C", str(checkout), "fetch", "--quiet", "--depth", "1",
         "--no-tags", "origin", f"refs/tags/{TAG}"])
    run(["git", "-C", str(checkout), "checkout", "--detach", "--quiet",
         "FETCH_HEAD"])
    actual_commit = run(["git", "-C", str(checkout), "rev-parse", "HEAD"],
                        stdout=subprocess.PIPE).stdout.decode().strip()
    print(f"Tag {TAG} resolved to {actual_commit}")
    if actual_commit != TAG_COMMIT:
        raise RuntimeError(
            f"Official tag unexpectedly moved: expected {TAG_COMMIT}."
        )

    print("\n=== Compare archive entries to tracked tag files ===", flush=True)
    # Inspect member names before extracting anything: never trust
    # a tarball filename without checking absolute/parent paths.
    archive_members: dict[str, tarfile.TarInfo] = {}
    with tarfile.open(archive, "r:gz") as tf:
        for member in tf.getmembers():
            name = PurePosixPath(member.name)
            if name.is_absolute() or ".." in name.parts:
                raise RuntimeError(f"Unsafe release archive member: {member.name}")
            if len(name.parts) == 0 or name.parts[0] != "getfem-5.5":
                raise RuntimeError(
                    f"Unexpected GetFEM top-level archive member: {member.name}"
                )
            if member.isfile():
                relative = str(PurePosixPath(*name.parts[1:]))
                if relative in archive_members:
                    raise RuntimeError(f"Duplicate archive file: {relative}")
                archive_members[relative] = member

        tracked_bytes = run(
            ["git", "-C", str(checkout), "ls-files", "-z"],
            stdout=subprocess.PIPE,
        ).stdout
        tracked = [os.fsdecode(p) for p in tracked_bytes.split(b"\0") if p]
        common = 0
        missing: list[str] = []
        different: list[str] = []
        missing_core: list[str] = []
        different_core: list[str] = []
        for relative in tracked:
            path = checkout / relative
            info = archive_members.get(relative)
            if info is None:
                missing.append(relative)
                if is_core(relative) and path.is_file() and not path.is_symlink():
                    missing_core.append(relative)
                continue
            if not path.is_file() or path.is_symlink():
                continue
            stream = tf.extractfile(info)
            if stream is None:
                raise RuntimeError(f"Could not read tar member: {relative}")
            with stream:
                original = stream.read()
            common += 1
            if original != path.read_bytes():
                different.append(relative)
                if is_core(relative):
                    different_core.append(relative)

    print(f"Git-tracked paths in tag:     {len(tracked)}")
    print(f"Regular files in tarball:    {len(archive_members)}")
    print(f"Common byte-compared files:  {common}")
    print(f"Git paths absent from tar:   {len(missing)} "
          "(some VCS/generated/dist omissions expected)")
    print(f"Different common files:      {len(different)}")
    print(f"Missing core source paths:   {len(missing_core)}")
    print(f"Different core source paths: {len(different_core)}")
    if different:
        print("DIFFERENT examples (first 25):")
        for p in different[:25]:
            print("  " + p)
    if missing_core:
        print("MISSING CORE examples (first 25):")
        for p in missing_core[:25]:
            print("  " + p)

    if not common or different_core or missing_core:
        print("FAIL: release archive core files are not byte-identical to "
              "the pinned Git tag. Do not change the Docker SHA256.")
        return 2

    print("PASS: all shared core source files match the official Git tag.")
    if different:
        print("NOTE: some shared non-core files differ; review above before "
              "claiming a fully identical release tree.")

    print("\n=== Approved patches, sequential disposable application ===",
          flush=True)
    extraction.mkdir()
    # The tar member names were checked above. Extraction occurs only
    # into our disposable temporary workdir, never into the SCRIMP tree.
    run(["tar", "--no-same-owner", "-xzf", str(archive), "-C",
         str(extraction)])
    src = extraction / "getfem-5.5"
    patches = sorted((REPO / "scripts/patches/getfem-5.5").glob("*.patch"))
    if len(patches) != 5:
        raise RuntimeError(f"Expected five approved GetFEM patches; got "
                           f"{len(patches)}: {[p.name for p in patches]}")
    for patch in patches:
        print(f"  {patch.name}", flush=True)
        with patch.open("rb") as f:
            subprocess.run(["patch", "--batch", "--forward", "-p1",
                            "--reject-file=-", "-d", str(src)],
                           stdin=f, check=True, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT)
    print("PASS: approved 0000..0004 patches apply sequentially to "
          "the Savannah GetFEM tarball.")
    print("\nSOURCE-AUDIT PASS: safe to REVIEW the observed archive digest "
          "for a new strict release pin.")
    print("No Docker image or SCRIMP file was changed.")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path,
                        help="optional predownloaded Savannah getfem-5.5.tar.gz")
    args = parser.parse_args()
    print("Official Git reference:", GIT_MIRROR, TAG, TAG_COMMIT)
    print("Pinned historical Fossies digest:", ORIGINAL_FOSSIES_SHA)
    print("Savannah digest observed on SCRIMP host:", KNOWN_SAVANNAH_SHA)
    print("Changes to SCRIMP Docker pin: NONE")

    with tempfile.TemporaryDirectory(prefix="scrimp-getfem55-audit-") as temp:
        root = Path(temp)
        archive = root / "getfem-5.5.tar.gz"
        if args.archive is not None:
            shutil.copyfile(args.archive.expanduser().resolve(), archive)
        else:
            print("Downloading:", SOURCE_URL, flush=True)
            run(["curl", "--fail", "--location", "--retry", "3",
                 "--output", str(archive), SOURCE_URL])

        actual_digest = digest(archive)
        print("Archive bytes:", archive.stat().st_size)
        print("Archive SHA256:", actual_digest)
        if actual_digest != KNOWN_SAVANNAH_SHA:
            print("FAIL: observed source archive SHA256 changed again. "
                  "Review the mirror, do not pin an unexpected artifact.")
            return 2
        with gzip.open(archive, "rb") as gz:
            # Complete decompression tests gzip CRC and verifies that
            # the archive isn't merely a successful HTTP response.
            while gz.read(1024 * 1024):
                pass
        print("PASS: expected Savannah SHA and gzip integrity.")
        return audit(archive, root / "tag", root / "extracted")


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, subprocess.CalledProcessError,
            EOFError, tarfile.TarError) as err:
        print(f"ERROR: GetFEM release audit failed: {err}", file=sys.stderr)
        sys.exit(2)
