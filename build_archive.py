#!/usr/bin/env python3
"""Rebuild the certificate manifest and ZIP with reproducible metadata."""

from __future__ import annotations

import hashlib
import stat
import zipfile
from pathlib import Path


PACKAGE = Path(__file__).resolve().parent
ARCHIVE = PACKAGE.parent / f"{PACKAGE.name}.zip"
MANIFEST = PACKAGE / "SHA256SUMS"
FIXED_TIME = (1980, 1, 1, 0, 0, 0)
EXCLUDED_PARTS = {"__pycache__"}
PROOF_ARTIFACTS = (
    "evidence/launch/all_s_s4_coupling_identity_exact.py",
    "evidence/launch/all_s_s4_explicit_counterexample_checker.py",
    "evidence/launch/all_s_s4_intrinsic_lift_exact.py",
    "evidence/launch/all_s_s4_local_majorants_exact.py",
    "evidence/launch/all_s_s4_scalar_hopf_vector_exact.py",
    "evidence/launch/all_s_s4_scalar_to_raw_tangent_exact.py",
    "evidence/launch/all_s_s4_transverse_hopf_exact.py",
)


def payload_files() -> list[Path]:
    files = []
    for path in PACKAGE.rglob("*"):
        if not path.is_file() or path == MANIFEST:
            continue
        if any(part in EXCLUDED_PARTS for part in path.relative_to(PACKAGE).parts):
            continue
        files.append(path)
    return sorted(files, key=lambda path: path.relative_to(PACKAGE).as_posix())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_manifest() -> None:
    files = [PACKAGE / relative for relative in PROOF_ARTIFACTS]
    lines = [
        f"{sha256(path)}  {path.relative_to(PACKAGE).as_posix()}"
        for path in files
    ]
    MANIFEST.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def zip_info(name: str, mode: int, is_directory: bool = False) -> zipfile.ZipInfo:
    if is_directory and not name.endswith("/"):
        name += "/"
    info = zipfile.ZipInfo(name, FIXED_TIME)
    info.create_system = 3
    kind = stat.S_IFDIR if is_directory else stat.S_IFREG
    info.external_attr = (kind | mode) << 16
    info.compress_type = zipfile.ZIP_DEFLATED
    return info


def build_archive(files: list[Path]) -> None:
    all_files = sorted(files + [MANIFEST], key=lambda path: path.relative_to(PACKAGE).as_posix())
    directories = {PACKAGE.name}
    for path in all_files:
        relative = path.relative_to(PACKAGE)
        for parent in relative.parents:
            if parent != Path("."):
                directories.add(f"{PACKAGE.name}/{parent.as_posix()}")

    temporary = ARCHIVE.with_suffix(".zip.tmp")
    with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for directory in sorted(directories):
            archive.writestr(zip_info(directory, 0o755, is_directory=True), b"")
        for path in all_files:
            relative = path.relative_to(PACKAGE).as_posix()
            mode = 0o755 if path.suffix == ".sh" or path.name == Path(__file__).name else 0o644
            archive.writestr(zip_info(f"{PACKAGE.name}/{relative}", mode), path.read_bytes())
    temporary.replace(ARCHIVE)


def main() -> None:
    files = payload_files()
    write_manifest()
    build_archive(files)
    print(f"wrote {MANIFEST}")
    print(f"wrote {ARCHIVE}")


if __name__ == "__main__":
    main()
