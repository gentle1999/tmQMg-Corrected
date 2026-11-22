"""Build the SDF-only release from verified shards and include attribution."""

import argparse
import gzip
import hashlib
import io
import json
import shutil
import subprocess
import tarfile
from pathlib import Path
from tempfile import NamedTemporaryFile, TemporaryDirectory

from dataset import ROOT

NOTICES = (
    "sources/DATA-LICENSE.md",
    "sources/CC-BY-4.0.txt",
    "sources/tmQMg-LICENSE.txt",
    "sources/xyz2mol_tm-LICENSE.txt",
)


def add(archive, path, content):
    entry = tarfile.TarInfo(path)
    entry.size = len(content)
    archive.addfile(entry, io.BytesIO(content))


def build_sdf(manifest):
    output = ROOT / "dist/tmqmg-corrected.tar.gz"
    output.parent.mkdir(exist_ok=True)
    with NamedTemporaryFile(dir=output.parent, delete=False) as tmp:
        temporary = Path(tmp.name)
    try:
        with (
            temporary.open("wb") as raw,
            gzip.GzipFile(
                filename="", fileobj=raw, mode="wb", compresslevel=6, mtime=0
            ) as compressed,
            tarfile.open(fileobj=compressed, mode="w|") as archive,
        ):
            sums = []
            for part in manifest["parts"]:
                name = part["sdf"]
                content = (ROOT / name).read_bytes()
                digest = hashlib.sha256(content).hexdigest()
                if digest != manifest["sha256"][name]:
                    raise ValueError(f"Unexpected SHA-256 for {name}")
                add(archive, name, content)
                sums.append(f"{digest}  {name}\n")
            add(archive, "SHA256SUMS", "".join(sums).encode())
            add(archive, "VERSION", (manifest["version"] + "\n").encode())
            for name in NOTICES:
                add(archive, name, (ROOT / name).read_bytes())
        temporary.replace(output)
    finally:
        temporary.unlink(missing_ok=True)
    print(f"{output} ({manifest['targets']:,} graphs)")
    return output


def digest(path):
    sha = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            sha.update(chunk)
    return sha.hexdigest()


def build_bundle(manifest, sdf):
    version = manifest["version"]
    stem = f"tmqmg-corrected-{version}"
    tag = f"v{version}"
    git = ["git", "-C", str(ROOT)]
    commit = subprocess.check_output(
        [*git, "rev-parse", f"refs/tags/{tag}^{{commit}}"], text=True
    ).strip()
    subprocess.run([*git, "diff", "--quiet", commit, "--"], check=True)
    log_dir = sdf.parent / f"tmqmg-corrected-gxtb-logs-{version}"
    target = sdf.parent / stem
    with TemporaryDirectory(dir=sdf.parent) as tmp:
        bundle = Path(tmp) / stem
        bundle.mkdir()
        snapshot = bundle / f"{stem}-snapshot.tar.gz"
        subprocess.run(
            [
                *git,
                "archive",
                "--format=tar.gz",
                f"--prefix={stem}/",
                "-o",
                str(snapshot),
                tag,
            ],
            check=True,
        )
        structure_file = bundle / f"{stem}-sdf.tar.gz"
        shutil.copyfile(sdf, structure_file)
        for name, expected in manifest["gxtb_log_files_sha256"].items():
            path = bundle / name
            shutil.copyfile(log_dir / name, path)
            if digest(path) != expected:
                raise ValueError(f"Changed release file: {name}")
        log_name = log_dir.name
        (bundle / "README.md").write_text(
            f"# tmQMg-Corrected {version}\n\nGit tag `{tag}`, commit `{commit}`.\n\n"
            f"- [Lewis graphs in SDF]({structure_file.name})\n"
            f"- [Dataset, evidence and reproduction scripts]({snapshot.name})\n"
            f"- [Original g-xTB logs]({log_name}.zip)\n"
            f"- [Log index]({log_name}-index.csv)\n\n"
            "Run `sha256sum -c SHA256SUMS` in this directory to verify the files.\n"
        )
        sums = "".join(f"{digest(p)}  {p.name}\n" for p in sorted(bundle.iterdir()))
        (bundle / "SHA256SUMS").write_text(sums)
        if target.exists():
            if {p.name for p in target.iterdir()} != {
                p.name for p in bundle.iterdir()
            } or any(digest(target / p.name) != digest(p) for p in bundle.iterdir()):
                raise ValueError(f"Existing version bundle differs: {target}")
        else:
            bundle.rename(target)
    print(target)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--bundle",
        action="store_true",
        help="assemble the tagged snapshot, SDF, original logs and index",
    )
    manifest = json.loads((ROOT / "manifest.json").read_bytes())
    args = parser.parse_args()
    sdf = build_sdf(manifest)
    if args.bundle:
        build_bundle(manifest, sdf)


if __name__ == "__main__":
    main()
