"""Package the original converged g-xTB logs for a tmQMg-Corrected release."""

import argparse
import csv
import hashlib
import io
import json
import re
import sqlite3
import zlib
from collections import Counter
from contextlib import closing
from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

from dataset import ROOT


def sha(data):
    return hashlib.sha256(data).hexdigest()


def release_runs(manifest):
    expected = {}
    for part in manifest["parts"]:
        path = part["metadata"]
        digest = hashlib.sha256()
        with (ROOT / path).open("rb") as stream:
            for line in stream:
                digest.update(line)
                row = json.loads(line)
                for run in row["gxtb"]:
                    key = (row["id"], run["cohort"])
                    if key in expected or run["status"] != "converged":
                        raise ValueError(f"Unexpected release record: {key}")
                    expected[key] = (
                        row["xyz_sha256"],
                        run["wbo_sha256"],
                        run["options"],
                        run["total_energy_Eh"],
                        run["homo_lumo_gap_eV"],
                        run["charge"],
                        run["uhf"],
                    )
        if digest.hexdigest() != manifest["sha256"][path]:
            raise ValueError(f"Changed release metadata: {path}")
    if (len(expected), len({case for case, _ in expected})) != (
        manifest["gxtb_runs"],
        manifest["gxtb_cases"],
    ):
        raise ValueError("Run inventory differs from the release manifest")
    return expected


def original_logs(main_dir, additional_dir, edakup_dir):
    for directory in (main_dir, additional_dir):
        for db in sorted(directory.glob("gxtb-*.sqlite")):
            shard = int(db.stem.split("-")[1])
            cohort = (
                "additional_control"
                if directory == additional_dir
                else ("control" if shard >= 100 else "main")
            )
            # These are frozen SQLite shards. Reject outstanding WAL data.
            wal = db.with_name(db.name + "-wal")
            if wal.exists() and wal.stat().st_size:
                raise ValueError(f"Uncheckpointed database: {db}")
            with closing(
                sqlite3.connect(
                    db.resolve().as_uri() + "?mode=ro&immutable=1", uri=True
                )
            ) as connection:
                rows = connection.execute(
                    "SELECT id, xyz_sha256, charge, uhf, log_zlib, wbo_zlib FROM runs "
                    "WHERE status = 'converged' ORDER BY id"
                )
                for case, xyz_hash, charge, uhf, log, wbo in rows:
                    yield (
                        case,
                        cohort,
                        db.name,
                        xyz_hash,
                        charge,
                        uhf,
                        zlib.decompress(log),
                        zlib.decompress(wbo),
                    )
    yield (
        "EDAKUP",
        "main",
        "EDAKUP/xtb.log",
        sha((edakup_dir / "EDAKUP.xyz").read_bytes()),
        0,
        0,
        (edakup_dir / "xtb.log").read_bytes(),
        (edakup_dir / "wbo").read_bytes(),
    )


def check_log(case, log):
    command = re.search(rb"program call\s+: [^\r\n]*? (\S+\.xyz) (--gxtb[^\r\n]+)", log)
    energy = re.search(rb"\| TOTAL ENERGY\s+([-\d.]+) Eh", log)
    gap = re.search(rb":: HOMO-LUMO gap\s+([-\d.]+) eV", log)
    if not (command and energy and gap and b"normal termination of xtb" in log):
        raise ValueError(f"{case}: incomplete original log")
    if Path(command[1].decode()).name != f"{case}.xyz":
        raise ValueError(f"{case}: different XYZ filename in log")
    return command[2].decode().strip(), float(energy[1]), float(gap[1])


def add(archive, name, content):
    entry = ZipInfo(name)  # Fixed timestamp for reproducible archives.
    archive.writestr(entry, content, compress_type=ZIP_DEFLATED, compresslevel=6)


def write_bundle(directory, stem, version, expected, originals):
    zip_path = directory / (stem + ".zip")
    index_path = directory / (stem + "-index.csv")
    index = []
    with ZipFile(zip_path, "w") as archive:
        for case, cohort, shard, xyz, charge, uhf, log, wbo in originals:
            key = (case, cohort)
            wbo_hash = sha(wbo)
            options, energy, gap = check_log(case, log)
            if expected.pop(key, None) != (
                xyz,
                wbo_hash,
                options,
                energy,
                gap,
                charge,
                uhf,
            ):
                raise ValueError(
                    f"{key}: duplicate, unexpected or inconsistent original"
                )
            name = f"logs/{cohort}/{case}.log"
            add(archive, name, log)
            index.append((case, cohort, name, sha(log), shard, xyz, wbo_hash))
        if expected:
            raise ValueError(
                f"{len(expected)} missing originals, including {next(iter(expected))}"
            )
        with index_path.open("w", newline="", encoding="utf-8") as out:
            writer = csv.writer(out, lineterminator="\n")
            writer.writerow(
                (
                    "id",
                    "cohort",
                    "log_path",
                    "log_sha256",
                    "source_shard",
                    "xyz_sha256",
                    "wbo_sha256",
                )
            )
            writer.writerows(sorted(index))
        add(archive, "index.csv", index_path.read_bytes())
        add(archive, "VERSION", (version + "\n").encode())

    # Read back every compressed member before exposing the bundle.
    with ZipFile(zip_path) as archive:
        if (
            archive.read("VERSION").decode().strip() != version
            or archive.read("index.csv") != index_path.read_bytes()
        ):
            raise ValueError("Packaged metadata differs from the release")
        if len(archive.infolist()) != len(index) + 2:
            raise ValueError("Unexpected ZIP inventory")
        for row in csv.DictReader(io.StringIO(archive.read("index.csv").decode())):
            if sha(archive.read(row["log_path"])) != row["log_sha256"]:
                raise ValueError(f"Archived log checksum mismatch: {row['log_path']}")
    with (directory / "SHA256SUMS").open("w") as sums:
        for path in (zip_path, index_path):
            digest = hashlib.sha256()
            with path.open("rb") as stream:
                for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                    digest.update(chunk)
            print(f"{digest.hexdigest()}  {path.name}", file=sums)
    return Counter(row[1] for row in index)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--main-dir", type=Path, required=True, help="main and control SQLite shards"
    )
    parser.add_argument(
        "--additional-dir",
        type=Path,
        required=True,
        help="additional control SQLite shards",
    )
    parser.add_argument(
        "--edakup-dir", type=Path, required=True, help="original EDAKUP job directory"
    )
    parser.add_argument("--output-dir", type=Path, default=ROOT / "dist")
    args = parser.parse_args()
    manifest = json.loads((ROOT / "manifest.json").read_bytes())
    version = manifest.get("gxtb_log_version", manifest["version"])
    stem = f"tmqmg-corrected-gxtb-logs-{version}"
    args.output_dir.mkdir(parents=True, exist_ok=True)
    target = args.output_dir / stem
    if target.exists():
        parser.error(
            f"Bundle already exists at {target}. Choose a fresh --output-dir to rebuild."
        )
    with TemporaryDirectory(dir=args.output_dir) as tmp:
        bundle = Path(tmp) / stem
        bundle.mkdir()
        counts = write_bundle(
            bundle,
            stem,
            version,
            release_runs(manifest),
            original_logs(args.main_dir, args.additional_dir, args.edakup_dir),
        )
        bundle.rename(target)
    print(f"{target}: {sum(counts.values()):,} verified logs, {dict(counts)}")


if __name__ == "__main__":
    main()
