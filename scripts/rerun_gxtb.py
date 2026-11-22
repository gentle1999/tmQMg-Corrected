"""Repeat archived g-xTB single points on the original XYZ and compare the outputs."""

import argparse
import csv
import hashlib
import itertools
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from tempfile import TemporaryDirectory

from dataset import ROOT


def run_case(task):
    executable, environment, case, xyz, run = task
    with TemporaryDirectory(prefix="tmqmg-gxtb-") as folder:
        work = Path(folder)
        (work / "input.xyz").write_bytes(xyz)
        options = run["options"].split()
        process = subprocess.run(
            [executable, "input.xyz", *options],
            check=False,
            cwd=work,
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        log = process.stdout
        if process.returncode or b"normal termination of xtb" not in log:
            raise ValueError(
                f"{case} ({run['cohort']}): calculation failed\n{log[-2000:].decode(errors='replace')}"
            )
        raw = (work / "wbo").read_bytes()
    pairs = {
        (int(i), int(j)): float(v) for i, j, v in map(bytes.split, raw.splitlines())
    }
    expected = {(i, j): v for i, j, v in run["wbo_pairs"]}
    digest = hashlib.sha256(raw).hexdigest()
    if pairs.keys() != expected.keys():
        raise ValueError(f"{case} ({run['cohort']}): WBO pair inventory differs")
    max_delta = max((abs(v - expected[k]) for k, v in pairs.items()), default=0)
    energy = re.search(rb"\| TOTAL ENERGY\s+([-\d.]+) Eh", log)
    gap = re.search(rb":: HOMO-LUMO gap\s+([-\d.]+) eV", log)
    if energy is None or gap is None:
        raise ValueError(
            f"{case} ({run['cohort']}): missing total energy or HOMO-LUMO gap"
        )
    return (
        case,
        run["cohort"],
        digest,
        digest == run["wbo_sha256"],
        max_delta,
        float(energy[1]) - run["total_energy_Eh"],
        float(gap[1]) - run["homo_lumo_gap_eV"],
    )


def tasks(part, executable, environment, requested):
    with (
        (ROOT / part["metadata"]).open() as meta,
        (ROOT / part["xyz"]).open("rb") as xyz,
    ):
        for row in map(json.loads, meta):
            count = next(xyz)
            case = row["id"]
            if next(xyz).split()[0] != f"id={case}".encode():
                raise ValueError(f"{case}: XYZ and metadata out of order")
            raw = count + b"\n" + b"".join(itertools.islice(xyz, int(count)))
            if requested and case not in requested:
                continue
            if hashlib.sha256(raw).hexdigest() != row["xyz_sha256"]:
                raise ValueError(f"{case}: XYZ digest")
            for run in row["gxtb"]:
                yield executable, environment, case, raw, run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--xtb", type=Path, required=True, help="g-xTB v2.0.0 xtb executable"
    )
    parser.add_argument(
        "--jobs", type=int, default=1, help="concurrent single-threaded runs"
    )
    parser.add_argument("--part", type=int, choices=range(1, 51), metavar="1..50")
    parser.add_argument(
        "ids", nargs="*", help="specific CSD IDs (default: all archived runs)"
    )
    args = parser.parse_args()
    requested = {i.upper() for i in args.ids}
    completed = set()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    manifest = json.loads((ROOT / "manifest.json").read_bytes())
    executable = args.xtb.resolve()
    if (
        hashlib.sha256(executable.read_bytes()).hexdigest()
        != manifest["gxtb_method"]["executable_sha256"]
    ):
        parser.error("executable differs from the pinned g-xTB build")
    env = dict(
        os.environ,
        OMP_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
        OMP_STACKSIZE="1G",
        XTBPATH=str(executable.parents[1] / "share/xtb"),
    )
    out = csv.writer(sys.stdout, lineterminator="\n")
    out.writerow(
        (
            "id",
            "cohort",
            "wbo_sha256",
            "wbo_bytes_match",
            "max_abs_wbo_delta",
            "energy_delta_Eh",
            "gap_delta_eV",
        )
    )
    parts = (
        manifest["parts"] if args.part is None else [manifest["parts"][args.part - 1]]
    )
    with ThreadPoolExecutor(args.jobs) as pool:
        for part in parts:
            if requested and not any(
                part["first"] <= case <= part["last"] for case in requested
            ):
                continue
            pending = tasks(part, executable, env, requested)
            while batch := list(itertools.islice(pending, args.jobs * 4)):
                for result in pool.map(run_case, batch):
                    out.writerow(result)
                    completed.add(result[0])
                sys.stdout.flush()

    if missing := requested - completed:
        parser.error(f"No archived run in the selection for {sorted(missing)}")


if __name__ == "__main__":
    main()
