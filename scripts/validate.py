"""Verify the tmQMg-Corrected release, optionally against the pinned upstream files."""

import argparse
import base64
import csv
import hashlib
import json
import math
import zlib
from collections import Counter, defaultdict
from contextlib import nullcontext
from itertools import islice
from pathlib import Path
from zipfile import ZipFile

from dataset import ROOT
from exclusions import (
    BRIDGE_AUDIT,
    INDEX_FIELDS,
    boron_assessment,
    bridge_assessment,
    element_counts,
    formula,
    index_entry,
    load_bridge_audit,
    spin_reassessment,
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def file_digest(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_upstream(path, digest):
    require(file_digest(path) == digest, f"Unexpected source CSV digest: {path}")
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return {
            r["id"]: (r["smiles"], int(float(r["charge"])))
            for r in csv.DictReader(stream)
        }


def check_xyz(stream, row, archive):
    case = row["id"]
    count = stream.readline()
    comment = stream.readline().decode().strip()
    expected = f"id={case} charge={row['charge']} multiplicity={row['multiplicity']}"
    require(comment == expected, f"{case}: XYZ header differs from metadata")
    atoms = list(islice(stream, int(count)))
    require(len(atoms) == int(count), f"{case}: incomplete XYZ")
    original = count + b"\n" + b"".join(atoms)
    require(
        hashlib.sha256(original).hexdigest() == row["xyz_sha256"], f"{case}: XYZ digest"
    )
    if archive is not None:
        require(
            archive.read(f"xyz/{case}.xyz") == original, f"{case}: upstream XYZ differs"
        )
    return atoms


NONMETALS = {
    "H",
    "He",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Po",
    "At",
    "Rn",
}


def controls(raw):
    return {
        key: row["q"]
        for key, row in raw.items()
        if row["cases"] >= 15 and row["bonds"] >= 40
    }


def assessment(value, reference, other=()):
    if value is None:
        return "Unreported"
    if reference is None:
        return "Sparse controls"
    if not reference[0] <= value <= reference[-1]:
        return "Outside range"
    if not reference[1] <= value <= reference[-2] and any(
        lo <= value <= hi for lo, hi in other
    ):
        return "Closer to another bond class"
    return "Within range"


def check_run(row, run):
    case = row["id"]
    require(
        run["status"] == "converged"
        and run["charge"] == row["charge"]
        and run["uhf"] == row["multiplicity"] - 1,
        f"{case}: archived single-point state",
    )
    # The pinned Fortran build prints the archived WBO range with 15 significant
    # digits. Keeping row order and trailing spaces reproduces the original file.
    lines = []
    for i, j, value in run["wbo_pairs"]:
        decimals = 15 if value < 1 else 14
        lines.append(f"{i:12d}{j:12d}{value:19.{decimals}f}     \n")
    digest = hashlib.sha256("".join(lines).encode("ascii")).hexdigest()
    require(digest == run["wbo_sha256"], f"{case}: {run['cohort']} raw WBO digest")


def check_bonds(row, atoms, edges, calibration):
    geometry, wbo, alternatives, threshold = calibration
    case = row["id"]
    audit = row["wbo_analysis"]
    run = next(g for g in row["gxtb"] if g["cohort"] == audit["cohort"])
    values = {tuple(sorted((i, j))): v for i, j, v in run["wbo_pairs"]}
    symbols = [a[0].decode() for a in atoms]
    positions = [tuple(map(float, a[1:])) for a in atoms]
    expected = {
        (i, j)
        for (i, j), order in edges.items()
        if order in (1, 2, 3, 4)
        and symbols[i - 1] in NONMETALS
        and symbols[j - 1] in NONMETALS
    }
    counts: Counter[str] = Counter()
    for i, j, pair, order, control, length, length_label, value, value_label in audit[
        "bonds"
    ]:
        require((i, j) in expected, f"{case}: surplus bond {i}-{j}")
        expected.remove((i, j))
        distance = math.dist(positions[i - 1], positions[j - 1])
        raw = values.get((i, j))
        key = f"{pair}/{control}"
        require(
            edges[i, j] == order
            and pair == "-".join(sorted((symbols[i - 1], symbols[j - 1])))
            and round(distance, 6) == length
            and (None if raw is None else round(raw, 8)) == value
            and assessment(distance, geometry.get(key)) == length_label
            and assessment(raw, wbo.get(key), alternatives.get(key, ())) == value_label,
            f"{case}: bond evidence {i}-{j}",
        )
        counts[value_label] += 1
    require(
        not expected and counts == Counter(audit["bond_counts"]),
        f"{case}: bond inventory",
    )
    contacts = sorted(
        [i, j, round(v, 8)]
        for (i, j), v in values.items()
        if v > threshold
        and (i, j) not in edges
        and symbols[i - 1] in NONMETALS
        and symbols[j - 1] in NONMETALS
    )
    require(
        threshold == audit["contact_threshold"]
        and contacts == sorted(audit["unlisted_contacts"]),
        f"{case}: unlisted contacts",
    )


def check_sdf(stream, row, atoms, calibration):
    case = row["id"]
    require(
        stream.readline().decode().rstrip() == case,
        f"{case}: SDF title or record order",
    )
    record = []
    for line in stream:
        if line == b"$$$$\n":
            break
        record.append(line)
    else:
        raise ValueError(f"{case}: incomplete SDF")
    start = record.index(b"M  V30 BEGIN ATOM\n") + 1
    end = record.index(b"M  V30 END ATOM\n", start)
    require(end - start == len(atoms), f"{case}: atom count")
    for number, (xyz_line, sdf_line) in enumerate(
        zip(atoms, record[start:end], strict=True), start=1
    ):
        symbol, *coords = xyz_line.split()
        fields = sdf_line.split()
        require(
            int(fields[2]) == number
            and fields[3] == symbol
            and list(map(float, fields[4:7])) == list(map(float, coords)),
            f"{case}: atom {number} differs from XYZ",
        )
    start = record.index(b"M  V30 BEGIN BOND\n") + 1
    end = record.index(b"M  V30 END BOND\n", start)
    edges = {}
    for line in record[start:end]:
        order, i, j = map(int, line.split()[3:6])
        edges[tuple(sorted((i, j)))] = order
    check_bonds(row, [line.split() for line in atoms], edges, calibration)
    props = {
        line.strip()[4:-1].decode(): record[i + 1].decode().rstrip("\n")
        for i, line in enumerate(record)
        if line.startswith(b">  <")
    }
    native = zlib.decompress(
        base64.b64decode(props["NATIVE_GRAPH_ZLIB_BASE64"], validate=True)
    )
    require(
        hashlib.sha256(native).hexdigest() == row["graph_sha256"],
        f"{case}: native graph digest",
    )
    require(
        props.get("SOURCE_NATIVE_GRAPH_SHA256")
        == row.get("source_native_graph_sha256"),
        f"{case}: source native graph digest differs between SDF and metadata",
    )
    for prop, key in (
        ("GRAPH_SOURCE", "graph_source"),
        ("SELECTION_REASON", "selection_reason"),
        ("UPSTREAM_COMPARISON", "original_graph_comparison"),
        ("SPIN_REVIEW", "spin_review"),
        ("CASE_ID", "id"),
        ("GRAPH_SHA256", "graph_sha256"),
        ("XYZ_SHA256", "xyz_sha256"),
        ("EVIDENCE_TIER", "tier"),
        ("INPUT_CHARGE", "charge"),
        ("INPUT_MULTIPLICITY", "multiplicity"),
    ):
        require(
            props[prop] == " ".join(str(row[key]).splitlines()),
            f"{case}: {prop} differs from metadata",
        )


def check_exclusion(row, atoms, part, index, bridges):
    case = row["id"]
    assessment = (row.get("exclusion_evidence") or {}).get("assessment")
    if assessment is None:
        raise ValueError(f"{case}: missing exclusion assessment")
    if not isinstance(assessment, dict):
        raise TypeError(f"{case}: exclusion assessment must be a JSON object")
    require(
        assessment["reason"]
        and assessment["basis"]
        and assessment["sources"]
        and row["selection_reason"]
        == assessment["reason"] + ". " + assessment["basis"],
        f"{case}: incomplete exclusion reason, basis or references",
    )
    require(
        index.pop(case, None) == index_entry(row, part),
        f"{case}: exclusion index mismatch",
    )
    tokens = [line.decode().split() for line in atoms]
    counts = Counter(atom[0] for atom in tokens)
    require(formula(counts) == assessment["input_formula"], f"{case}: input formula")
    if row["exclusion"] == "Missing input hydrogens":
        expected = element_counts(assessment["source_formula"])
        delta = assessment["missing_H"]
        require(
            delta > 0
            and not counts - expected
            and expected - counts == {"H": delta}
            and counts["H"] == assessment["input_H"]
            and expected["H"] == assessment["expected_H"],
            f"{case}: hydrogen evidence inconsistent with XYZ",
        )
    elif row["exclusion"] == "Boron hydrogen bridge":
        require(
            assessment == bridge_assessment(bridges[case]),
            f"{case}: bridge evidence mismatch",
        )
    elif row["exclusion"] == "Boron cluster":
        reason, basis, measured = boron_assessment(case, tokens)
        require(
            (assessment["reason"], assessment["basis"]) == (reason, basis)
            and all(assessment[k] == v for k, v in measured.items()),
            f"{case}: boron evidence inconsistent with coordinates",
        )
    elif row["exclusion"] == "Questionable electronic state":
        if "spin_reassessment" in row["exclusion_evidence"]:
            require(
                assessment == spin_reassessment(row, counts),
                f"{case}: alternative spin assessment differs from its evidence",
            )
            return
        require(
            assessment["input_multiplicity"] == row["multiplicity"]
            and assessment["reported_state"]
            and any(
                url.startswith("https://doi.org/") for url in assessment["sources"]
            ),
            f"{case}: missing input state or literature evidence",
        )
    else:
        raise ValueError(f"{case}: unknown exclusion")


class ControlFit:
    """Rebuild both empirical calibrations and their independent holdout counts."""

    def __init__(self):
        self.distances = defaultdict(list)
        self.wbos = defaultdict(list)
        self.distance_cases: Counter[str] = Counter()
        self.wbo_cases: Counter[str] = Counter()
        self.max_contacts = []
        self.holdout = []

    def add(self, row, atoms, run):
        case = row["id"]
        holdout = (
            int(hashlib.sha256(("holdout" + case).encode()).hexdigest()[:6], 16) % 5
            == 0
        )
        positions = [tuple(map(float, line.split()[1:])) for line in atoms]
        symbols = [line.split()[0].decode() for line in atoms]
        pairs = {tuple(sorted((i, j))): v for i, j, v in run["wbo_pairs"]}
        lengths, wbos, bonds = set(), set(), []
        audit = row.get("wbo_analysis") or row["previous_graph_review"]["wbo_analysis"]
        for i, j, pair, _, control, *_ in audit["bonds"]:
            key = f"{pair}/{control}"
            value = pairs.pop((i, j), None)
            if holdout:
                bonds.append((key, value))
            else:
                self.distances[key].append(
                    math.dist(positions[i - 1], positions[j - 1])
                )
                lengths.add(key)
                if value is not None:
                    self.wbos[key].append(value)
                    wbos.add(key)
        maximum = max(
            (
                v
                for (i, j), v in pairs.items()
                if symbols[i - 1] in NONMETALS and symbols[j - 1] in NONMETALS
            ),
            default=0.0,
        )
        if holdout:
            self.holdout.append((bonds, maximum))
        else:
            self.max_contacts.append(maximum)
            self.distance_cases.update(lengths)
            self.wbo_cases.update(wbos)

    @staticmethod
    def quantile(sorted_values, fraction):
        rank = (len(sorted_values) - 1) * fraction
        low = math.floor(rank)
        a, b = sorted_values[low], sorted_values[math.ceil(rank)]
        return a + (rank - low) * (b - a)

    def verify(self, manifest, calibration):
        spec = manifest["wbo_calibration"]
        for observations, counts, archived in (
            (self.distances, self.distance_cases, manifest["geometry_calibration"]),
            (self.wbos, self.wbo_cases, spec["ranges"]),
        ):
            require(observations.keys() == archived.keys(), "Control classes differ")
            for key, sample in observations.items():
                sample.sort()
                expected = archived[key]
                require(
                    len(sample) == expected["bonds"]
                    and counts[key] == expected["cases"],
                    f"{key}: control sample size",
                )
                quantiles = [self.quantile(sample, f) for f in spec["range_quantiles"]]
                require(
                    len(quantiles) == len(expected["q"])
                    and all(
                        math.isclose(a, b, abs_tol=1e-10, rel_tol=0)
                        for a, b in zip(quantiles, expected["q"], strict=True)
                    ),
                    f"{key}: control quantiles",
                )
        _, wbo, alternatives, threshold = calibration
        require(
            len(self.max_contacts) == spec["fit_cases"]
            and len(self.holdout) == spec["validation_cases"],
            "Control partition sizes",
        )
        require(
            math.isclose(
                self.quantile(sorted(self.max_contacts), 0.99),
                threshold,
                rel_tol=0,
                abs_tol=1e-12,
            ),
            "Control off-graph contact threshold",
        )
        labels = {
            "Within range": "supported",
            "Outside range": "outside_empirical_range",
            "Closer to another bond class": "alternative_order_favored",
            "Sparse controls": "sparse_reference",
            "Unreported": "not_reported",
        }
        totals: Counter[str] = Counter()
        for bonds, maximum in self.holdout:
            statuses = [
                labels[assessment(v, wbo.get(k), alternatives.get(k, ()))]
                for k, v in bonds
            ]
            totals.update(statuses)
            if set(statuses) == {"supported"} and maximum <= threshold:
                totals["all_bonds_supported_and_no_extra"] += 1
        require(totals == Counter(spec["holdout_results"]), "Holdout assessment counts")


def validate(manifest, source_rows, archive):
    totals: Counter[str] = Counter()
    tiers: Counter[str] = Counter()
    exclusions: Counter[str] = Counter()
    wbo = controls(manifest["wbo_calibration"]["ranges"])
    alternatives = {
        key: [
            (q[2], q[4])
            for other, q in wbo.items()
            if other != key and other.split("/")[0] == key.split("/")[0]
        ]
        for key in wbo
    }
    calibration = (
        controls(manifest["geometry_calibration"]),
        wbo,
        alternatives,
        manifest["wbo_calibration"]["max_unlisted_threshold"],
    )
    control_fit = ControlFit()
    files = set()
    with (ROOT / manifest["exclusion_index"]).open(
        encoding="utf-8", newline=""
    ) as stream:
        table = csv.DictReader(stream)
        require(
            tuple(table.fieldnames or ()) == INDEX_FIELDS,
            "Unexpected exclusion index fields",
        )
        index = {}
        for entry in table:
            require(
                entry["id"] not in index, f"Duplicate exclusion index ID: {entry['id']}"
            )
            index[entry["id"]] = entry
    require(
        file_digest(ROOT / BRIDGE_AUDIT) == manifest["boron_bridge_audit_sha256"],
        "Bridge index digest",
    )
    require(manifest["boron_bridge_audit"] == BRIDGE_AUDIT, "Bridge index path")
    bridges = load_bridge_audit()
    unseen_bridge_ids = set(bridges)
    previous = ""
    for number, part in enumerate(manifest["parts"], 1):
        for key, extension in (("xyz", "xyz"), ("sdf", "sdf"), ("metadata", "jsonl")):
            name = part[key]
            directory = "data" if key == "sdf" else "evidence/inputs"
            require(
                name == f"{directory}/tmqmg-corrected-{number:02}.{extension}",
                f"Unexpected {name}",
            )
            require(
                file_digest(ROOT / name) == manifest["sha256"][name], f"{name}: SHA-256"
            )
            files.add(name)
        inputs = targets = 0
        with (
            (ROOT / part["xyz"]).open("rb") as xyz,
            (ROOT / part["sdf"]).open("rb") as sdf,
            (ROOT / part["metadata"]).open(encoding="utf-8") as metadata,
        ):
            for line in metadata:
                row = json.loads(line)
                case = row["id"]
                unseen_bridge_ids.discard(case)
                require(case > previous, f"{case}: duplicate or unsorted ID")
                if not inputs:
                    require(case == part["first"], f"{case}: first ID")
                previous = case
                inputs += 1
                atoms = check_xyz(xyz, row, archive)
                if source_rows is not None:
                    require(
                        source_rows.pop(case)
                        == (row["original_smiles"] or "", row["charge"]),
                        f"{case}: upstream SMILES or charge differs",
                    )
                totals["gxtb_runs"] += len(row["gxtb"])
                totals["gxtb_cases"] += bool(row["gxtb"])
                for run in row["gxtb"]:
                    check_run(row, run)
                    if run["cohort"] == "control":
                        control_fit.add(row, atoms, run)
                if row["benchmark"]:
                    require(
                        row["exclusion"] is None, f"{case}: target has an exclusion"
                    )
                    require(
                        case not in bridges
                        or bridges[case]["evidence"] == "geometry_only",
                        f"{case}: boron hydrogen bridge cannot be a scoring target",
                    )
                    check_sdf(sdf, row, atoms, calibration)
                    tiers[str(row["tier"])] += 1
                    targets += 1
                else:
                    require(
                        row["tier"] is None,
                        f"{case}: excluded entry has an evidence tier",
                    )
                    check_exclusion(row, atoms, part, index, bridges)
                    exclusions[row["exclusion"]] += 1
            require(
                not xyz.read(1) and not sdf.read(1), f"Group {number}: surplus records"
            )
        require(
            (previous, inputs, targets)
            == (part["last"], part["inputs"], part["targets"]),
            f"Group {number}: ID range or counts differ from manifest",
        )
        totals["inputs"] += inputs
        totals["targets"] += targets
    control_fit.verify(manifest, calibration)
    require(len(manifest["parts"]) == 50, "Expected 50 groups")
    require(files == set(manifest["sha256"]), "File list differs from checksum list")
    require(
        files
        == {
            p.relative_to(ROOT).as_posix()
            for directory in ("data", "evidence/inputs")
            for p in (ROOT / directory).iterdir()
        },
        "Unexpected files in data/ or evidence/inputs/",
    )
    require(not unseen_bridge_ids, "Bridge index contains unknown case IDs")
    require(not index, "Surplus entries in the exclusion index")
    require(not source_rows, "Unused upstream CSV entries")
    require(dict(totals) == {k: manifest[k] for k in totals}, "Release totals differ")
    require(
        tiers == manifest["tiers"] and exclusions == manifest["exclusions"],
        "Evidence tier or exclusion counts differ",
    )
    print(
        f"Verified {len(files)} files, {totals['inputs']:,} inputs, "
        f"{totals['targets']:,} targets and {sum(exclusions.values()):,} exclusions."
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--upstream-csv", type=Path, help="Pinned tmQMg properties and targets CSV"
    )
    parser.add_argument(
        "--upstream-xyz", type=Path, help="Pinned tmQMg XYZ ZIP archive"
    )
    args = parser.parse_args()
    manifest = json.loads((ROOT / "manifest.json").read_bytes())
    rows = None
    if args.upstream_csv:
        rows = read_upstream(args.upstream_csv, manifest["source_csv_sha256"])
    if args.upstream_xyz:
        require(
            file_digest(args.upstream_xyz) == manifest["source_xyz_zip_sha256"],
            f"Unexpected source XYZ archive digest: {args.upstream_xyz}",
        )
    with ZipFile(args.upstream_xyz) if args.upstream_xyz else nullcontext() as archive:
        validate(manifest, rows, archive)


if __name__ == "__main__":
    main()
