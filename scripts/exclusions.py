"""Reconstruct the case evidence for withheld reference graphs from XYZ and cited records."""

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from collections import Counter
from itertools import combinations, islice
from pathlib import Path
from typing import TypedDict

from dataset import ROOT, explicit_h_parser

CSD_URL = "https://github.com/jensengroup/xyz2mol_tm/blob/db2734b084045ef895528c7a3d5040d38da23983/SMILES_csvs/csd_smiles.csv"
CSD_SHA256 = "f9b5916ff14520a05e4c549f2140d6ae42cd92d6fc712037cf100c8eccdef991"
INDEX = "evidence/tables/excluded-records.csv"
INDEX_FIELDS = (
    "id",
    "exclusion",
    "reason",
    "input_formula",
    "basis",
    "sources",
    "xyz",
    "metadata",
)


def csd_rows(path, *, expected_sha256=CSD_SHA256):
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected_sha256:
        raise ValueError("Unexpected CSD transcription version")
    with path.open(encoding="utf-8-sig", newline="") as stream:
        yield from enumerate(csv.DictReader(stream), 2)


def element_counts(formula):
    return Counter(
        {s: int(n or 1) for s, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula)}
    )


def formula(counts):
    order = sorted(counts)
    if "C" in counts:
        order = [
            "C",
            *(["H"] if "H" in counts else []),
            *sorted(counts.keys() - {"C", "H"}),
        ]
    return "".join(s + (str(counts[s]) if counts[s] > 1 else "") for s in order)


class MetalBridgeExample(TypedDict, total=False):
    BHM_example: list[float]


class BoronGeometry(MetalBridgeExample):
    boron_atoms: list[int]
    nearest_BB: list[float]
    BB_contacts: list[list[float]]
    BBB_triangle: list[int]
    BHB_contacts: list[list[float]]


def boron_geometry(atoms) -> BoronGeometry:
    # All indices refer to the original XYZ, starting at 1. Cutoffs select
    # contacts for inspection, without assigning covalent bonds or bond orders.
    pos = {i: tuple(map(float, a[1:])) for i, a in enumerate(atoms, 1)}
    boron = [i for i, a in enumerate(atoms, 1) if a[0] == "B"]
    pairs = sorted((math.dist(pos[i], pos[j]), i, j) for i, j in combinations(boron, 2))
    contacts = [(i, j, round(d, 4)) for d, i, j in pairs if d <= 2.10]
    neighbors: dict[int, set[int]] = {i: set() for i in boron}
    for i, j, _ in contacts:
        neighbors[i].add(j)
        neighbors[j].add(i)
    triad = next(
        (
            [i, j, k]
            for i, j, _ in contacts
            for k in sorted(neighbors[i] & neighbors[j])
        ),
        [],
    )
    bridges = []
    for h, a in enumerate(atoms, 1):
        if a[0] == "H":
            nearby = [(b, math.dist(pos[h], pos[b])) for b in boron]
            for (i, d1), (j, d2) in combinations(
                [(b, d) for b, d in nearby if d <= 1.50], 2
            ):
                bridges.append([i, h, j, round(d1, 4), round(d2, 4)])
    d, i, j = pairs[0]
    return BoronGeometry(
        boron_atoms=boron,
        nearest_BB=[i, j, round(d, 4)],
        BB_contacts=[list(pair) for pair in contacts],
        BBB_triangle=triad,
        BHB_contacts=bridges,
    )


def boron_assessment(case, atoms):
    g = boron_geometry(atoms)
    i, j, d = g["nearest_BB"]
    basis = f"XYZ contains {len(g['boron_atoms'])} B atoms. The shortest contact is B{i}···B{j} at {d:.4f} Å."
    reason = "Multicentre boron framework outside the two-centre Lewis graph convention"
    if g["BBB_triangle"]:
        a, b, c = g["BBB_triangle"]
        basis += f" B{a}, B{b} and B{c} are pairwise within 2.10 Å. The structure has {len(g['BB_contacts'])} B···B contacts in this range."
    if g["BHB_contacts"]:
        first_b, bridge_h, second_b, d1, d2 = g["BHB_contacts"][0]
        basis += f" H{bridge_h} is {d1:.4f} Å from B{first_b} and {d2:.4f} Å from B{second_b}."
    if not (g["BBB_triangle"] or g["BHB_contacts"]):
        if case in ("NAHNOB", "QEJYAH"):
            metal = {"NAHNOB": "Sc", "QEJYAH": "Y"}[case]
            reason = "Multicentre bonding in a B–H–metal bridge"
            basis += f" The original same-ID CSD connection table shows a B–H–{metal} bridge."
            b, h, m = {"NAHNOB": (11, 12, 1), "QEJYAH": (6, 7, 1)}[case]
            assert [atoms[i - 1][0] for i in (b, h, m)] == ["B", "H", metal]
            pb, ph, pm = (tuple(map(float, atoms[i - 1][1:])) for i in (b, h, m))
            bh, mh = math.dist(pb, ph), math.dist(pm, ph)
            g["BHM_example"] = [b, h, m, round(bh, 4), round(mh, 4)]
            basis += (
                f" In XYZ, B{b}–H{h} is {bh:.4f} Å and H{h}···{metal}{m} is {mh:.4f} Å."
            )
        elif case == "FIGROC":
            basis += " The original same-ID CSD connection table shows a B/C cage coordinated to Ni."
        elif case == "YIBNIG":
            basis += " The same-ID table and XYZ support a six-vertex B₂C₄ cage. B2 lies 1.68–1.88 Å from basal B3, C4, C5, C6 and C7."
        else:
            raise ValueError(f"{case}: no reviewed evidence for this boron environment")
    return reason, basis, g


BRIDGE_AUDIT = "evidence/tables/boron-hydrogen-bridges.jsonl"


def load_bridge_audit():
    audit = {}
    with (ROOT / BRIDGE_AUDIT).open() as stream:
        for row in map(json.loads, stream):
            case = row["id"]
            if case in audit or row["evidence"] not in {
                "csd",
                "xyz_wbo",
                "geometry_only",
            }:
                raise ValueError(f"{case}: invalid or duplicate bridge audit entry")
            audit[case] = row
    return audit


def bridge_assessment(r):
    sites = r["xyz_sites"]
    if r["evidence"] == "csd":
        refs = [r["csd_row"]]
        kinds = Counter("B–H–" + s[3] for s in r["csd_sites"])
        basis = (
            "The original same-ID CSD graph records "
            + ", ".join(
                f"{n} {s} bridge" + ("s" if n != 1 else "") for s, n in kinds.items()
            )
            + "."
        )
        if element_counts(r["csd_component_formula"]) != element_counts(
            r["input_formula"]
        ):
            basis += f" The source component formula, {r['csd_component_formula']}, differs from the input."
        if not sites:
            basis += " The optimized XYZ falls outside the distance window. The correspondence of the bridge to the source structure remains uncertain."
    elif r["evidence"] == "xyz_wbo":
        refs = ["https://github.com/grimme-lab/g-xtb/releases/tag/v2.0.0"]
        sites = [s for s in sites if all(v is not None and v >= 0.10 for v in s[6:8])]
        basis = "The input XYZ and the g-xTB single point at the input electronic state support a boron hydrogen bridge."
    else:
        raise ValueError(f"{r['id']}: geometric proximity alone is insufficient")
    if sites:
        b, h, x, s, d1, d2, w1, w2 = sites[0]
        basis += (
            f" At B{b}···H{h}···{s}{x}, the two distances are {d1:.4f} and {d2:.4f} Å."
        )
        if w1 is not None and w2 is not None:
            basis += f" The corresponding WBOs are {w1:.3f} and {w2:.3f}."
    return {
        "reason": "A bridging hydrogen precludes a unique two-centre Lewis reference",
        "basis": basis,
        "input_formula": r["input_formula"],
        "sources": refs,
        "bridge_evidence": r["evidence"],
    }


def spin_reassessment(row, counts):
    """Summarize the alternative state and its evidence for a withdrawn graph."""
    review = row["exclusion_evidence"]["spin_reassessment"]
    kind = review["evidence_type"]
    candidate = review["candidate_state"]
    refs = []
    if kind == "published calculation":
        records = row["spin_evidence"]["published_spin_records_json"]
        multiplicity = int(row["spin_evidence"]["computed_multiplicity"])
        if not records or candidate != {3: "Triplet", 5: "Quintet"}[multiplicity]:
            raise ValueError(f"{row['id']}: inconsistent candidate spin")
        for match in records:
            entry = match["record"]
            if (
                entry["refcode"] != row["id"]
                or int(entry["total_charge"]) != row["charge"]
                or int(entry["spin_multiplicity"]) != multiplicity
                or match["input_xyz_sha256"] != row["xyz_sha256"]
                or not match["all_nearest_H_sites_match"]
            ):
                raise ValueError(f"{row['id']}: mismatched published spin record")
            refs.append("https://doi.org/" + match["data_doi"])
        datasets = ", ".join(match["source"] for match in records)
        basis = (
            f"The published multiplicity is {multiplicity} in {datasets} for the same CSD ID, "
            "charge and hydrogen sites at the source geometry."
        )
        for match in records:
            if "splitting" in match["record"]:
                basis += (
                    " The published vertical spin splitting is "
                    f"{match['record']['splitting']} kcal/mol."
                )
        if not any("splitting" in match["record"] for match in records):
            basis += " The supplementary record gives the multiplicity alone."
        if "measurement" in review:
            basis += " " + review["measurement"]
            refs.append(review["measurement_source"])
        reason = f"{candidate} reported for a singlet input"
    else:
        reason = "Alternative spin state remains unresolved"
        basis = review["basis"]
        if kind == "fixed-geometry calculation":
            fixed = row["case_evidence"]["fixed_xyz_review"][0]
            basis += " " + fixed["electronic_state_evidence"]
        refs = review["sources"]
    return {
        "reason": reason,
        "basis": f"Input multiplicity {row['multiplicity']}. {candidate}. {basis}",
        "input_formula": formula(counts),
        "sources": list(dict.fromkeys(refs)),
        "input_multiplicity": row["multiplicity"],
        "candidate_state": candidate,
        "evidence_type": kind,
    }


def assess(row, atoms, csd, bridges):
    case = row["id"]
    counts = Counter(a[0] for a in atoms)
    raw = row.get("exclusion_evidence") or {}
    refs = []
    source_row = csd.get(case)
    if source_row is not None:
        refs.append(f"{CSD_URL}#L{source_row[0]}")
    if row["exclusion"] == "Boron hydrogen bridge":
        return bridge_assessment(bridges[case])
    if row["exclusion"] == "Missing input hydrogens":
        if case == "EYINEF":
            source_formula = "C18H18HgN2O4"
            refs.extend(
                [
                    "https://doi.org/10.1107/S0108270104014131",
                    "https://www.crystallography.net/cod/2014185.cif",
                ]
            )
            source = "The published complex formula, excluding lattice water,"
        else:
            if source_row is None:
                raise ValueError(f"{case}: missing same-ID composition")
            source_formula = source_row[1]
            source = "The metal component in the same-ID CSD transcription"
        expected = element_counts(source_formula)
        delta = expected["H"] - counts["H"]
        if expected - counts != {"H": delta} or counts - expected or delta <= 0:
            raise ValueError(f"{case}: source and input disagree beyond missing H")
        reason = "Input XYZ lacks hydrogens present in the source component"
        basis = (
            f"{source} is {formula(expected)}. XYZ gives {formula(counts)} with the same heavy atoms and {delta} fewer H atom"
            + ("." if delta == 1 else "s.")
        )
        detail = {
            "source_formula": source_formula,
            "input_H": counts["H"],
            "expected_H": expected["H"],
            "missing_H": delta,
        }
    elif row["exclusion"] == "Questionable electronic state":
        if "spin_reassessment" in raw:
            return spin_reassessment(row, counts)
        reason = "Input singlet conflicts with the published spin assignment"
        state = raw.get("literature_multiplicity", raw.get("correct_multiplicity"))
        spin_basis = raw.get("multiplicity_basis") or raw.get("basis")
        if case == "BALFAY":
            spin_basis = "Variable-temperature EPR and the crystal structure support a triplet for the iron trinitrosyl complex of the same composition."
        elif case == "JETGIY":
            reason = "Input singlet uncertain for a composition reported as high-spin"
            state = "High-spin Fe(II). The multiplicity of this XYZ remains unresolved"
            spin_basis = "The report describes high-spin Fe(II) for this composition. A reconstruction from the published coordinates matches 29 heavy atoms with an RMSD of 0.393 Å and the positions of 14 listed H atoms. Independent same-ID crystal confirmation is unavailable."
        if not spin_basis or state is None:
            raise ValueError(f"{case}: missing spin evidence")
        basis = f"Input multiplicity {row['multiplicity']}. {spin_basis.rstrip('.')}."
        refs = [
            u
            for u in row["source_references"]
            if "doi.org/" in u and not u.endswith(")")
        ] + refs
        detail = {
            "input_multiplicity": row["multiplicity"],
            "reported_state": str(state),
        }
    elif row["exclusion"] == "Boron cluster":
        reason, basis, detail = boron_assessment(case, atoms)
    else:
        raise ValueError(f"{case}: unknown exclusion")
    refs.extend(
        u for u in row["source_references"] if u != CSD_URL and not u.endswith(")")
    )
    return dict(
        reason=reason,
        basis=basis,
        input_formula=formula(counts),
        sources=list(dict.fromkeys(refs)),
        **detail,
    )


def index_entry(row, part):
    a = row["exclusion_evidence"]["assessment"]
    return {
        "id": row["id"],
        "exclusion": row["exclusion"],
        "reason": a["reason"],
        "input_formula": a["input_formula"],
        "basis": a["basis"],
        "sources": " ".join(a["sources"]),
        "xyz": part["xyz"],
        "metadata": part["metadata"],
    }


def check_hydrogen_sites(case, atoms, sites, smiles_by_field):
    from rdkit import Chem

    positions = [tuple(map(float, atom[1:])) for atom in atoms]
    heavy = [i for i, atom in enumerate(atoms) if atom[0] != "H"]
    attached = Counter(
        min(heavy, key=lambda k: math.dist(positions[k], positions[i]))
        for i, atom in enumerate(atoms)
        if atom[0] == "H"
    )
    for site in sites:
        i = site["xyz_idx"]
        if atoms[i][0] != site["element"] or attached[i] != site["input_H"]:
            raise ValueError(f"{case}: input hydrogen count at XYZ site {i}")
    indexed = [site for site in sites if "source_idx" in site]
    if not indexed:
        return "XYZ only"
    for field, smiles in smiles_by_field.items():
        params = explicit_h_parser(sanitize=field != "smiles_csd_api")
        mol = Chem.MolFromSmiles(smiles, params) if smiles else None
        if mol is None:
            continue
        mol.UpdatePropertyCache(strict=False)
        for site in indexed:
            j = site["source_idx"]
            if j >= mol.GetNumAtoms():
                break
            atom = mol.GetAtomWithIdx(j)
            if (
                atom.GetSymbol() != site["element"]
                or atom.GetTotalNumHs(includeNeighbors=True) != site["source_H"]
            ):
                break
        else:
            return field
    raise ValueError(
        f"{case}: source indices or site hydrogen counts differ from the pinned CSD"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csd-csv", required=True, type=Path)
    parser.add_argument(
        "--check-sites",
        action="store_true",
        help="verify archived hydrogen-site counts against XYZ and CSD using RDKit",
    )
    args = parser.parse_args()
    csd, smiles = {}, {}
    for line, record in csd_rows(args.csd_csv):
        case = record["IDs"]
        csd[case] = (line, record["formula_heaviest_fragment"])
        if args.check_sites:
            smiles[case] = {
                field: record[field]
                for field in ("smiles_csd_api", "smiles_csd_api_fixed")
            }
    site_counts: Counter[str] = Counter()
    bridges = load_bridge_audit()
    for part in json.loads((ROOT / "manifest.json").read_bytes())["parts"]:
        with (
            (ROOT / part["metadata"]).open() as meta,
            (ROOT / part["xyz"]).open() as xyz,
        ):
            for row in map(json.loads, meta):
                n = int(next(xyz))
                if next(xyz).split()[0] != f"id={row['id']}":
                    raise ValueError("XYZ and metadata identifiers differ")
                atoms = list(islice(xyz, n))
                if not row["benchmark"]:
                    atoms = [s.split() for s in atoms]
                    a = assess(row, atoms, csd, bridges)
                    sites = row["exclusion_evidence"].get("site_mapping")
                    if (
                        args.check_sites
                        and row["exclusion"] == "Missing input hydrogens"
                        and sites
                    ):
                        field = check_hydrogen_sites(
                            row["id"], atoms, sites, smiles.get(row["id"], {})
                        )
                        site_counts[field] += 1
                    print(
                        json.dumps(
                            {"id": row["id"], "assessment": a}, ensure_ascii=False
                        )
                    )

    if args.check_sites:
        print(f"Hydrogen-site checks by source: {dict(site_counts)}", file=sys.stderr)


if __name__ == "__main__":
    main()
