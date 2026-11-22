"""Recompute boron exclusion evidence from XYZ and the pinned CSD CSV."""

import argparse
import csv
import json
from collections import Counter
from itertools import combinations, islice
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from rdkit import Chem

from dataset import ROOT, TABLES, TM, explicit_h_parser
from exclusions import BRIDGE_AUDIT, CSD_URL, csd_rows, element_counts, formula

PT = Chem.GetPeriodicTable()
DIBORON_RING = Chem.MolFromSmarts("[#5]1~[#6]~[#6]~[#5]~[#6]~[#6]~1")
BHB_CLASS = "Framework with B–H–B bridges"
BHM_CLASS = "B–N–B ligands with hydrogen bridges to the metal and no B–B bonds"


def source_bridges(sources):
    params = explicit_h_parser(sanitize=False)
    found = {}
    for case, (line, smiles, component_formula) in sources.items():
        if "B" not in smiles or "H" not in smiles or "[" not in smiles:
            continue
        mol = Chem.MolFromSmiles(smiles, params)
        if mol is None:
            continue
        sites = []
        for h in mol.GetAtoms():
            if h.GetAtomicNum() != 1:
                continue
            for first, second in combinations(h.GetNeighbors(), 2):
                a, b = (
                    (second, first) if second.GetAtomicNum() == 5 else (first, second)
                )
                if a.GetAtomicNum() == 5 and (
                    b.GetAtomicNum() == 5 or b.GetAtomicNum() in TM
                ):
                    sites.append(
                        [a.GetIdx() + 1, h.GetIdx() + 1, b.GetIdx() + 1, b.GetSymbol()]
                    )
        if sites:
            found[case] = {
                "csd_row": f"{CSD_URL}#L{line}",
                "csd_sites": sites,
                "csd_component_formula": component_formula,
            }
    return found


def geometry(row, atoms):
    symbols = [a[0] for a in atoms]
    xyz = np.array([a[1:] for a in atoms], dtype=float)
    borons = [i for i, s in enumerate(symbols) if s == "B"]
    metals = [i for i, s in enumerate(symbols) if PT.GetAtomicNumber(s) in TM]
    runs = [
        g
        for g in row["gxtb"]
        if g["status"] == "converged"
        and g["charge"] == row["charge"]
        and g["uhf"] == row["multiplicity"] - 1
    ]
    analysis = (
        row.get("wbo_analysis")
        or row.get("previous_graph_review", {}).get("wbo_analysis")
        or {}
    )
    cohort = analysis.get("cohort")
    run = next((g for g in runs if g["cohort"] == cohort), runs[0] if runs else None)
    wbo = {
        tuple(sorted((i, j))): v
        for i, j, v in (run.get("wbo_pairs", []) if run else [])
    }
    sites = []
    for i, s in enumerate(symbols):
        if s != "H":
            continue
        distances = np.linalg.norm(xyz - xyz[i], axis=1)
        near = [j for j in borons if distances[j] <= 1.75]
        for j in near:
            for k in metals + [b for b in near if b > j]:
                limit = (
                    1.75 if symbols[k] == "B" else PT.GetRcovalent(symbols[k]) + 0.95
                )
                if distances[k] <= limit:
                    sites.append(
                        [
                            j + 1,
                            i + 1,
                            k + 1,
                            symbols[k],
                            round(distances[j], 4),
                            round(distances[k], 4),
                            wbo.get(tuple(sorted((j + 1, i + 1)))),
                            wbo.get(tuple(sorted((k + 1, i + 1)))),
                        ]
                    )
    return sites, run.get("cohort") if run else None


def supported(site):
    return all(v is not None and v >= 0.10 for v in site[6:8])


def scan(csd):
    for part in json.loads((ROOT / "manifest.json").read_bytes())["parts"]:
        with (
            (ROOT / part["xyz"]).open() as xyz,
            (ROOT / part["metadata"]).open("rb") as meta,
        ):
            for raw in meta:
                n = int(next(xyz))
                case = next(xyz).split()[0][3:]
                lines = list(islice(xyz, n))
                if case not in csd and not any(
                    line.split(None, 1)[0] == "B" for line in lines
                ):
                    continue
                row = json.loads(raw)
                if row["id"] != case:
                    raise ValueError(case)
                atoms = [line.split() for line in lines]
                sites, cohort = geometry(row, atoms)
                if sites or case in csd:
                    source = csd.get(case, {})
                    evidence = (
                        "csd"
                        if source
                        else "xyz_wbo"
                        if any(map(supported, sites))
                        else "geometry_only"
                    )
                    yield dict(
                        id=case,
                        evidence=evidence,
                        input_formula=formula(Counter(a[0] for a in atoms)),
                        **source,
                        xyz_sites=sites,
                        wbo_cohort=cohort,
                    )


def graph_evidence(smiles):
    params = explicit_h_parser(sanitize=False)
    mol = Chem.MolFromSmiles(smiles, params)
    if mol is None:
        return None
    z = [a.GetAtomicNum() for a in mol.GetAtoms()]
    neighbors = [{a.GetIdx() for a in atom.GetNeighbors()} for atom in mol.GetAtoms()]
    boron = {i for i, element in enumerate(z) if element == 5}
    # Triangles include at least one B, and all other vertices are B or C.
    faces = {
        tuple(sorted((b, i, j)))
        for b in boron
        for i, j in combinations(neighbors[b], 2)
        if z[i] in (5, 6) and z[j] in (5, 6) and j in neighbors[i]
    }
    edge_use = Counter(edge for face in faces for edge in combinations(face, 2))
    bhb = bhm = bnb = 0
    for i, adjacent in enumerate(neighbors):
        n_b = len(boron & adjacent)
        if z[i] == 1:
            bhb += n_b >= 2
            bhm += n_b >= 1 and any(z[a] in TM for a in adjacent)
        elif z[i] == 7:
            bnb += n_b >= 2
    # Broad topology screen for EDAKUP-like rings. No aromaticity or charge assumption.
    diboron_rings = 0
    for ring in mol.GetSubstructMatches(DIBORON_RING):
        members = set(ring)
        if all(len(neighbors[i] & members) == 2 for i in ring) and all(
            not (neighbors[i] & boron) for i in members & boron
        ):
            diboron_rings += 1
    return {
        "para_C4B2_rings_without_BB": diboron_rings,
        "borons": len(boron),
        "BB_edges": sum(len(neighbors[i] & boron) for i in boron) // 2,
        "BBB_triangles": sum(all(i in boron for i in face) for face in faces),
        "BC_faces": len(faces),
        "shared_face_edges": sum(count >= 2 for count in edge_use.values()),
        "BHB_hydrogens": bhb,
        "BHM_hydrogens": bhm,
        "BNB_nitrogens": bnb,
    }


def classify(raw, xyz):
    if raw is None:
        if xyz["BBB_triangle"]:
            return "Triangular B contacts without an original same-ID connection table"
        return "No original same-ID connection table available for verification"
    if raw["shared_face_edges"]:
        return "Boron framework with edge-sharing triangular faces"
    if raw["BHB_hydrogens"]:
        return BHB_CLASS
    if not raw["BB_edges"] and raw["BNB_nitrogens"] and raw["BHM_hydrogens"]:
        return BHM_CLASS
    return "Same-ID graph requires review against other criteria"


def hydrogen_inventory(smiles, *, sanitize):
    params = explicit_h_parser(sanitize=sanitize)
    mol = Chem.MolFromSmiles(smiles, params)
    mol.UpdatePropertyCache(strict=sanitize)
    atoms = list(Chem.AddHs(mol).GetAtoms())
    borons = [a for a in atoms if a.GetAtomicNum() == 5]
    (metal,) = [a for a in atoms if a.GetAtomicNum() in TM]

    def count_h(atom):
        return sum(h.GetAtomicNum() == 1 for h in atom.GetNeighbors())

    return {
        "formula": formula(Counter(a.GetSymbol() for a in atoms)),
        "H_per_B": sorted(count_h(b) for b in borons),
        "metal_H": count_h(metal),
        "metal_formal_charge": metal.GetFormalCharge(),
        "total_charge": Chem.GetFormalCharge(mol),
    }


def bridge_review(rec, csd_smiles, component_formula, upstream, raw, xyz):
    original = rec["original_smiles"]
    label = hydrogen_inventory(original, sanitize=True)
    source = hydrogen_inventory(csd_smiles, sanitize=False)
    components = original.split(".")
    return {
        "xyz_formula": xyz["input_formula"],
        "upstream_formula": label["formula"],
        "csd_component_formula": formula(element_counts(component_formula)),
        "input_charge": rec["charge"],
        "upstream_charge": label["total_charge"],
        "isolated_H_plus": components.count("[H+]"),
        "isolated_BH2_plus": components.count("[BH2+]"),
        "upstream_BB_edges": upstream["BB_edges"],
        "csd_BB_edges": raw["BB_edges"],
        "upstream_BHB": upstream["BHB_hydrogens"],
        "csd_BHB": raw["BHB_hydrogens"],
        "upstream_BHM": upstream["BHM_hydrogens"],
        "csd_BHM": raw["BHM_hydrogens"],
        "upstream_H_per_B": label["H_per_B"],
        "csd_H_per_B": source["H_per_B"],
        "upstream_metal_H": label["metal_H"],
        "upstream_metal_formal_charge": label["metal_formal_charge"],
        "xyz_BHB_contacts": xyz["BHB_contacts"],
        "xyz_BHM_example": xyz.get("BHM_example", []),
        "upstream_smiles": original,
    }


def write_cluster_reviews(sources, output):
    rows = []
    for part in json.loads((ROOT / "manifest.json").read_bytes())["parts"]:
        with (ROOT / part["metadata"]).open() as stream:
            for rec in map(json.loads, stream):
                if rec["exclusion"] != "Boron cluster":
                    continue
                case = rec["id"]
                xyz = rec["exclusion_evidence"]["assessment"]
                line, smiles, component_formula = sources.get(case, (None, "", ""))
                raw = (
                    graph_evidence(smiles)
                    if smiles and smiles != "API_smiles_missing"
                    else None
                )
                upstream = (
                    graph_evidence(rec["original_smiles"])
                    if rec["original_smiles"]
                    else None
                )
                n_b = len(xyz["boron_atoms"])
                if raw is not None and raw["borons"] != n_b:
                    raise ValueError(f"{case}: different B counts in CSD and input")
                degree = Counter(i for a, b, _ in xyz["BB_contacts"] for i in (a, b))
                category = classify(raw, xyz)
                source_url = f"{CSD_URL}#L{line}" if line else ""
                review = ""
                if category in (BHB_CLASS, BHM_CLASS):
                    detail = bridge_review(
                        rec, smiles, component_formula, upstream, raw, xyz
                    )
                    review = json.dumps(detail, separators=(",", ":"))
                rows.append(
                    {
                        "id": case,
                        "classification": category,
                        "B_atoms": n_b,
                        "xyz_BB_contacts": len(xyz["BB_contacts"]),
                        "xyz_max_BB_neighbors": max(degree.values(), default=0),
                        "xyz_BBB_triangle": int(bool(xyz["BBB_triangle"])),
                        "xyz_BHB_contacts": len(xyz["BHB_contacts"]),
                        "csd_raw_status": (
                            "readable"
                            if raw is not None
                            else "no_row"
                            if line is None
                            else "missing"
                            if smiles in ("", "API_smiles_missing")
                            else "unreadable"
                        ),
                        "csd_row": source_url,
                        "csd_graph": json.dumps(raw, separators=(",", ":"))
                        if raw is not None
                        else "",
                        "upstream_graph": json.dumps(upstream, separators=(",", ":"))
                        if upstream is not None
                        else "",
                        "bridge_review": review,
                    }
                )
    with (output / "boron-exclusion-audit.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(dict(Counter(r["classification"] for r in rows)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csd-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=TABLES)
    args = parser.parse_args()
    sources = {
        row["IDs"]: (line, row["smiles_csd_api"], row["formula_heaviest_fragment"])
        for line, row in csd_rows(args.csd_csv)
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=args.output_dir) as directory:
        output = Path(directory)
        with (output / Path(BRIDGE_AUDIT).name).open("w") as stream:
            for record in scan(source_bridges(sources)):
                print(
                    json.dumps(record, ensure_ascii=False, separators=(",", ":")),
                    file=stream,
                )
        write_cluster_reviews(sources, output)
        for table in output.iterdir():
            table.replace(args.output_dir / table.name)


if __name__ == "__main__":
    main()
