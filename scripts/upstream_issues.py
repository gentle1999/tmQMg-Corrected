"""Audit upstream labels against input constraints and supported reference differences."""

import csv
import json
from collections import Counter
from functools import lru_cache
from itertools import combinations
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import NamedTuple

from molgr.fallback.utils.consts import (
    METAL_VALENCE_AVAILABLE_MINOR,
    METAL_VALENCE_AVAILABLE_PRIOR,
)
from rdkit import Chem

from dataset import ROOT, TABLES, TM, explicit_h_parser
from ligand_edits import classify, heavy_skeleton

PARAMS = explicit_h_parser(sanitize=False)
ELEMENTS = Chem.GetPeriodicTable()


def xyz_records(stream):
    for count in stream:
        fields = dict(item.split("=", 1) for item in next(stream).split())
        atoms = Counter(next(stream).split()[0] for _ in range(int(count)))
        yield fields, atoms


def formula(atoms):
    return " ".join(f"{element}{n}" for element, n in sorted(atoms.items()) if n)


def electrons(composition, charge):
    return sum(ELEMENTS.GetAtomicNumber(s) * n for s, n in composition.items()) - charge


class MetalState(NamedTuple):
    formal_charge: int
    parent_h: int
    terminal_neutral_h: int
    terminal_hydride: int
    adjusted_state: int


def metal_state(atom):
    """Reassign neutral terminal M–H to H− while preserving the pair charge."""
    charge = atom.GetFormalCharge()
    parent_h = atom.GetTotalNumHs(includeNeighbors=False)
    terminal = Counter(
        h.GetFormalCharge()
        for h in atom.GetNeighbors()
        if h.GetAtomicNum() == 1 and h.GetDegree() == 1
    )
    if terminal.keys() - {0, -1}:
        raise ValueError(f"Unsupported terminal H charge on {atom.GetSymbol()}")
    return MetalState(
        charge, parent_h, terminal[0], terminal[-1], charge + parent_h + terminal[0]
    )


def has_only_donor_bonds(metal):
    return metal.GetDegree() > 0 and all(
        bond.GetBondType() == Chem.BondType.DATIVE
        and bond.GetEndAtomIdx() == metal.GetIdx()
        for bond in metal.GetBonds()
    )


def zero_valence_evidence(record, old, new, before, after):
    # Covalent M–ligand bonds can give a neutral formal charge even for an
    # oxidized metal. Compare ionic ligand assignments in the dative convention.
    if not (
        before.formal_charge == before.adjusted_state == 0
        and after.formal_charge == after.adjusted_state != 0
    ):
        return ""
    source = record["graph_source"]
    if record["tier"] not in (1, 2) or source not in (
        "Primary report",
        "Original CSD connection table for the same ID",
        "Matches the revised CSD connectivity for the same ID",
    ):
        return ""
    if not all(has_only_donor_bonds(metal) for metal in (old, new)):
        return ""
    charge = record["charge"]
    if any(
        Chem.GetFormalCharge(metal.GetOwningMol()) != charge for metal in (old, new)
    ):
        return ""
    return (
        f"{old.GetSymbol()} 0 → {after.adjusted_state:+d}, input charge {charge:+d}, "
        f"nonmetal charge sum {charge:+d} → {charge - after.formal_charge:+d}, "
        f"{source} (tier {record['tier']})"
    )


def audit_metal_pair(record, upstream, reference, writer):
    centers = [
        [a for a in mol.GetAtoms() if a.GetAtomicNum() in TM]
        for mol in (upstream, reference)
    ]
    if any(len(atoms) != 1 for atoms in centers):
        return "", ""
    (old,), (new,) = centers
    element = old.GetSymbol()
    if element != new.GetSymbol():
        return "", ""
    allowed = set(METAL_VALENCE_AVAILABLE_PRIOR[element]) | set(
        METAL_VALENCE_AVAILABLE_MINOR[element]
    )
    allowed_text = "|".join(map(str, sorted(allowed)))
    before, after = metal_state(old), metal_state(new)
    hit = before.adjusted_state not in allowed and after.adjusted_state in allowed
    zero = zero_valence_evidence(record, old, new, before, after)
    writer.writerow(
        (record["id"], element, allowed_text, *before, *after, int(hit), zero)
    )
    abnormal = (
        f"{element} {before.adjusted_state}/{after.adjusted_state} (MolGR: {allowed_text})"
        if hit
        else ""
    )
    return abnormal, zero


def label_constraints(mol, record, xyz_atoms):
    atoms = Counter(a.GetSymbol() for a in mol.GetAtoms())
    atoms["H"] += sum(a.GetTotalNumHs() for a in mol.GetAtoms())
    charge = Chem.GetFormalCharge(mol)
    input_charge = record["charge"]
    old_e = electrons(atoms, charge)
    xyz_e = electrons(xyz_atoms, input_charge)
    unpaired = record["multiplicity"] - 1
    state = []
    if charge != input_charge:
        state.append(f"charge {charge}/{input_charge}")
    if old_e != xyz_e:
        state.append(f"electrons {old_e}/{xyz_e}")
    if old_e < unpaired or (old_e - unpaired) % 2:
        state.append(f"electrons {old_e}, input multiplicity {unpaired + 1}")
    composition = (
        "" if +atoms == xyz_atoms else f"{formula(atoms)} / {formula(xyz_atoms)}"
    )
    return "; ".join(state), composition


def path_key(ranks, i, j, k):
    return (min(ranks[i], ranks[k]), ranks[j], max(ranks[i], ranks[k]))


def alternating_sites(mol, ranks):
    sites = set()
    for atom in mol.GetAtoms():
        opposite = [
            a.GetIdx()
            for a in atom.GetNeighbors()
            if a.GetFormalCharge() * atom.GetFormalCharge() < 0
        ]
        for i, k in combinations(opposite, 2):
            sites.add(path_key(ranks, i, atom.GetIdx(), k))
    return sites


DIENE = Chem.MolFromSmarts("[*]=[*]-[*]=[*]")
AROMATIC = Chem.MolFromSmarts("[a]:[a]:[a]")


def conjugated_sites(mol, ranks):
    sites = set()
    for query in (DIENE, AROMATIC):
        for match in mol.GetSubstructMatches(query, uniquify=False, maxMatches=0):
            for i, j, k in zip(match, match[1:], match[2:], strict=False):
                if all(mol.GetAtomWithIdx(a).GetFormalCharge() == 0 for a in (i, j, k)):
                    sites.add(path_key(ranks, i, j, k))
    return sites


@lru_cache(maxsize=4096)
def fragment(smiles):
    mol = Chem.MolFromSmiles(smiles, PARAMS)
    mol.UpdatePropertyCache(strict=False)
    return mol


def ranks(mol):
    graph, _ = heavy_skeleton(mol)
    return list(Chem.CanonicalRankAtoms(graph, breakTies=False))


def alternating_charge_hits(edit):
    alternation = []
    for group in edit["matched_skeletons"]:
        before, after = set(group["upstream_forms"]), set(group["sdf_forms"])
        # Discard ambiguous pairings between different copies of the same skeleton.
        if len(before) != 1 or len(after) != 1:
            continue
        left, right = before.pop(), after.pop()
        a, b = fragment(left), fragment(right)
        if Chem.GetFormalCharge(a) == Chem.GetFormalCharge(b) and sum(
            x.GetTotalNumHs() for x in a.GetAtoms()
        ) == sum(x.GetTotalNumHs() for x in b.GetAtoms()):
            continue
        if "+" in left and "-" in left:
            sites = alternating_sites(a, ranks(a))
            if sites:
                right_ranks = ranks(b)
                if (sites - alternating_sites(b, right_ranks)) & conjugated_sites(
                    b, right_ranks
                ):
                    alternation.append(group["skeleton"])
    return "|".join(alternation)


@lru_cache(maxsize=8192)
def charge_and_h(form):
    mol = Chem.MolFromSmiles(form, sanitize=False)
    mol.UpdatePropertyCache(strict=False)
    return Chem.GetFormalCharge(mol), sum(a.GetTotalNumHs() for a in mol.GetAtoms())


def add_ligand_issue_types(row, edit):
    connections = set()
    if edit["unmatched_skeletons"]:
        connections.add("heavy_atom_skeleton")
    bonds = set()
    for group in edit["matched_skeletons"]:
        kinds = {change["kind"] for change in group["changes"]}
        if "attached_H" in kinds and not row["composition"]:
            connections.add("hydrogen_site")
        if kinds & {"bond_order", "ring_order"}:
            old, new = (
                Counter(map(charge_and_h, group[k]))
                for k in ("upstream_forms", "sdf_forms")
            )
            if old != new:
                bonds.add(group["skeleton"])
    for key, hits in (
        ("connectivity", connections),
        ("bond_order", bonds),
    ):
        if hits:
            row[key] = "|".join(sorted(hits))


ISSUE_TYPES: tuple[str, ...] = (
    "electronic_state",
    "composition",
    "connectivity",
    "bond_order",
    "alternating_charges",
    "metal_charge",
    "zero_valence",
)


def write_issues(output):
    manifest = json.loads((ROOT / "manifest.json").read_bytes())
    with (TABLES / "equivalence-rescore.csv").open() as stream:
        wanted = {
            r["id"]
            for r in csv.DictReader(stream)
            if r["recomputed"] == "Not equivalent"
        }
    classified = []
    coverage: Counter[str] = Counter()
    with (
        (output / "metal-valence-audit.csv").open("w", newline="") as audit,
        (output / "upstream-issue-types.csv").open("w", newline="") as types,
        (output / "ligand-edits.jsonl").open("w") as edits,
    ):
        writer = csv.writer(audit, lineterminator="\n")
        type_writer = csv.DictWriter(types, ("id", *ISSUE_TYPES), lineterminator="\n")
        type_writer.writeheader()
        state_columns = [
            f"{side}_{field}"
            for side in ("upstream", "reference")
            for field in MetalState._fields
        ]
        writer.writerow(
            ("id", "metal", "allowed_states", *state_columns, "hit", "zero_valence")
        )
        for part in manifest["parts"]:
            with (
                (ROOT / part["metadata"]).open() as meta,
                (ROOT / part["xyz"]).open() as xyz,
                (ROOT / part["sdf"]).open("rb") as sdf,
            ):
                references = iter(
                    Chem.ForwardSDMolSupplier(sdf, sanitize=False, removeHs=False)
                )
                for line, (fields, atoms) in zip(meta, xyz_records(xyz), strict=True):
                    record = json.loads(line)
                    case = record["id"]
                    assert (case, record["charge"], record["multiplicity"]) == (
                        fields["id"],
                        int(fields["charge"]),
                        int(fields["multiplicity"]),
                    )
                    reference = next(references) if record["benchmark"] else None
                    if record["benchmark"] and (
                        reference is None or reference.GetProp("CASE_ID") != case
                    ):
                        raise ValueError(f"Missing or misaligned SDF for {case}")
                    if not record["original_smiles"]:
                        coverage["missing"] += 1
                        continue
                    old = Chem.MolFromSmiles(record["original_smiles"], PARAMS)
                    if old is None:
                        coverage["unreadable"] += 1
                        continue
                    coverage["readable"] += 1
                    old.UpdatePropertyCache(strict=False)
                    state, composition = label_constraints(old, record, atoms)
                    row = {"electronic_state": state, "composition": composition}
                    if reference is not None:
                        reference.UpdatePropertyCache(strict=False)
                        row["metal_charge"], row["zero_valence"] = audit_metal_pair(
                            record, old, reference, writer
                        )
                    if case in wanted:
                        if reference is None:
                            raise ValueError(
                                f"Missing reference for nonmatching case {case}"
                            )
                        edit = {"id": case, **classify(old, reference)}
                        classified.append(case)
                        print(
                            json.dumps(edit, ensure_ascii=False, separators=(",", ":")),
                            file=edits,
                        )
                        row["alternating_charges"] = alternating_charge_hits(edit)
                        add_ligand_issue_types(row, edit)
                    if any(row.values()):
                        type_writer.writerow({"id": case, **row})
                if next(references, False) is not False:
                    raise ValueError(f"Extra SDF record in {part['sdf']}")

    if classified != sorted(wanted):
        raise ValueError(
            "Expected one ordered ligand comparison for every nonmatching case"
        )
    print(dict(coverage))


def main():
    with TemporaryDirectory(dir=TABLES) as directory:
        output = Path(directory)
        write_issues(output)
        for table in output.iterdir():
            table.replace(TABLES / table.name)


if __name__ == "__main__":
    main()
