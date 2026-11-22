"""Render paired Lewis graphs on the released 3D conformer with rdkit-dof."""

import csv
import io
import json
import xml.etree.ElementTree as ET

import networkx as nx
import numpy as np
from rdkit import Chem
from rdkit_dof import DofDrawSettings, MolsToDofSvgAnimation, MolToDofImage

from dataset import ROOT, explicit_h_parser

EXAMPLES = {
    "ABEVAG": (
        "Graph and XYZ hydrogen counts",
        "5 × CH₃–CH₂–NH · 30 H in total",
        "5 × CH₃–C≡N · 15 H, matching XYZ",
        (34, -113, -175),
    ),
    "DEFWIY": (
        "C–C bonds in the two alkynyl ligands",
        "2 × C–C · donor carbons −3",
        "2 × C≡C · donor carbons −1",
        (0, 25, 4),
    ),
    "SAZFEF": (
        "Ligand bond orders and Ni formal charge",
        "Two C=S bonds, central C–C, Ni(0)",
        "Two C–S⁻ bonds, central C=C, Ni(+2)",
        (-113, 180, 180),
    ),
    "RAWVAN": (
        "Connectivity of the perchlorate counterion",
        "One Cl–O bond · Ag(+7)",
        "Four Cl–O bonds · Ag(+1)",
        (-50, 20, 0),
    ),
    "JEGLOA": (
        "Methyl groups bonded to selenium",
        "One methyl group detached from Se",
        "Both methyl groups bonded to Se",
        (-133, 60, 143),
    ),
    "ERULUB": (
        "Two oxygen labels shift the total charge by two",
        "Two S=O bonds · charge +2 · 300 electrons",
        "Two S–O⁻ bonds · charge 0 · 302 electrons",
        (-116, 10, 0),
    ),
    "CUDBOU": (
        "Alternating charges around a six-membered carbon ring",
        "C10⁻, C12⁺, C14⁻, C16⁺",
        "Four neutral C atoms · alternating bonds",
        (-178, 3, 0),
    ),
}
# One-based XYZ indices for the labeled insets.
DETAIL_SITES = {
    "CUDBOU": ("Ring detail", (7, 8, 10, 12, 14, 16)),
    "ERULUB": ("S and the two O atoms", (2, 3, 4)),
    "DEFWIY": ("Au–alkynyl detail", (14, 13, 1, 38, 39)),
}
SETTINGS = DofDrawSettings(
    env_file=None, min_alpha=0.58, atom_colors={1: (0.30, 0.35, 0.39)}
)
SVG = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG)


def element(name, **attrs):
    return ET.Element(f"{{{SVG}}}{name}", {k: str(v) for k, v in attrs.items()})


def label(x, y, value, size, weight="normal", color="#1f3645"):
    obj = element("text", x=x, y=y, fill=color)
    obj.set(
        "font-family",
        "DejaVu Sans, Arial, sans-serif",
    )
    obj.set("font-size", str(size))
    obj.set("font-weight", weight)
    obj.text = value
    return obj


def skeleton(mol):
    graph: nx.Graph[int] = nx.Graph()
    for a in mol.GetAtoms():
        graph.add_node(a.GetIdx(), element=a.GetAtomicNum())
    for b in mol.GetBonds():
        graph.add_edge(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
    return graph


def read_examples():
    cases = {}
    for part in json.loads((ROOT / "manifest.json").read_text())["parts"]:
        if not any(part["first"] <= case <= part["last"] for case in EXAMPLES):
            continue
        with (ROOT / part["metadata"]).open("rb") as stream:
            records = {
                d["id"]: d
                for line in stream
                if (d := json.loads(line))["id"] in EXAMPLES
            }
        if not records:
            continue
        with (ROOT / part["sdf"]).open("rb") as stream:
            for record in stream:
                name = record.strip().decode()
                if name in records:
                    block = [record]
                    for line in stream:
                        block.append(line)
                        if line == b"$$$$\n":
                            break
                    mol = next(
                        Chem.ForwardSDMolSupplier(
                            io.BytesIO(b"".join(block)), sanitize=False, removeHs=False
                        )
                    )
                    cases[name] = records[name], mol
    if cases.keys() != EXAMPLES.keys():
        raise ValueError("Missing example")
    return cases


def order(mol, i, j):
    bond = mol.GetBondBetweenAtoms(i, j)
    return bond.GetBondTypeAsDouble() if bond else 0


def highlight_site(case, pair, before, after):
    return {
        "ABEVAG": pair == ["C", "N"] and (before, after) == (1, 3),
        "DEFWIY": pair == ["C", "C"] and (before, after) == (1, 3),
        "SAZFEF": pair in (["C", "S"], ["C", "C"]) and {before, after} == {1, 2},
        "JEGLOA": pair == ["C", "Se"] and (before, after) == (0, 1),
        "RAWVAN": pair == ["Cl", "O"],
        "ERULUB": pair == ["O", "S"] and (before, after) == (2, 1),
        "CUDBOU": False,  # Select the charged ring atoms below.
    }[case]


def align(case, old, new):
    before, after = skeleton(old), skeleton(new)
    if case == "JEGLOA":
        after.remove_edge(
            *next(
                (i, j)
                for i, j in after.edges
                if {
                    new.GetAtomWithIdx(i).GetSymbol(),
                    new.GetAtomWithIdx(j).GetSymbol(),
                }
                == {"C", "Se"}
            )
        )
    if case == "RAWVAN":
        # Upstream [O-]Cl retains one terminal oxygen of unknown XYZ identity.
        # Use the lowest-index terminal O as one reproducible correspondence.
        cl = next(a.GetIdx() for a in new.GetAtoms() if a.GetSymbol() == "Cl")
        neighbors = list(after[cl])
        retained = min(i for i in neighbors if after.degree(i) == 1)
        after.remove_edges_from((cl, i) for i in neighbors if i != retained)
    matcher = nx.isomorphism.GraphMatcher(after, before, node_match=lambda a, b: a == b)
    mapping = min(
        matcher.isomorphisms_iter(),
        key=lambda m: sum(
            order(old, m[i], m[j]) != order(new, i, j) for i, j in after.edges
        ),
    )
    positions = new.GetConformer()
    if not positions.Is3D():
        raise ValueError(f"{case}: expected released 3D coordinates")
    conf = Chem.Conformer(old.GetNumAtoms())
    conf.Set3D(True)
    for i, j in mapping.items():
        conf.SetAtomPosition(j, positions.GetAtomPosition(i))
    old.RemoveAllConformers()
    old.AddConformer(conf)
    return mapping


def hydrogen_neighbors(atom):
    return [a.GetIdx() for a in atom.GetNeighbors() if a.GetAtomicNum() == 1]


def upstream_on_xyz(old, new, target, mapping):
    # Materialize only hydrogen sites present in the source XYZ. Excess SMILES
    # valence stays on the parent atom and never receives a position.
    hydrogens_by_parent = {}
    for i, j in mapping.items():
        xyz_idx = new.GetAtomWithIdx(i).GetIntProp("xyz_index") - 1
        actual_h = hydrogen_neighbors(target.GetAtomWithIdx(xyz_idx))
        if actual_h:
            hydrogens_by_parent[j] = actual_h
    mol = Chem.RWMol(
        Chem.AddHs(old, onlyOnAtoms=list(hydrogens_by_parent), addCoords=False)
    )
    conf = mol.GetConformer()
    xyz = target.GetConformer()
    excess = []
    for j, xyz_neighbors in hydrogens_by_parent.items():
        added_h = hydrogen_neighbors(mol.GetAtomWithIdx(j))
        if len(added_h) < len(xyz_neighbors):
            raise ValueError(f"Ambiguous H correspondence at SMILES atom {j + 1}")
        for idx, xyz_idx in zip(added_h, xyz_neighbors, strict=False):
            conf.SetAtomPosition(idx, xyz.GetAtomPosition(xyz_idx))
        excess.extend(added_h[len(xyz_neighbors) :])
    for idx in sorted(excess, reverse=True):
        mol.RemoveAtom(idx)
    mol.UpdatePropertyCache(strict=False)
    for atom in mol.GetAtoms():
        count = atom.GetNumImplicitHs()
        if count:
            subscript = f"<sub>{count}</sub>" if count > 1 else ""
            atom.SetProp("atomLabel", f"{atom.GetSymbol()}H{subscript}")
    if mol.GetNumAtoms() != target.GetNumAtoms():
        raise ValueError("Each displayed atom must have an XYZ position")
    return mol


def viewing_frames(mol, target, angles):
    coordinates = mol.GetConformer().GetPositions()
    heavy = [a.GetIdx() for a in target.GetAtoms() if a.GetAtomicNum() != 1]
    center = target.GetConformer().GetPositions()[heavy].mean(axis=0)
    x, y, z = (coordinates - center).T
    yaw, pitch, roll = np.deg2rad(angles)
    x, z = np.cos(yaw) * x + np.sin(yaw) * z, -np.sin(yaw) * x + np.cos(yaw) * z
    y, z = np.cos(pitch) * y - np.sin(pitch) * z, np.sin(pitch) * y + np.cos(pitch) * z
    x, y = np.cos(roll) * x - np.sin(roll) * y, np.sin(roll) * x + np.cos(roll) * y
    frames = []
    for phase in np.linspace(0, 2 * np.pi, 24, endpoint=False):
        angle = np.deg2rad(22) * np.sin(phase)
        positions = np.column_stack(
            (
                x * np.cos(angle) + z * np.sin(angle),
                y,
                -x * np.sin(angle) + z * np.cos(angle),
            )
        )
        frames.append(positioned_copy(mol, positions))
    return frames


def positioned_copy(mol, positions):
    frame = Chem.Mol(mol)
    for idx, point in enumerate(positions):
        frame.GetConformer().SetAtomPosition(idx, point)
    return frame


def site_detail(frame, indices, color, labels=None, width=235):
    # Preserve isolated atoms, including disconnected upstream oxygens.
    detail = Chem.RWMol(frame)
    for index, name in (labels or {}).items():
        atom = detail.GetAtomWithIdx(index)
        charge = atom.GetFormalCharge()
        if charge:
            magnitude = str(abs(charge)) if abs(charge) > 1 else ""
            sign = "+" if charge > 0 else "-"
            atom.SetProp("atomLabel", f"{name}<sup>{magnitude}{sign}</sup>")
        else:
            atom.SetProp("atomLabel", name)
    for idx in reversed(range(frame.GetNumAtoms())):
        if idx not in indices:
            detail.RemoveAtom(idx)
    for atom in detail.GetAtoms():
        atom.SetNoImplicit(True)
    detail.UpdatePropertyCache(strict=False)
    return ET.fromstring(
        MolToDofImage(
            detail,
            size=(width, 125),
            return_image=False,
            settings=SETTINGS,
            highlightAtoms=list(range(detail.GetNumAtoms())),
            highlightBonds=list(range(detail.GetNumBonds())),
            highlightColor=color,
            prepareMolsBeforeDrawing=False,
            clearBackground=False,
            maxFontSize=19,
        )
    )


def draw(case, metadata, target):
    for atom in target.GetAtoms():
        atom.SetIntProp("xyz_index", atom.GetIdx() + 1)
    target.UpdatePropertyCache(strict=False)
    h = Chem.RemoveHsParameters()
    h.updateExplicitCount = True  # type: ignore[assignment]
    new = Chem.RemoveHs(target, h, sanitize=False)
    params = explicit_h_parser(sanitize=False)
    old = Chem.MolFromSmiles(metadata["original_smiles"], params)
    old.UpdatePropertyCache(strict=False)
    old = Chem.RemoveHs(old, h, sanitize=False)
    for mol in (old, new):
        Chem.RemoveStereochemistry(mol)
        mol.UpdatePropertyCache(strict=False)
    mapping = align(case, old, new)
    upstream = upstream_on_xyz(old, new, target, mapping)
    implicit_h = sum(a.GetNumImplicitHs() for a in upstream.GetAtoms())
    if implicit_h != (15 if case == "ABEVAG" else 0):
        raise ValueError(f"{case}: unexpected upstream implicit H count {implicit_h}")
    rows = []
    atoms: list[set[int]] = [set(), set()]
    bonds: list[list[int]] = [[], []]
    wbo = {tuple(sorted(row[:2])): row for row in metadata["wbo_analysis"]["bonds"]}
    for bond in new.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        u, v = mapping[i], mapping[j]
        pair = sorted(new.GetAtomWithIdx(k).GetSymbol() for k in (i, j))
        before, after = order(old, u, v), order(new, i, j)
        if not highlight_site(case, pair, before, after):
            continue
        atoms[0].update([u, v])
        a, b = sorted(new.GetAtomWithIdx(k).GetIntProp("xyz_index") for k in (i, j))
        atoms[1].update([a - 1, b - 1])
        if before:
            bonds[0].append(old.GetBondBetweenAtoms(u, v).GetIdx())
        bonds[1].append(target.GetBondBetweenAtoms(a - 1, b - 1).GetIdx())
        row = wbo[a, b]
        rows.append(
            [
                case,
                a,
                b,
                row[2],
                int(before),
                int(after),
                row[5],
                row[7],
            ]
        )
    expected_rows = {
        "ABEVAG": 5,
        "DEFWIY": 2,
        "SAZFEF": 3,
        "JEGLOA": 1,
        "RAWVAN": 4,
        "ERULUB": 2,
        "CUDBOU": 0,
    }
    assert len(rows) == expected_rows[case], (case, rows)
    for i, j in mapping.items():
        if (
            new.GetAtomWithIdx(i).GetFormalCharge()
            != old.GetAtomWithIdx(j).GetFormalCharge()
        ):
            atoms[0].add(j)
            atoms[1].add(new.GetAtomWithIdx(i).GetIntProp("xyz_index") - 1)
    panels = []
    details = []
    for k, mol in enumerate((upstream, target)):
        color = (0.68, 0.27, 0.07, 1) if k == 0 else (0.015, 0.44, 0.39, 1)
        frames = viewing_frames(mol, target, EXAMPLES[case][3])
        if case == "RAWVAN":
            cl = next(a for a in new.GetAtoms() if a.GetSymbol() == "Cl")
            site = [cl, *cl.GetNeighbors()]
            indices = {
                mapping[a.GetIdx()] if k == 0 else a.GetIntProp("xyz_index") - 1
                for a in site
            }
            details.append(site_detail(frames[0], indices, color))
        if case in DETAIL_SITES:
            _, xyz_indices = DETAIL_SITES[case]
            site_labels = {}
            for atom in new.GetAtoms():
                xyz_id = atom.GetIntProp("xyz_index")
                if xyz_id in xyz_indices:
                    idx = mapping[atom.GetIdx()] if k == 0 else xyz_id - 1
                    site_labels[idx] = (
                        atom.GetSymbol()
                        if case == "DEFWIY"
                        else f"{atom.GetSymbol()}{xyz_id}"
                    )
            width = 325 if case == "DEFWIY" else 235
            details.append(
                site_detail(frames[0], site_labels.keys(), color, site_labels, width)
            )
        svg = MolsToDofSvgAnimation(
            frames,
            duration=[850] + [110] * (len(frames) - 1),
            size=(460, 295),
            return_image=False,
            settings=SETTINGS,
            highlightAtomLists=[sorted(atoms[k])] * len(frames),
            highlightBondLists=[bonds[k]] * len(frames),
            highlightColor=color,
            padding=0.05,
            maxFontSize=17,
            bondLineWidth=2.2,
            flagCloseContactsDist=0,
            clearBackground=False,
            prepareMolsBeforeDrawing=False,
        )
        panels.append(ET.fromstring(svg))
    return panels, sorted(rows), details


def draw_examples():
    cases = read_examples()
    figures = ROOT / "evidence" / "figures"
    tables = ROOT / "evidence" / "tables"
    figures.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)
    all_rows = []
    for case, (title, left, right, _) in EXAMPLES.items():
        metadata, mol = cases[case]
        panels, changes, details = draw(case, metadata, mol)
        all_rows.extend(changes)
        height = 528 if details else 398
        canvas = element("svg", width=960, height=height, viewBox=f"0 0 960 {height}")
        canvas.append(element("rect", width=960, height=height, fill="white"))
        canvas.append(element("rect", x=0, y=0, width=960, height=40, fill="#eaf0f2"))
        canvas.append(
            label(
                16,
                27,
                f"{case}    {title}    ·    Evidence tier {metadata['tier']}",
                20,
                "bold",
            )
        )
        for panel, x, name, caption, color in zip(
            panels,
            (10, 490),
            ("Upstream tmQMg", "tmQMg-Corrected"),
            (left, right),
            ("#914218", "#007362"),
            strict=True,
        ):
            canvas.append(label(x + 12, 66, name, 18, "bold", color))
            panel.set("x", str(x))
            panel.set("y", "70")
            canvas.append(panel)
            text = label(x + 230, height - 15, caption, 18, "bold", color)
            text.set("text-anchor", "middle")
            canvas.append(text)
        for x, detail in zip((10, 490), details, strict=False):
            caption = (
                DETAIL_SITES[case][0] if case in DETAIL_SITES else "Cl and four O atoms"
            )
            canvas.append(label(x + 22, 387, caption, 15))
            detail.set("x", str(x + (122 if case == "DEFWIY" else 155)))
            detail.set("y", "360")
            canvas.append(detail)
        ET.ElementTree(canvas).write(
            figures / f"{case.lower()}-3d.svg", encoding="utf-8", xml_declaration=True
        )
    with (tables / "example-bonds.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(
            [
                "id",
                "xyz_i",
                "xyz_j",
                "element_pair",
                "upstream_order",
                "released_order",
                "distance_A",
                "gxtb_WBO",
            ]
        )
        writer.writerows(all_rows)
