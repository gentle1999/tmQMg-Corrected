"""Rebuild the two README overview plots and the seven paired molecular diagrams."""

import csv
import io
import json
from collections import Counter

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image
from rdkit import Chem
from rdkit.Chem import rdCoordGen, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from dataset import ROOT, TABLES, explicit_h_parser
from examples import draw_examples

TYPES = {
    "Boron cluster": "Boron clusters or hydrogen bridges",
    "Boron hydrogen bridge": "Boron clusters or hydrogen bridges",
    "Missing input hydrogens": "Atoms missing from XYZ",
    "Questionable electronic state": "Questionable electronic state",
}
CASES = (
    (
        "TIKBIC",
        "Boron cluster",
        "The graph detaches BH$_2^+$.\nH35 lies 1.278 and 1.468 Å\nfrom the two B atoms.",
    ),
    (
        "FOJFOZ",
        "Missing input hydrogens",
        "XYZ has 0 H. The source has 6 H\nin Hg–CH$_3$ and O–CH$_3$. The upstream\ngraph splits into four fragments.",
    ),
    (
        "BALFAY",
        "Questionable electronic state",
        "Singlet input, reported triplet.\nThe label shows Fe$^{4+}$ and three NO$^-$.\nFormal charges cannot establish spin.",
    ),
)
INK, MUTED, ACCENT = "#193744", "#455e69", "#08786e"
ORANGE = "#ab5527"


def upstream_image(case, smiles):
    params = explicit_h_parser(sanitize=False)
    mol = Chem.MolFromSmiles(smiles, params)
    mol.UpdatePropertyCache(strict=False)
    if case == "TIKBIC":
        # CoordGen lays out the boron framework on its undirected skeleton. The drawing
        # below retains every bond type and charge from the original SMILES.
        layout = Chem.Mol(mol)
        for bond in layout.GetBonds():
            if bond.GetBondType() == Chem.BondType.DATIVE:
                bond.SetBondType(Chem.BondType.SINGLE)
        rdCoordGen.AddCoords(layout)
        mol.AddConformer(layout.GetConformer())
    else:
        rdDepictor.Compute2DCoords(mol)
    highlight = {
        "TIKBIC": {"B"},
        "FOJFOZ": {"Hg"},
        "BALFAY": {"Fe"},
    }[case]
    sites = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in highlight]
    drawer = rdMolDraw2D.MolDraw2DCairo(760, 390)
    opts = drawer.drawOptions()
    opts.prepareMolsBeforeDrawing = False  # type: ignore[assignment]
    opts.setHighlightColour((0.98, 0.83, 0.69))
    opts.maxFontSize = 32  # type: ignore[assignment]
    opts.bondLineWidth = 2.5  # type: ignore[assignment]
    drawer.DrawMolecule(mol, highlightAtoms=sites, highlightBonds=[])
    drawer.FinishDrawing()
    return Image.open(io.BytesIO(drawer.GetDrawingText()))


def draw_input_overview():
    manifest = json.loads((ROOT / "manifest.json").read_bytes())
    with (ROOT / manifest["exclusion_index"]).open() as stream:
        index = {row["id"]: row for row in csv.DictReader(stream)}
    counts = Counter(TYPES[row["exclusion"]] for row in index.values())
    total, retained = manifest["inputs"], manifest["targets"]
    if retained + sum(counts.values()) != total:
        raise ValueError("Exclusion counts do not partition the input collection")

    fig = plt.figure(figsize=(12, 5.8), facecolor="white")

    def label(x, y, message, size: float = 12, color=INK, **kw):
        fig.text(x, y, message, fontsize=size, color=color, va="center", **kw)

    label(0.025, 0.946, "Input collection", 17, weight="bold")
    label(0.407, 0.946, "Excluded structures", 17, weight="bold")
    label(0.665, 0.946, "Upstream SMILES", 12, MUTED)
    label(0.025, 0.822, f"{total:,} inputs", 21, weight="bold")
    bar = fig.add_axes((0.025, 0.716, 0.325, 0.052))
    bar.barh([0], [retained], color=ACCENT)
    bar.barh([0], [len(index)], left=retained, color=ORANGE)
    bar.set(xlim=(0, total), ylim=(-0.4, 0.4))
    bar.set_axis_off()
    label(
        0.025,
        0.670,
        f"Released graphs   {retained:,}   ({retained / total:.2%})",
        13,
        ACCENT,
    )
    label(
        0.025,
        0.608,
        f"Excluded     {len(index):,}   ({len(index) / total:.2%})",
        13,
        ORANGE,
    )
    label(0.025, 0.480, "Primary reason for exclusion", 14, weight="bold")
    for i, (_, reason, _) in enumerate(CASES):
        name = TYPES[reason]
        n = counts[name]
        y = 0.380 - i * 0.121
        label(0.025, y, name)
        label(0.349, y, f"{n:,}", ha="right", weight="bold")
        ax = fig.add_axes((0.025, y - 0.054, 0.324, 0.020))
        ax.barh([0], [n], color=ORANGE)
        ax.set(xlim=(0, max(counts.values())), ylim=(-0.4, 0.4))
        ax.set_axis_off()

    short_names = {
        "Boron cluster": "B–H–B bridges",
        "Missing input hydrogens": "Missing H atoms",
        "Questionable electronic state": "Spin assignment",
    }
    for i, (case, reason, caption) in enumerate(CASES):
        row = index[case]
        if row["exclusion"] != reason:
            raise ValueError(f"{case}: unexpected primary exclusion")
        with (ROOT / row["metadata"]).open() as stream:
            record = next(r for line in stream if (r := json.loads(line))["id"] == case)
        y = 0.742 - i * 0.285
        ax = fig.add_axes((0.383, y - 0.127, 0.285, 0.250))
        ax.imshow(upstream_image(case, record["original_smiles"]))
        ax.set_axis_off()
        label(
            0.677, y + 0.071, f"{case}   ·   {short_names[reason]}", 13, weight="bold"
        )
        label(0.677, y - 0.032, caption, 11.5, MUTED, linespacing=1.50)
        if i < 2:
            fig.add_artist(
                Line2D(
                    [0.391, 0.977],
                    [y - 0.144] * 2,
                    transform=fig.transFigure,
                    color="#dce4e6",
                    lw=0.8,
                )
            )
    output = ROOT / "evidence/figures/input-overview.png"
    fig.savefig(output, dpi=180)
    plt.close(fig)


LABELS = {
    "bond_order": "Intraligand bond order",
    "connectivity": "Ligand or counterion connectivity",
    "metal_charge": "Anomalous metal valence",
    "zero_valence": "Incorrect zero-valent metal",
    "electronic_state": "Charge or electron count inconsistent with input",
    "composition": "Atomic composition inconsistent with XYZ",
    "alternating_charges": "Alternating charges in a conjugated system",
}


def draw_issue_distribution():
    with (TABLES / "equivalence-rescore.csv").open() as stream:
        retained = {row["id"] for row in csv.DictReader(stream)}
    if len(retained) != json.loads((ROOT / "manifest.json").read_bytes())["targets"]:
        raise ValueError("Retained ID count differs from the release manifest")

    hits: dict[str, set[str]] = {key: set() for key in LABELS}
    with (TABLES / "upstream-issue-types.csv").open() as stream:
        rows = csv.DictReader(stream)
        if set(rows.fieldnames or ()) != {"id", *LABELS}:
            raise ValueError("Unexpected issue categories")
        for row in rows:
            if row["id"] in retained:
                for key, ids in hits.items():
                    if row[key]:
                        ids.add(row["id"])

    keys = sorted(hits, key=lambda key: len(hits[key]), reverse=True)
    counts = [len(hits[key]) for key in keys]
    affected = len(set.union(*hits.values()))
    fig, ax = plt.subplots(figsize=(11.8, 4.7))
    fig.subplots_adjust(left=0.42, right=0.90, top=0.78, bottom=0.10)
    fig.text(
        0.025,
        0.89,
        "Upstream graph annotation issues",
        size=18,
        weight="bold",
        color="#193744",
    )
    fig.text(
        0.025,
        0.81,
        f"{affected:,} of {len(retained):,} retained complexes meet at least one criterion. Categories overlap.",
        size=11,
        color="#455e69",
    )
    bars = ax.barh(range(len(keys)), counts, height=0.53, color="#08786e")
    ax.set_yticks(range(len(keys)), [LABELS[key] for key in keys], fontsize=11)
    ax.invert_yaxis()
    ax.set_xlim(left=0)
    ax.xaxis.set_major_formatter("{x:,.0f}")
    ax.bar_label(bars, labels=[f"{n:,}" for n in counts], padding=5, fontsize=11)
    ax.set_axisbelow(True)
    ax.grid(axis="x", color="#dee5e6", linewidth=0.7)
    ax.spines[:].set_visible(False)
    ax.tick_params(length=0, colors="#193744")
    fig.savefig(ROOT / "evidence/figures/upstream-issue-distribution.png", dpi=180)
    plt.close(fig)


def main():
    (ROOT / "evidence/figures").mkdir(exist_ok=True)
    draw_input_overview()
    draw_issue_distribution()
    draw_examples()


if __name__ == "__main__":
    main()
