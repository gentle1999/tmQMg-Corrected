"""Compare the metal-free ligand graphs of all recomputed nonmatching complexes."""

from collections import Counter, defaultdict

from rdkit import Chem

from dataset import TM

# Compare the loss and gain of two written forms on the same metal-free skeleton
# and the same set of symmetry classes. An optional key restricts the whole fragment.
FORMS = {
    "carbonyl_from_anionic_oxygen": (None, "[#6]-[O-]", "[#6]=[O+0]"),
    "imine_from_anionic_nitrogen": (None, "[#6]-[N-]", "[#6]=[N+0]"),
    "thiocarbonyl_to_thiolate": (None, "[#6]=[S+0]", "[#6]-[S-]"),
    "ncn_charge_separation": (None, "[N+0]-[C+0]-[N+0]", "[N+]=[C-]-[N+0]"),
    "borohydride_from_bh3": ("B", "[#5H3+0]", "[#5H4-]"),
    "borohydride_from_bh2": ("B", "[#5H2+]", "[#5H4-]"),
    "dioxygen": ("OO", "[O-]-[O-]", "[O+0]=[O+0]"),
    "nitrosyl": ("NO", "[N-]=[O+0]", "[N+0]#[O+]"),
    "dithiolene": (None, "[S+0]=[#6]-[#6]=[S+0]", "[S-]-[#6]=[#6]-[S-]"),
    "azopyridine": (
        "C1CCC(NNC2CCCCN2)CC1",
        "[#6]-[N-]-[#7]=[#6]",
        "[#6]-[N+0]=[N+0]-[#6]",
    ),
    "dimethylglyoximate": (
        "CC(NO)C(C)NO",
        "[O+0]=[#7]-[#6]=[#6]-[#7]=[O+0]",
        "[O-]-[#7]=[#6]-[#6]=[#7]-[O-]",
    ),
    "nitrate": ("ON(O)O", "[O+0]-[N+0](-[O+0])-[O+0]", "[O-]-[N+](-[O-])=[O+0]"),
    "ortho_quinonediimine": (
        "NC1CCCCC1N",
        "[N-]-[#6]~[#6]-[N-]",
        "[N+0]=[#6]~[#6]=[N+0]",
    ),
    "acetonitrile": ("CCN", "[#6]-[#6]-[#7]", "[#6]-[#6]#[#7]"),
    "dimethyl_sulfoxide": (
        "CS(C)O",
        "[#6]-[S+0](-[#6])=[O+0]",
        "[#6]-[S+](-[#6])-[O-]",
    ),
    "azide": ("NNN", "[N-2]-[N+]#[N+0]", "[N-]=[N+]=[N-]"),
}
QUERIES = {
    name: (key, Chem.MolFromSmarts(a), Chem.MolFromSmarts(b))
    for name, (key, a, b) in FORMS.items()
}


def heavy_skeleton(mol, *, include_metals=False):
    graph = Chem.RWMol()
    indices = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1 or (
            not include_metals and atom.GetAtomicNum() in TM
        ):
            continue
        node = Chem.Atom(atom.GetAtomicNum())
        node.SetIsotope(atom.GetIsotope())
        node.SetNoImplicit(True)
        indices[atom.GetIdx()] = graph.AddAtom(node)
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if i in indices and j in indices:
            graph.AddBond(indices[i], indices[j], Chem.BondType.SINGLE)
    graph.UpdatePropertyCache(strict=False)
    return graph, list(indices)


def profiles(mol):
    graph, indices = heavy_skeleton(mol)
    members: list[tuple[int, ...]] = []
    fragments = Chem.GetMolFrags(
        graph, asMols=True, sanitizeFrags=False, fragsMolAtomMapping=members
    )
    groups = defaultdict(list)
    centers: Counter[tuple[str, tuple[str, ...]]] = Counter()
    locations = {}
    for fragment, positions in zip(fragments, members, strict=True):
        key = Chem.MolToSmiles(fragment)
        ranks = list(Chem.CanonicalRankAtoms(fragment, breakTies=False))
        Chem.GetSymmSSSR(fragment)
        props: defaultdict[tuple[str, str], Counter[int | str]] = defaultdict(Counter)
        ids = [indices[i] for i in positions]
        charge = 0
        for i, atom in enumerate(fragment.GetAtoms()):
            orig = mol.GetAtomWithIdx(ids[i])
            site = f"{atom.GetSymbol()}{ranks[i]}"
            props["atom_charge", site][orig.GetFormalCharge()] += 1
            props["attached_H", site][orig.GetTotalNumHs(includeNeighbors=True)] += 1
            locations[ids[i]] = (key, ranks[i])
            charge += orig.GetFormalCharge()
            neighbors = tuple(sorted(a.GetSymbol() for a in atom.GetNeighbors()))
            centers[atom.GetSymbol(), neighbors] += 1
        props["net_charge", ""][charge] = 1
        for bond in fragment.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            raw = mol.GetBondBetweenAtoms(ids[i], ids[j])
            pair = sorted(
                f"{fragment.GetAtomWithIdx(k).GetSymbol()}{ranks[k]}" for k in (i, j)
            )
            kind = "ring_order" if bond.IsInRing() else "bond_order"
            props[kind, "-".join(pair)][str(raw.GetBondType())] += 1
        groups[key].append((props, [i + 1 for i in ids]))
    forms: defaultdict[tuple[str, str, tuple[int, ...]], Counter[int]] = defaultdict(
        Counter
    )
    for name, (required, source, target) in QUERIES.items():
        if required is not None and required not in groups:
            continue
        for side, query in enumerate((source, target)):
            for match in mol.GetSubstructMatches(query, uniquify=True, maxMatches=0):
                key = locations[match[0]][0]
                if required is None or key == required:
                    site_ranks = tuple(sorted(locations[idx][1] for idx in match))
                    forms[name, key, site_ranks][side] += 1
    return groups, centers, forms


def delta(before, after):
    return {
        "before": dict(sorted((before - after).items())),
        "after": dict(sorted((after - before).items())),
    }


def compare_groups(left, right):
    # Cancel identical profiles before pooling equivalent copies or recording unmatched fragments.
    residual = list(right)
    remaining = []
    for item in left:
        match = next(
            (i for i, other in enumerate(residual) if item[0] == other[0]), None
        )
        if match is None:
            remaining.append(item)
        else:
            residual.pop(match)
    return remaining, residual


def pooled(items):
    result: defaultdict[tuple[str, str], Counter[int | str]] = defaultdict(Counter)
    for props, _ in items:
        for key, values in props.items():
            result[key].update(values)
    return result


def fragment_labels(mol):
    # Keep the H count on each original atom when the metal and H nodes disappear.
    labels = []
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol().lower() if atom.GetIsAromatic() else atom.GetSymbol()
        count = atom.GetTotalNumHs(includeNeighbors=True)
        hydrogen = {0: "", 1: "H"}.get(count, f"H{count}")
        charge = atom.GetFormalCharge()
        sign = {-1: "-", 0: "", 1: "+"}.get(charge, f"{charge:+d}")
        labels.append(f"[{atom.GetIsotope() or ''}{symbol}{hydrogen}{sign}]")
    return labels


def written_fragments(mol, labels, groups):
    return [
        Chem.MolFragmentToSmiles(
            mol,
            atomsToUse=[i - 1 for i in ids],
            atomSymbols=labels,
            isomericSmiles=False,
        )
        for _, ids in groups
    ]


def fragment_audit(old, new, group):
    # Keep the H count assigned in the complete complex before removing the metal.
    def states(mol, copies):
        values = []
        for ids in copies:
            atoms = [mol.GetAtomWithIdx(i - 1) for i in ids]
            values.append(
                (
                    sum(a.GetFormalCharge() for a in atoms),
                    sum(a.GetTotalNumHs(includeNeighbors=True) for a in atoms),
                )
            )
        return sorted(values)

    before = states(old, group["upstream_atoms"])
    after = states(new, group["sdf_atoms"])
    result = {"upstream_charge_H": before, "sdf_charge_H": after}
    if [q for q, _ in before] != [q for q, _ in after]:
        return {**result, "status": "charge_changed"}
    if before != after:
        return {**result, "status": "hydrogen_changed"}
    # With more than one distinct residual form we do not guess a correspondence.
    left, right = set(group["upstream_forms"]), set(group["sdf_forms"])
    status = "undetermined"
    if len(left) == len(right) == 1:
        a, b = (Chem.MolFromSmiles(next(iter(forms))) for forms in (left, right))
        if a is not None and b is not None:
            a_key, b_key = Chem.MolToSmiles(a), Chem.MolToSmiles(b)
            if a_key == b_key:
                status = "normalized_identity"
            else:
                supplier = Chem.ResonanceMolSupplier(
                    a, Chem.ALLOW_CHARGE_SEPARATION, 100000
                )
                if any(
                    m is not None and Chem.MolToSmiles(m) == b_key for m in supplier
                ):
                    status = "resonance_witness"
    return {**result, "status": status}


def classify(old, new):
    def inventory(mol):
        return Counter(
            (a.GetAtomicNum(), a.GetIsotope())
            for a in mol.GetAtoms()
            if a.GetAtomicNum() != 1 and a.GetAtomicNum() not in TM
        )

    if inventory(old) != inventory(new):
        raise ValueError(
            "Ligand connectivity comparison requires conserved heavy atoms"
        )
    left, left_centers, old_forms = profiles(old)
    right, right_centers, new_forms = profiles(new)
    matched = []
    unmatched = []
    unpaired_charges: list[defaultdict[str, Counter[int]]] = [
        defaultdict(Counter),
        defaultdict(Counter),
    ]
    old_labels, new_labels = fragment_labels(old), fragment_labels(new)
    for skeleton in sorted(left.keys() | right.keys()):
        a, b = compare_groups(left[skeleton], right[skeleton])
        if not a and not b:
            continue
        info = {
            "skeleton": skeleton,
            "upstream_atoms": [ids for _, ids in a],
            "sdf_atoms": [ids for _, ids in b],
            "upstream_forms": written_fragments(old, old_labels, a),
            "sdf_forms": written_fragments(new, new_labels, b),
        }
        if len(a) != len(b):
            for items, totals, molecule, field in zip(
                (a, b),
                unpaired_charges,
                (old, new),
                ("upstream_charges", "sdf_charges"),
                strict=True,
            ):
                fragment_charges = []
                for _, ids in items:
                    atoms = [molecule.GetAtomWithIdx(i - 1) for i in ids]
                    fragment_charges.append([atom.GetFormalCharge() for atom in atoms])
                    for atom in atoms:
                        totals[atom.GetSymbol()][atom.GetFormalCharge()] += 1
                info[field] = fragment_charges
            unmatched.append(info)
            continue
        first, second = pooled(a), pooled(b)
        changes = []
        for kind, site in sorted(first):
            x, y = first[kind, site], second[kind, site]
            if x == y:
                continue
            # Written aromatic vs localized ring orders do not establish a chemical edit.
            change_kind = kind
            if kind in {"bond_order", "ring_order"}:
                labels = x.keys() | y.keys()
                if "AROMATIC" in labels:
                    change_kind = "aromatic_encoding"
                elif not labels <= {"SINGLE", "DOUBLE", "TRIPLE", "QUADRUPLE"}:
                    change_kind = "other_bond_encoding"
            changes.append({"kind": change_kind, "site": site, **delta(x, y)})
        if changes:
            matched.append({**info, "changes": changes})
    centers = []
    if unmatched:
        for symbol in sorted(
            {key[0] for key in left_centers.keys() | right_centers.keys()}
        ):
            a = Counter(
                {",".join(k): v for (s, k), v in left_centers.items() if s == symbol}
            )
            b = Counter(
                {",".join(k): v for (s, k), v in right_centers.items() if s == symbol}
            )
            if a != b:
                centers.append({"element": symbol, **delta(a, b)})
    families = set()
    audit_families = defaultdict(set)
    for site in old_forms.keys() & new_forms.keys():
        name, key, _ = site
        a, b = old_forms[site], new_forms[site]
        if len(left[key]) == len(right[key]) and a[0] > b[0] and a[1] < b[1]:
            families.add(name)
            if name in {"imine_from_anionic_nitrogen", "carbonyl_from_anionic_oxygen"}:
                audit_families[key].add(name)
    for group in matched:
        if group["skeleton"] in audit_families:
            group["CN_CO_audit"] = {
                "families": sorted(audit_families.pop(group["skeleton"])),
                **fragment_audit(old, new, group),
            }
    if audit_families:
        raise ValueError("Missing matched ligand for a C–N/C–O audit")
    transitions = {r["element"]: (set(r["before"]), set(r["after"])) for r in centers}
    for name, element, before, after, required in (
        ("perchlorate", "Cl", "O", "O,O,O,O", "OCl(O)(O)O"),
        ("selenoether", "Se", "C", "C,C", None),
        ("carbon_fluorine", "F", "", "C", None),
        ("tetrafluoroborate", "B", "F,F,F", "F,F,F,F", "FB(F)(F)F"),
        ("hexafluorophosphate", "P", "F,F,F,F,F", "F,F,F,F,F,F", "FP(F)(F)(F)(F)F"),
    ):
        if transitions.get(element) == ({before}, {after}) and (
            required is None or len(right[required]) > len(left[required])
        ):
            families.add(name)
    x, y = unpaired_charges
    charge_changes = {
        s: delta(x[s], y[s]) for s in sorted(x.keys() | y.keys()) if x[s] != y[s]
    }
    kinds = {change["kind"] for group in matched for change in group["changes"]}
    categories = []
    if charge_changes or "atom_charge" in kinds:
        categories.append("atom_charge")
    if kinds & {"bond_order", "ring_order"}:
        categories.append("bond_order")
    if unmatched:
        categories.append("connectivity")
    return {
        "families": sorted(families),
        "matched_skeletons": matched,
        "unmatched_skeletons": unmatched,
        "changed_neighbors": centers,
        "unpaired_atom_charge_deltas": charge_changes,
        "categories": categories or ["no_confirmed_change_in_three_categories"],
    }
