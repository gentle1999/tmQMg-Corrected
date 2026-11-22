"""Write explicit-H SMILES in the reference XYZ atom-number space."""

from collections import defaultdict

from rdkit import Chem

from dataset import TM, explicit_h_parser
from ligand_edits import heavy_skeleton


def fragments(mol, *, include_metals):
    skeleton, indices = heavy_skeleton(mol, include_metals=include_metals)
    sites: list[tuple[int, ...]] = []
    pieces = Chem.GetMolFrags(
        skeleton, asMols=True, sanitizeFrags=False, fragsMolAtomMapping=sites
    )
    return [
        (piece, [indices[i] for i in group])
        for piece, group in zip(pieces, sites, strict=True)
    ]


def match_fragments(left, right, mapping, *, exact):
    # Largest fragments first. Already assigned atoms are excluded from later matches.
    occupied = set(mapping.values())
    for piece, sites in sorted(left, key=lambda item: -len(item[1])):
        if any(i in mapping for i in sites):
            continue
        for target, candidates in right:
            if len(sites) > len(candidates) or all(j in occupied for j in candidates):
                continue
            if exact and (len(sites), piece.GetNumBonds()) != (
                len(candidates),
                target.GetNumBonds(),
            ):
                continue
            opts = Chem.SubstructMatchParameters()
            if occupied.intersection(candidates):
                opts.setExtraFinalCheck(
                    lambda _, match, sites=candidates: all(
                        sites[k] not in occupied for k in match
                    )
                )
            match = target.GetSubstructMatch(piece, opts)
            if match:
                paired = [candidates[k] for k in match]
                mapping.update(zip(sites, paired, strict=True))
                occupied.update(paired)
                break


def topology_mapping(upstream, reference):
    mapping: dict[int, int] = {}
    match_fragments(
        fragments(upstream, include_metals=True),
        fragments(reference, include_metals=True),
        mapping,
        exact=True,
    )
    left = fragments(upstream, include_metals=False)
    right = fragments(reference, include_metals=False)
    match_fragments(left, right, mapping, exact=True)
    match_fragments(left, right, mapping, exact=False)
    reverse = {j: i for i, j in mapping.items()}
    match_fragments(right, left, reverse, exact=False)
    mapping = {i: j for j, i in reverse.items()}

    for metal in {a.GetAtomicNum() for a in upstream.GetAtoms()} & TM:
        old = [a.GetIdx() for a in upstream.GetAtoms() if a.GetAtomicNum() == metal]
        new = [a.GetIdx() for a in reference.GetAtoms() if a.GetAtomicNum() == metal]
        if len(old) == len(new) == 1:
            mapping[old[0]] = new[0]

    def hydrogens(mol, parents):
        groups = defaultdict(list)
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 1:
                neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
                if all(i in parents for i in neighbors):
                    key = (
                        atom.GetIsotope(),
                        tuple(sorted(parents[i] for i in neighbors)),
                    )
                    groups[key].append(atom.GetIdx())
        return groups

    old_h = hydrogens(upstream, mapping)
    new_h = hydrogens(reference, {i: i for i in mapping.values()})
    for key, sites in old_h.items():
        mapping.update(zip(sites, new_h[key], strict=False))
    return mapping


def numbered_smiles(mol, numbers):
    for atom, number in zip(mol.GetAtoms(), numbers, strict=True):
        atom.SetAtomMapNum(number)
    try:
        return Chem.MolToSmiles(mol, canonical=False, allHsExplicit=True)
    except RuntimeError as error:
        if "bad bond stereo" not in str(error):
            raise
        # Some native E/Z labels lack stereo-defining atoms. Reconstruct their
        # SMILES directions from the neighboring bonds of the released graph.
        for bond in mol.GetBonds():
            if (
                bond.GetStereo() in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ)
                and not bond.GetStereoAtoms()
            ):
                bond.SetStereo(Chem.BondStereo.STEREONONE)
        Chem.SetBondStereoFromDirections(mol)
        return Chem.MolToSmiles(mol, canonical=False, allHsExplicit=True)


def mapped_pair(smiles, reference, witness):
    count = reference.GetNumAtoms()
    right = numbered_smiles(reference, range(1, count + 1))
    if not smiles:
        return "", right

    params = explicit_h_parser(sanitize=True)
    source = Chem.MolFromSmiles(smiles, params)
    if source is None:
        params = explicit_h_parser(sanitize=False)
        source = Chem.MolFromSmiles(smiles, params)
        if source is None:
            return "", right
        source.UpdatePropertyCache(strict=False)
    # Parsing for the historic score can collapse written H atoms and change the
    # atom order. In that case, align the H-preserving graph independently.
    written_h = any(a.GetAtomicNum() == 1 for a in source.GetAtoms())
    upstream = Chem.AddHs(source)
    if witness and not written_h and len(witness) == upstream.GetNumAtoms():
        mapping = dict(enumerate(witness))
    else:
        mapping = topology_mapping(upstream, reference)
    if len(set(mapping.values())) != len(mapping):
        raise ValueError("Non-injective atom correspondence")
    extra = count
    numbers = []
    for atom in upstream.GetAtoms():
        i = atom.GetIdx()
        if i in mapping:
            other = reference.GetAtomWithIdx(mapping[i])
            if (atom.GetAtomicNum(), atom.GetIsotope()) != (
                other.GetAtomicNum(),
                other.GetIsotope(),
            ):
                raise ValueError("Atom correspondence changes element or isotope")
            numbers.append(mapping[i] + 1)
        else:
            extra += 1
            numbers.append(extra)
    return numbered_smiles(upstream, numbers), right
