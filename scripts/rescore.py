"""Recompare released upstream SMILES and native reference graphs with MolGR."""

import argparse
import base64
import csv
import hashlib
import json
import multiprocessing as mp
import sys
import zlib
from importlib.metadata import version
from pathlib import Path

from molgr.utils.equivalence import evaluate_equivalence
from rdkit import Chem

from dataset import ROOT
from exclusions import csd_rows
from mapped_smiles import mapped_pair

LABELS = {
    "equivalent": "Equivalent",
    "not_equivalent": "Not equivalent",
    "inconclusive": "Inconclusive",
}


def decode_reference(encoded: str) -> Chem.Mol:
    binary = zlib.decompress(base64.b64decode(encoded, validate=True))
    # RDKit accepts its native binary serialization, but the stub declares str.
    return Chem.Mol(binary)  # type: ignore[call-overload]


def compare(smiles, encoded, options):
    if not smiles:
        return "Missing label", "Missing upstream SMILES", ()
    upstream = Chem.MolFromSmiles(smiles)
    if upstream is None:
        return "Unreadable label", "RDKit could not parse upstream SMILES", ()
    reference = decode_reference(encoded)
    result = evaluate_equivalence(Chem.AddHs(upstream), reference, **options)
    return LABELS[result.decision.value], result.reason, result.atom_mapping


def worker(pipe, options):
    with pipe:
        pipe.send("ready")
        while True:
            smiles, encoded = pipe.recv()
            pipe.send(compare(smiles, encoded, options))


def comparisons(records, policy, *, include_smiles=False):
    # Native C++ can keep running past a Python alarm. A timed-out worker
    # is terminated and the following molecule gets a fresh process.
    options = dict(policy)
    timeout = options.pop("comparison_timeout_seconds")
    context = mp.get_context("spawn")
    process = None
    pipe = None
    try:
        for case, published, smiles, encoded in records:
            if process is None or pipe is None:
                pipe, remote = context.Pipe()
                process = context.Process(target=worker, args=(remote, options))
                process.start()
                remote.close()
                pipe.recv()  # Import time is outside the per-comparison budget.
            pipe.send((smiles, encoded))
            if pipe.poll(timeout):
                try:
                    label, reason, witness = pipe.recv()
                except EOFError as error:
                    raise RuntimeError(
                        f"{case}: MolGR worker exited without a result"
                    ) from error
            else:
                process.kill()
                process.join()
                pipe.close()
                pipe = None
                process = None
                label, reason, witness = (
                    "Inconclusive",
                    f"Comparison timed out after {timeout} s",
                    (),
                )
            row = (case, published, label, reason)
            if include_smiles:
                reference = decode_reference(encoded)
                try:
                    row += mapped_pair(smiles, reference, witness)
                except Exception as error:
                    raise ValueError(
                        f"{case}: could not export mapped SMILES"
                    ) from error
            yield row
    finally:
        if process is not None:
            if process.is_alive():
                process.kill()
            process.join()
        if pipe is not None:
            pipe.close()


def records_for(parts, requested):
    found = set()
    for part in parts:
        if requested and not any(
            part["first"] <= case <= part["last"] for case in requested
        ):
            continue
        with (ROOT / part["sdf"]).open("rb") as stream:
            for mol in Chem.ForwardSDMolSupplier(
                stream, sanitize=False, removeHs=False
            ):
                if mol is None:
                    raise ValueError(f"Unreadable reference in {part['sdf']}")
                case = mol.GetProp("CASE_ID")
                if requested and case not in requested:
                    continue
                found.add(case)
                yield (
                    case,
                    mol.GetProp("UPSTREAM_COMPARISON"),
                    mol.GetProp("UPSTREAM_SMILES")
                    if mol.HasProp("UPSTREAM_SMILES")
                    else "",
                    mol.GetProp("NATIVE_GRAPH_ZLIB_BASE64"),
                )
    if missing := requested - found:
        raise ValueError(f"No scoring reference for {', '.join(sorted(missing))}")


def source_records(parts, requested, csd_csv, tmqmg_csv=None):
    sources = {row["IDs"]: (line, row) for line, row in csd_rows(csd_csv)}
    tmqmg_sources = None
    fields = ("smiles_csd_api", "smiles_csd_api_fixed", "huckel_smiles_CSD_xyz")
    for part in parts:
        if requested and not any(
            part["first"] <= case <= part["last"] for case in requested
        ):
            continue
        references = records_for([part], set())
        with (ROOT / part["metadata"]).open() as meta:
            for row in map(json.loads, meta):
                if not row["benchmark"]:
                    continue
                case, _, _, encoded = next(references)
                if case != row["id"]:
                    raise ValueError(f"{case}: reference and metadata out of order")
                if requested and case not in requested:
                    continue
                evidence = row.get("source_evidence") or {}
                ref = evidence.get("matched_source_ref")
                support_field = (row.get("csd_support") or {}).get("field")
                legacy_smiles = evidence.get("matched_source_smiles")
                if not (ref or support_field or legacy_smiles):
                    continue
                if ref and evidence["source_table"] == "tmqmg_smiles":
                    if tmqmg_sources is None:
                        path = tmqmg_csv or csd_csv.with_name("tmqmg_smiles.csv")
                        digest = "1db086d9398dd65cf7de44f3fc36d3087b433ef9bafcc3317b081142a2d8423b"
                        tmqmg_sources = {
                            r["IDs"]: (i, r)
                            for i, r in csd_rows(path, expected_sha256=digest)
                        }
                    line, source = tmqmg_sources[case]
                    label = "tmqmg_smiles." + ref["field"]
                else:
                    if ref and evidence["source_table"] != "csd_smiles":
                        raise ValueError(f"{case}: unknown source table")
                    line, source = sources[case]
                    label = ""
                if ref:
                    smiles = source[ref["field"]]
                    if (
                        line != ref["row"]
                        or hashlib.sha256(smiles.encode()).hexdigest() != ref["sha256"]
                    ):
                        raise ValueError(f"{case}: source row or SMILES digest differs")
                else:
                    smiles = legacy_smiles or source[support_field]
                if not label:
                    label = "|".join(
                        name for name in fields if source.get(name) == smiles
                    )
                yield case, label or "archived source SMILES", smiles, encoded


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ids", nargs="*", help="CSD IDs (default: all scoring targets)")
    parser.add_argument(
        "--csd-csv",
        type=Path,
        help="compare documented source SMILES using the pinned CSD CSV",
    )
    parser.add_argument(
        "--tmqmg-csv",
        type=Path,
        help="pinned tmqmg_smiles.csv (default: beside --csd-csv)",
    )
    args = parser.parse_args()
    if args.tmqmg_csv and not args.csd_csv:
        parser.error("--tmqmg-csv requires --csd-csv")
    manifest = json.loads((ROOT / "manifest.json").read_bytes())
    requested = {case.upper() for case in args.ids}
    parts = manifest["parts"]
    writer = csv.writer(sys.stdout, lineterminator="\n")
    if args.csd_csv:
        records = source_records(parts, requested, args.csd_csv, args.tmqmg_csv)
        writer.writerow(("id", "source_fields", "recomputed", "reason"))
    else:
        records = records_for(parts, requested)
        writer.writerow(
            (
                "id",
                "published",
                "recomputed",
                "reason",
                "upstream_mapped_smiles",
                "reference_mapped_smiles",
            )
        )
    seen = set()
    for row in comparisons(
        records, manifest["benchmark"], include_smiles=not args.csd_csv
    ):
        writer.writerow(row)
        seen.add(row[0])
    if missing := requested - seen:
        raise ValueError(f"No comparison available for {', '.join(sorted(missing))}")
    print(f"MolGR {version('molgr')}: {len(seen)} comparisons.", file=sys.stderr)


if __name__ == "__main__":
    main()
