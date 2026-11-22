# Data format and reproduction

The manifest lists 50 groups in ascending CSD ID order. Group prefixes run from `tmqmg-corrected-01` to `tmqmg-corrected-50`. The first 47 groups contain 1,491 inputs each and the last three contain 1,490 each. Join files by ID. XYZ and metadata cover every input, while SDF covers the retained subset. [File manifest](../manifest.json)

| File | Identifier | Contents |
| --- | --- | --- |
| `evidence/inputs/*.xyz` | `id` on line two | Original elements in original order, followed by Cartesian coordinates in Å. A comment such as `id=ABEVAG charge=0 multiplicity=1` identifies the input state. Restore an empty comment to reproduce the upstream XYZ bytes. |
| `data/*.sdf` | Title and `CASE_ID` | V3000 graph with explicit H atoms, coordinates, integer bonds, formal charges, radicals and coordination bonds. |
| `evidence/inputs/*.jsonl` | `id` | One record per ID with input state, upstream SMILES, scoring eligibility, evidence tier, annotation or exclusion rationale, sources and available calculations. |
| `manifest.json` | Paths and `sha256` | Group boundaries and counts, file digests, upstream input hashes, calculation settings and frozen control distributions. |

## Evidence table guide

Join the tables by CSD `id`. The two graph comparison tables use different sources. The upstream comparison measures agreement with tmQMg. The source comparison checks the released graph against its documented structural source.

| Table | Coverage | When to use it |
| --- | ---: | --- |
| [excluded-records.csv](tables/excluded-records.csv) | 2,322 excluded IDs | Look up why an input was excluded, follow its source and locate its XYZ and case record. Supplies the counts for the README input overview. |
| [boron-exclusion-audit.csv](tables/boron-exclusion-audit.csv) | 1,189 original boron exclusions | Distinguish cage and bridge topologies using XYZ contacts and same-ID source connectivity. For the 16 bridge classifications, the JSON field `bridge_review` gives the upstream and source hydrogen inventories and bond assignments. |
| [boron-hydrogen-bridges.jsonl](tables/boron-hydrogen-bridges.jsonl) | 565 candidates | Inspect atom indices, distances, WBOs and source bridge assignments. Of 516 supported cases, 354 were newly excluded and 162 already had another exclusion. The other 49 have geometric contacts alone. |
| [equivalence-rescore.csv](tables/equivalence-rescore.csv) | All 72,225 released graphs | Compare the upstream and reference explicit-H mapped SMILES, read the verdict and reason, or select the inequivalent pairs. |
| [upstream-issue-types.csv](tables/upstream-issue-types.csv) | 4,577 flagged IDs | Look up the seven issue categories. Joining to the retained IDs in the equivalence table gives the 4,560 complexes counted in the README chart. |
| [ligand-edits.jsonl](tables/ligand-edits.jsonl) | 4,915 inequivalent pairs | Locate changes in ligand skeletons, bonds, attached hydrogens and formal charges using atom indices and the two ligand forms. |
| [metal-valence-audit.csv](tables/metal-valence-audit.csv) | 71,985 metal-centre comparisons | Recheck the two metal-valence categories using formal charges, terminal hydrides, allowed valences and source support. Includes pairs that passed the screen. |
| [source-comparisons.csv](tables/source-comparisons.csv) | 62,409 IDs with source SMILES | Check which source field was compared with the released graph and whether the comparison established equivalence. |
| [example-bonds.csv](tables/example-bonds.csv) | 17 bonds in six examples | Verify the illustrated bond assignments using XYZ indices, both written bond orders, distances and WBOs. CUDBOU illustrates alternating atomic charges. |
| [example-gxtb-rerun.csv](tables/example-gxtb-rerun.csv) | Seven README examples | Verify that fresh single points reproduce the archived WBO files, energies and orbital gaps. |

## Verification of the 2026-11-20 release

| Check | Scope | Result |
| --- | ---: | --- |
| Input identity | 74,547 IDs | Original XYZ bytes, upstream SMILES and input charges match the pinned tmQMg files. All 72,225 SDF records preserve XYZ atom order and coordinates. Each native graph matches its SHA-256 digest. |
| Exclusions | 2,322 IDs | Rebuilt every assessment from its recorded inputs and sources. The boron topology and hydrogen-bridge audits reproduce. For the 492 H-deficient IDs with recorded site mappings, the nearest heavy atom in XYZ confirms the input H count at all 890 sites. Source indices are available for 888 sites across 491 IDs. Element and H count agree at every indexed site, using original CSD connectivity for 453 IDs and revised CSD connectivity for 38. |
| Archived single points | 72,628 runs on 72,581 IDs | During release preparation we checked the extracted options, energy and orbital gap against the original logs and the WBO pairs against the original WBO files. We also verified the input geometry hashes, charges and UHF settings. The public validator reconstructs the original WBO file bytes from every full-precision pair list and verifies its digest. The cohorts contain 61,737 main runs, 2,200 controls and 8,691 additional controls. |
| Control distributions | 2,200 control IDs | Refit the bond-length and WBO quantiles on 1,757 IDs and reapplied the rules to 443 held-out IDs. All class counts, quantiles, the nonbonded contact threshold and holdout counts reproduce. |
| Bond evidence | 3,935,782 nonmetal covalent bonds | Recomputed every bond length from XYZ and looked up its WBO in the archived single point. The written orders, assessments, per-ID counts and off-graph contacts agree. |
| Upstream equivalence | 72,225 retained IDs | 66,706 equivalent, 4,915 inequivalent, 360 inconclusive, 240 missing upstream labels and four unreadable labels. All seven issue categories and the paired illustrations were regenerated. [Per-ID results](tables/equivalence-rescore.csv) |
| Source graph equivalence | 62,409 retained IDs with a documented source SMILES | 62,351 equivalent and 58 inconclusive. Of the inconclusive cases, 36 lacked a resonance witness within the supported rules, 18 had unresolved radical differences and four reached the atom-mapping limit. [Source fields and verdicts](tables/source-comparisons.csv) |
| Independent g-xTB runs | Seven illustrated complexes | All seven fresh WBO files match their archived SHA-256 digests. The reported energies and orbital gaps also match at the precision in the logs. [Run results](tables/example-gxtb-rerun.csv) |

The electronic-state decisions cite published measurements or assignments. The exclusion rerun checks the recorded input and cited assessment. For JETGIY, the relationship between the optimized input and the reported high-spin crystal remains uncertain. The per-case records describe the limits of each source. [Exclusion evidence](tables/excluded-records.csv)

## Reference graph fields

`GRAPH_SOURCE`, `SOURCE_REFERENCES` and `SELECTION_REASON` describe the assignment and its sources. Read the JSONL for the scope of each source. `UPSTREAM_SMILES` is present where available. `UPSTREAM_COMPARISON` records the historical equivalence verdict. `EVIDENCE_TIER`, `SPIN_REVIEW` and the optional `EVIDENCE_NOTES` describe confidence and limitations. `INPUT_CHARGE` and `INPUT_MULTIPLICITY` specify the calculation input. `GXTB_SUMMARY` and `WBO_SUMMARY` summarize the selected calculation. For scoring, decode and decompress `NATIVE_GRAPH_ZLIB_BASE64` to recover the original RDKit graph with its atom-level electronic properties. `GRAPH_SHA256` identifies the decompressed graph bytes and `NATIVE_RDKIT_VERSION` records the serialization version.

HUPQUG, LACGEB and QOFJIE also had a spin-coupling note inside the serialized RDKit molecule. Translating that property changed their serialized digests. `SOURCE_NATIVE_GRAPH_SHA256` (JSONL `source_native_graph_sha256`) retains the previous digest for each. Restoring the original property value reproduced the previous bytes in all three cases. All 72,225 V3000 connection tables retain their pre-translation bytes.

In JSONL, filter on `benchmark` before scoring. `tier` describes a released graph, and `exclusion` identifies a withheld one. Depending on the case, `csd_support`, `source_evidence`, `literature_scope`, `site_evidence`, `reviewed_bonds`, and the hydrogen and spin review fields provide additional evidence. A link to the input dataset alone identifies the input. Consult the evidence tier and source description for support for the complete graph. “Not independently verified” in the spin review means the reference uses the input multiplicity.

## Excluded inputs

All 2,322 exclusions have `benchmark=false` and `tier=null`. The index links each ID to the relevant XYZ and JSONL. The assessment includes a reason, evidence, input formula and sources. Four historical exclusion labels map to the three categories on the README. [Exclusion index](tables/excluded-records.csv)

| Archival category | IDs | Evidence |
| --- | ---: | --- |
| Missing input hydrogens | 770 | The source and XYZ have identical heavy-element counts. For 769 IDs the comparison uses the metal-containing component in the same-ID CSD transcription. EYINEF uses the published complex formula, excluding lattice water, and COD 2014185. The assessment records `source_formula`, `input_H`, `expected_H` and `missing_H`. |
| Boron cluster | 1,189 | The geometry review records B atoms, B···B contacts and H atoms close to two B atoms. In 1,181 cases three B atoms are pairwise close, and a further four have B–H–B contacts. FIGROC and YIBNIG have evidence for boron–carbon cages. NAHNOB and QEJYAH have same-ID evidence for B–H–metal bridges. |
| Questionable electronic state | 9 | Eight singlet inputs conflict with a reported spin state. The report for the composition of JETGIY describes high-spin Fe(II), but independent same-ID crystal correspondence remains unresolved. Each record retains the input multiplicity, reported state, evidence and citations. |
| Boron hydrogen bridge | 354 | Of the previously eligible IDs, 310 have an explicit same-ID source bridge, and 44 have both short contacts and supporting WBOs from the same calculation. |

For example, ABATUS has 15 H atoms in XYZ and 17 in the source. Its site mapping identifies the deficit. FOQBOC has singlet input multiplicity, while the source reports a triplet–quintet transition for compound 9. EDAKUP is retained in group 09 with the upstream aromatic boron-ring representation and tier-3 electronic assignment. [CSD formulas and connectivity](https://github.com/jensengroup/xyz2mol_tm/blob/db2734b084045ef895528c7a3d5040d38da23983/SMILES_csvs/csd_smiles.csv) · [FOQBOC primary report](https://doi.org/10.1016/S0022-328X%2899%2900354-X) · [EDAKUP case metadata](inputs/tmqmg-corrected-09.jsonl)

### Boron and bridging hydrogens

The bridge rule concerns one H shared by two B atoms or by B and a metal. These multicentre bonds fall outside the two-centre Lewis convention used for this release. Terminal B–H groups, Si–H–metal contacts, C–H···metal contacts and shared protons have separate chemical interpretations. The geometry fields `nearest_BB` and `BB_contacts` give one-based XYZ B indices and distances, with a 2.10 Å cutoff for the contact list. `BBB_triangle` identifies a triple satisfying the cutoff on all sides. `BHB_contacts` records B, H and B indices with both distances ≤1.50 Å. These measurements locate candidate sites.

The full bridge scan found 472 IDs with a shared H in the original same-ID CSD graph. Seven of the 310 formerly eligible source-supported IDs fall outside the optimized-XYZ distance window. A further 44 IDs meet the geometry and WBO criteria without an explicit source bridge. For the geometry-and-WBO criterion, each B···H distance must be ≤1.75 Å. In B–H–metal contacts the H···metal distance must also be at most the RDKit metal covalent radius plus 0.95 Å. Both WBOs must reach 0.10 in the same converged calculation at the input charge and UHF setting. An explicit bridge in the original same-ID source graph independently satisfies the source criterion. Of the 516 supported bridge cases, 160 already had a boron exclusion and two had missing hydrogens. Another 49 IDs have geometric contacts alone, including 15 retained complexes. BAJHOK and SEPKAA have source component formulas that differ from the input, with geometry and both WBOs also supporting a bridge.

`csd_sites` uses one-based indices in the CSD SMILES. `xyz_sites` gives one-based B, H and partner indices, partner element, both distances in Å and both WBOs. JSON `null` marks missing data. The 354 withdrawn graphs retain their former tier and review in `previous_graph_review`. [Per-ID bridge evidence](tables/boron-hydrogen-bridges.jsonl)

A separate topology check on the 1,189 historical boron exclusions found 1,096 source graphs with edge-sharing B/C triangular faces, 14 with B–H–B bridges, and two with B–N–B ligands bridged to a metal. The other 77 lack a readable original same-ID connection table. Of these, 75 have a B atom within 2.10 Å of at least five other B atoms. DITMOM and ZITQOM each have a triangular B₃ contact and two H atoms close to pairs of B atoms. The 16 entries classified by bridge topology all have upstream labels with the correct composition and total charge but incorrect internal connectivity. The upstream TIKBIC label disconnects a BH₂⁺ fragment from the complex. VUGBII loses an edge of its B₃ triangle. NAHNOB and QEJYAH assign ligand hydrogens to Sc and Y, respectively. Each of the 16 rows includes the upstream SMILES, hydrogen inventories and distances in `bridge_review`, alongside the source link in `csd_row`. [Classifications and paired graph checks](tables/boron-exclusion-audit.csv) · [Ir borane](https://doi.org/10.1021/om049960a) · [W borane](https://doi.org/10.3390/inorganics11060248) · [Sc compound](https://doi.org/10.1039/C0CC02810J) · [Y compound](https://patents.google.com/patent/US20100168404A1/en)

Recompute both boron evidence tables from the pinned CSD CSV, original XYZ and case records. The exclusion script reads the bridge index and rebuilds the 2,322 assessments from coordinates, pinned source formulas and recorded case evidence. Its output pairs each ID with the assessment stored in `exclusion_evidence.assessment`. Run from the repository root.

```bash
uv run --locked python scripts/audit_boron.py --csd-csv /path/to/csd_smiles.csv --output-dir boron-recomputed
cmp boron-recomputed/boron-hydrogen-bridges.jsonl evidence/tables/boron-hydrogen-bridges.jsonl
cmp boron-recomputed/boron-exclusion-audit.csv evidence/tables/boron-exclusion-audit.csv
uv run --locked python scripts/exclusions.py --csd-csv /path/to/csd_smiles.csv --check-sites > exclusions-recomputed.jsonl
```

## Bond evidence

g-xTB v2.0.0 single points use fixed input coordinates and the recorded charge and UHF setting. All 72,225 retained IDs and 356 excluded IDs have converged calculations. Some IDs have several runs, for a total of 72,628. `gxtb[]` preserves each run's options, energy in Eh, HOMO–LUMO gap in eV, WBO file digest and reported atom-pair values, including pairs involving a metal. An empty list means no archived calculation. [g-xTB release](https://github.com/grimme-lab/g-xtb/releases/tag/v2.0.0)

The repository distributes the extracted results in `evidence/inputs/*.jsonl`. The versioned companion ZIP preserves the original console output, including iteration histories and diagnostic messages. To check a published bond assignment, use the XYZ and the full-precision WBO pairs. To obtain a fresh calculation for comparison, use the single-point command below. [Original g-xTB logs](#original-g-xtb-logs)

Each row of `wbo_analysis.bonds` gives two one-based XYZ indices, the element pair, written SDF bond order, control class, XYZ distance in Å, distance assessment, WBO and WBO assessment. A missing WBO is `null`. The full-precision value is in the selected cohort of `gxtb[].wbo_pairs`. `unlisted_contacts` records off-graph nonmetal pairs above `contact_threshold`. The hydrogen-site audit uses zero-based indices in `exclusion_evidence.site_mapping`. Its `xyz_idx` addresses XYZ and its optional `source_idx` addresses the source SMILES.

Only the 2,200 runs in the `control` cohort enter the calibration. The holdout consists of IDs for which the integer represented by the first six hex digits of SHA-256 of ASCII `holdout` followed by the ID is divisible by five. The validator refits the archived calibration and checks each held-out bond.

The geometry and WBO controls require at least 15 IDs and 40 bonds of a given element pair and bond class. The empirical range spans the 0.5th to 99.5th percentiles. Within this range, a WBO outside the 5th–95th percentile interval and inside the interquartile range of another supported bond class for the same element pair is flagged as closer to another bond class. The off-graph threshold (WBO > 0.3578996052233558) is the 99th percentile of the largest unlisted nonmetal-pair WBO per fit complex, using zero for a complex with no such reported pair. Aromatic environments use a separate control. The V3000 entry records the actual assigned bond order even when the control class is aromatic.

To repeat the independent calculations for the seven illustrations, use the pinned executable and its accompanying parameter files. Each job uses one thread and its own working directory. Omitting the IDs selects all archived runs. The release-wide bond analysis above uses the archived results.

```bash
python3 scripts/rerun_gxtb.py --xtb /path/to/xtb-6.7.1/bin/xtb --jobs 4 DEFWIY JEGLOA RAWVAN ERULUB ABEVAG SAZFEF CUDBOU > example-gxtb-rerun.csv
cmp example-gxtb-rerun.csv evidence/tables/example-gxtb-rerun.csv
```

## Release bundle

The local distribution for tag `v2026-11-20` is `dist/tmqmg-corrected-2026-11-20/`. Open its `README.md` for links to the SDF archive, the repository snapshot, the original log ZIP and the standalone log index. The snapshot includes all 50 data groups, the evidence tables, documentation, scripts and dependency lockfile. Its contents come from the tagged Git tree. Run `sha256sum -c SHA256SUMS` in the bundle directory before extracting the files.

To assemble the same bundle from a checkout of the tag, place the verified g-xTB ZIP and index in `dist/tmqmg-corrected-gxtb-logs-2026-11-20/`, then run `python3 scripts/package.py --bundle` from the repository root. Extracting just the SDF archive provides all 72,225 released graphs and the attribution notices.

### Publish a GitHub Release

The release workflow runs on a `vYYYY-MM-DD` tag push or by manual dispatch with an existing tag. It checks the tag against the manifest, validates the dataset, and builds the SDF and repository snapshot from the tagged commit. The original log ZIP and its index come from the local release bundle. [Release workflow](../.github/workflows/release.yml)

Once the workflow is on the default branch and the tag is on GitHub, stage the two log files in the release draft and start the publication run. A tag push that includes the workflow creates the draft automatically. The commands below create the draft if needed, upload the two log files and dispatch the workflow from the default branch. Run them from a clone of the destination repository that has the local bundle.

```bash
tag=v2026-11-20
bundle="dist/tmqmg-corrected-${tag#v}"
logs="tmqmg-corrected-gxtb-logs-${tag#v}"

if ! gh release view "$tag" > /dev/null; then
  gh release create "$tag" --draft --verify-tag --title "tmQMg-Corrected ${tag#v}" --notes "Dataset and reproduction files for $tag."
fi
gh release upload "$tag" "$bundle/$logs.zip" "$bundle/$logs-index.csv"
gh workflow run release.yml -f tag="$tag"
```

The workflow publishes after the log files match the manifest and GitHub’s asset digests match all six local bundle files. A failed check leaves the draft available for correction and retry. Runs for an already published tag leave its assets intact.

## Original g-xTB logs

The companion bundle for release `2026-11-20` contains `tmqmg-corrected-gxtb-logs-2026-11-20.zip` and `tmqmg-corrected-gxtb-logs-2026-11-20-index.csv`, with digests for both in `SHA256SUMS`. The ZIP also contains an identical copy of the CSV as `index.csv`. The release manifest pins both digests in `gxtb_log_files_sha256`. The ZIP contains 72,628 original logs for 72,581 distinct CSD IDs. Look up a run by (`id`, `cohort`). The three cohorts are `main` (61,737), `control` (2,200) and `additional_control` (8,691).

Each index row gives the ZIP member path, the SHA-256 of the original log bytes, the source shard within the cohort, and the SHA-256 digests of the input XYZ and corresponding WBO file. The original EDAKUP job supplies its standalone log. All other logs come from the archived SQLite shards. During packaging, every log was matched to its release record by ID and cohort and checked for normal termination, input filename, run options, total energy and orbital gap. The source charge, UHF setting, XYZ digest and WBO digest were also checked. Each log was then read back from the ZIP and verified against its indexed SHA-256.

To check the bundle and read one log, run these commands in the bundle directory.

```bash
sha256sum -c SHA256SUMS
unzip -p tmqmg-corrected-gxtb-logs-2026-11-20.zip logs/main/EDAKUP.log | less
```

To rebuild from the original calculation stores, run the following command at the repository root. It writes a versioned directory under `dist/`.

```bash
python3 scripts/package_gxtb_logs.py \
  --main-dir /path/to/tmqmg-xyz-gxtb-20261119 \
  --additional-dir /path/to/tmqmg-pending-evidence \
  --edakup-dir /path/to/edakup-reference
```

## Read an ID

Run `uv sync --locked` and use `uv run --locked python` to execute this example.

```python
import json
from pathlib import Path
from rdkit import Chem

case = "ABEVAG"
parts = json.loads(Path("manifest.json").read_bytes())["parts"]
part = next(p for p in parts if p["first"] <= case <= p["last"])
with open(part["metadata"], encoding="utf-8") as stream:
    row = next(r for line in stream if (r := json.loads(line))["id"] == case)

if row["benchmark"]:
    with open(part["sdf"], "rb") as stream:
        molecules = Chem.ForwardSDMolSupplier(stream, sanitize=False, removeHs=False)
        mol = next(m for m in molecules if m.GetProp("CASE_ID") == case)
        print(case, row["tier"], mol.GetProp("SELECTION_REASON"))
else:
    print(case, row["selection_reason"])
```

## Check the Python scripts

Run mypy, Ruff lint and format checks, and Pyright from the repository root with the locked tools and the optional MolGR scoring dependency. The project Ruff configuration includes import sorting, comprehension and loop checks, subprocess warnings, and exception type checks.

```bash
uv run --locked --group checks --extra scoring mypy scripts
uv run --locked --group checks --extra scoring ruff check scripts
uv run --locked --group checks --extra scoring ruff format --check scripts
uv run --locked --group checks --extra scoring pyright scripts
```

## Validate the release

The standard-library validator checks the 150 file digests, group boundaries, IDs, row counts, XYZ digests, SDF atom order and coordinates, native graph hash, evidence tiers and exclusion assessments. It also recomputes nonmetal bond lengths, WBO assignments and off-graph contacts from the released XYZ, SDF and JSONL files. Optionally supply both pinned upstream files for an exact comparison of charge, SMILES and original XYZ bytes. [Input sources](methods.md)

```bash
python3 scripts/validate.py
python3 scripts/validate.py --upstream-csv /path/to/tmQMg_properties_and_targets.csv --upstream-xyz /path/to/tmQMg_xyz.zip
```

## Graph equivalence

Run MolGR 0.1.10 on one or more IDs, or omit the IDs for the full release. The scorer expands the upstream SMILES with RDKit and restores the native reference graph. The manifest specifies its comparison settings. A 60-second timeout gives an inconclusive verdict and the scorer continues with the next ID. The CSV retains the historical label alongside the new verdict and reason. [Full rescore](tables/equivalence-rescore.csv) · [MolGR 0.1.10](https://pypi.org/project/molgr/0.1.10/)

The last two columns of the upstream comparison table, `upstream_mapped_smiles` and `reference_mapped_smiles`, keep every hydrogen as a separate mapped atom, including hydrides and hydrogens implied by the upstream SMILES. Reference atom-map numbers are the one-based XYZ indices. Shared numbers identify corresponding atoms. The export uses the MolGR correspondence when it covers every atom and the upstream SMILES has no separate H nodes. Otherwise it matches heavy-atom fragments by element, isotope and connectivity, then pairs hydrogens at corresponding sites. Symmetry-equivalent sites can admit several valid assignments. An upstream atom without an established counterpart receives a number above the XYZ atom count. Missing upstream labels leave that column empty. For the four labels that fail RDKit sanitization, the column shows the readable raw connectivity and the scoring verdict stays `Unreadable label`.

```bash
uv run --locked --extra scoring python scripts/rescore.py DEFWIY JEGLOA RAWVAN ERULUB ABEVAG SAZFEF CUDBOU
uv run --locked --extra scoring python scripts/rescore.py > equivalence.csv
uv run --locked --extra scoring python scripts/rescore.py --csd-csv /path/to/csd_smiles.csv > source-comparisons.csv
```

Source comparison reads the Jensen Group CSVs at the pinned revision. Place both files in the same directory, or pass `--tmqmg-csv /path/to/tmqmg_smiles.csv` alongside `--csd-csv`. For `source_evidence.matched_source_ref` and `alternative_source_ref`, resolve the record by CSD ID and `source_table`. Each reference gives the CSV row (header is row 1), field name and SHA-256 of the exact UTF-8 string. An alternative reference records a different proposed assignment. Records with only `csd_support.field` use that field in the pinned `csd_smiles.csv`. [Source CSVs](https://github.com/jensengroup/xyz2mol_tm/tree/db2734b084045ef895528c7a3d5040d38da23983/SMILES_csvs)

When requesting individual IDs with `--csd-csv`, choose IDs from `source-comparisons.csv`. Of the seven illustrations, ABEVAG, CUDBOU, DEFWIY and SAZFEF have entries. For NAPKEW the documented comparison input is `smiles_CSD_fix` in `tmqmg_smiles.csv`. It reports which pinned fields contain the exact string and the independent equivalence verdict. A source can establish composition or topology even where MolGR leaves a complete graph comparison inconclusive.

## Rebuild tables, figures and the distribution archive

The first command reads the released upstream equivalence verdicts and rebuilds the ligand edits, metal valence audit and ID-level issue table. The second redraws both distributions and the seven graph pairs and writes `example-bonds.csv`. The package command builds `dist/tmqmg-corrected.tar.gz` with 50 SDF files, a version, SHA-256 checksums and license notices. After extraction, run `sha256sum -c SHA256SUMS` in the extracted directory. [Counting criteria](results.md) · [Data license](../sources/DATA-LICENSE.md)

```bash
uv run --locked --extra scoring python scripts/upstream_issues.py
uv run --locked python scripts/figures.py
python3 scripts/package.py
```
