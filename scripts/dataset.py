"""Shared paths and transition-metal elements for dataset audits."""

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "evidence" / "tables"
TM = set(range(21, 31)) | set(range(39, 49)) | set(range(72, 81)) | {57, 71}


def explicit_h_parser(*, sanitize: bool) -> "Chem.SmilesParserParams":
    from rdkit import Chem

    params = Chem.SmilesParserParams()
    # The RDKit binary accepts bools. Its bundled property setters omit their types.
    params.removeHs = False  # type: ignore[assignment]
    params.sanitize = sanitize  # type: ignore[assignment]
    return params
