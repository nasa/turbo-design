"""
Utilities to read streamline-wise compressor data from ``E3_HPC_Overall_Data.xlsx``.

The workbook groups streamlines by component in one sheet. This helper reads the
first sheet by default and automatically splits blocks whenever the ``SL`` column
resets to 1 (12 streamlines per block). It then labels the blocks as inlet,
rotor1/stator1, …, up to rotor10/stator10, and ignores any turbine/OGV data.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import pandas as pd

def load_overall_data(
    workbook_path: Path | str | None = None,
    sheet_name: int | str | None = None,
    pt_scale: float | None = None,
    tt_scale: float | None = None,
) -> Dict[str, pd.DataFrame]:
    """
    Load the E3 HPC overall data sheet and split it into per-component blocks
    based on streamline resets.

    Args:
        workbook_path: Path to ``E3_HPC_Overall_Data.xlsx``. Defaults to the file
            next to this script.
        sheet_name: Excel sheet index or name to read. If ``None`` (default),
            the reader searches for the first sheet containing an ``SL`` column.
        pt_scale: Optional divisor applied to ``Inlet Pt``/``PT Exit`` columns to
            convert reported units to the desired basis (e.g., divide by 8.65 or
            10.399 to non-dimensionalize). If ``None``, no scaling is applied.
        tt_scale: Optional divisor applied to ``TT Inlet``/``TT Exit`` columns
            (e.g., divide by 469.5249331 or other reference temperature). If
            ``None``, no scaling is applied.

    Returns:
        Dict mapping component name (e.g., ``rotor1``) to a pandas DataFrame with
        the original columns (SL, Inlet Pt, PT Exit, TT Inlet, TT Exit, Ma Inlet,
        Ma Exit, Ts Exit / Ts Inlet, Ps Exit / Ps Inlet, ds, PR, OTAC TT, Radius, Beta).
    """
    if workbook_path is None:
        workbook_path = Path(__file__).resolve().parent / "E3_HPC_Overall_Data.xlsx"
    workbook_path = Path(workbook_path)

    # Choose sheet: auto-detect SL column if not provided or missing
    xl = pd.ExcelFile(workbook_path)
    sheet_to_use = sheet_name if sheet_name is not None else xl.sheet_names[0]

    def _contains_sl(name) -> bool:
        temp = xl.parse(name, nrows=1)
        return "SL" in temp.columns

    if sheet_name is None or not _contains_sl(sheet_to_use):
        target_sheet = None
        for name in xl.sheet_names:
            if _contains_sl(name):
                target_sheet = name
                break
        sheet_to_use = target_sheet if target_sheet is not None else xl.sheet_names[0]

    df = pd.read_excel(workbook_path, sheet_name=sheet_to_use, header=0)

    if "SL" not in df.columns:
        raise ValueError(f"No sheet with an 'SL' column found in {workbook_path}")

    # Identify block boundaries where SL resets to 1
    sl_series = df["SL"]
    reset_indices: List[int] = [i for i, val in enumerate(sl_series) if val == 1]
    reset_indices.append(len(df))  # sentinel for last block

    # Build labels: inlet, rotor1/stator1 ... rotor10/stator10
    labels = ["inlet"]
    for stage in range(1, 11):
        labels.extend([f"rotor{stage}", f"stator{stage}"])

    blocks: Dict[str, pd.DataFrame] = {}
    for block_idx in range(len(reset_indices) - 1):
        name = labels[block_idx] if block_idx < len(labels) else f"section{block_idx}"
        start = reset_indices[block_idx]
        end = reset_indices[block_idx + 1]
        block = df.iloc[start:end].copy()
        blocks[name] = block.reset_index(drop=True)

    # Optional scaling of pressures/temperatures
    if pt_scale is not None:
        for col in ["Inlet Pt", "PT Exit"]:
            if col in df.columns:
                for name, block in blocks.items():
                    block[col] = block[col] / pt_scale
    if tt_scale is not None:
        for col in ["TT Inlet", "TT Exit"]:
            if col in df.columns:
                for name, block in blocks.items():
                    block[col] = block[col] / tt_scale

    return blocks


def load_blade_counts(workbook_path: Path | str | None = None) -> Dict[str, int]:
    """
    Read number of vanes/blades per row from the ``Design Parameters`` sheet.

    Returns keys like ``igv``, ``rotor1``, ``stator1``, …, ``rotor10``, ``stator10``.
    """
    if workbook_path is None:
        workbook_path = Path(__file__).resolve().parent / "E3_HPC_Overall_Data.xlsx"
    workbook_path = Path(workbook_path)

    df = pd.read_excel(workbook_path, sheet_name="Design Parameters", header=None)
    stage_col = df.iloc[:, 9]
    count_col = df.iloc[:, 10]

    counts: Dict[str, int] = {}
    for stage, count in zip(stage_col, count_col):
        if pd.isna(stage) or pd.isna(count):
            continue
        name = str(stage).strip().lower().replace(" ", "")
        try:
            counts[name] = int(float(count))
        except (TypeError, ValueError):
            continue
    return counts


if __name__ == "__main__":
    data = load_overall_data()
    counts = load_blade_counts()
    for name, block in data.items():
        print(f"{name.upper()} (rows={len(block)})  blades={counts.get(name, 'n/a')}")
        print(block.head(3))
        print()
