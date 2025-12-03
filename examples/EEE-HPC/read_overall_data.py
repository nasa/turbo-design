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

import numpy as np
import pandas as pd

def load_overall_data(
    workbook_path: Path | str | None = None,
    sheet_name: int | str | None = None,
    pt_scale: float | None = None,
    tt_scale: float | None = None,
    loss_sheet_name: int | str | None = None,
    convert_units: bool = True,
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
        loss_sheet_name: Optional sheet containing loss/incidence/deviation. If
            provided, those columns will be merged into the returned blocks when
            present (columns containing ``loss``, ``incidence``, ``deviation``).
        convert_units: When True (default), convert Pt columns from psi to Pa and
            TT columns from Rankine to Kelvin.

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

    def _split_blocks(frame: pd.DataFrame) -> Dict[str, pd.DataFrame]:
        frame = frame.copy()
        frame["SL"] = pd.to_numeric(frame["SL"], errors="coerce")
        frame = frame[frame["SL"].notna()].reset_index(drop=True)

        sl_series = frame["SL"]
        reset_indices: List[int] = [i for i, val in enumerate(sl_series) if val == 1]
        labels = ["inlet"]
        for stage in range(1, 11):
            labels.extend([f"rotor{stage}", f"stator{stage}"])
        blocks: Dict[str, pd.DataFrame] = {}
        if not reset_indices:
            blocks["inlet"] = frame.copy().reset_index(drop=True)
            return blocks

        reset_indices.append(len(frame))
        for idx in range(len(reset_indices) - 1):
            start = reset_indices[idx]
            end = reset_indices[idx + 1]
            name = labels[idx] if idx < len(labels) else f"section{idx}"
            block = frame.iloc[start:end].copy().reset_index(drop=True)
            blocks[name] = block
        return blocks

    main_blocks = _split_blocks(df)

    # Optional merge of loss/incidence/deviation from another sheet
    if loss_sheet_name is not None:
        df_loss = pd.read_excel(workbook_path, sheet_name=loss_sheet_name, header=0)

        # Try to normalize SL column and rename loss/incidence/deviation
        if "SL" not in df_loss.columns:
            for col in df_loss.columns:
                if str(df_loss.iloc[0][col]).strip().lower() == "sl":
                    df_loss = df_loss.rename(columns={col: "SL"})
                    df_loss = df_loss.iloc[1:].reset_index(drop=True)
                    break
        else:
            if str(df_loss.iloc[0]["SL"]).strip().lower() == "sl":
                df_loss = df_loss.iloc[1:].reset_index(drop=True)

        rename_map = {}
        for col in df_loss.columns:
            col_l = str(col).lower()
            first_val = str(df_loss.iloc[0][col]).lower() if len(df_loss) else ""
            if "loss" in col_l or "loss" in first_val:
                rename_map[col] = "Loss"
            elif "incidence" in col_l or "incidence" in first_val:
                rename_map[col] = "Incidence"
            elif "deviation" in col_l or "deviation" in first_val:
                rename_map[col] = "Deviation"
        df_loss = df_loss.rename(columns=rename_map)

        if "SL" in df_loss.columns:
            loss_blocks = _split_blocks(df_loss)
            for name, lblock in loss_blocks.items():
                target = main_blocks.get(name)
                if target is None:
                    continue
                for col in ("Loss", "Incidence", "Deviation"):
                    if col in lblock.columns:
                        target[col] = pd.to_numeric(lblock[col], errors="coerce")

    # Identify block boundaries where SL resets to 1
    # Optional scaling of pressures/temperatures
    if pt_scale is not None:
        for col in ["Inlet Pt", "PT Exit"]:
            if col in df.columns:
                for name, block in main_blocks.items():
                    block[col] = block[col] / pt_scale
    if tt_scale is not None:
        for col in ["TT Inlet", "TT Exit"]:
            if col in df.columns:
                for name, block in main_blocks.items():
                    block[col] = block[col] / tt_scale

    if convert_units:
        psi_to_pa = 6894.757
        r_to_k = 5.0 / 9.0
        for name, block in main_blocks.items():
            for col in ["Inlet Pt", "PT Exit"]:
                if col in block.columns:
                    block[col] = pd.to_numeric(block[col], errors="coerce") * psi_to_pa
            for col in ["TT Inlet", "TT Exit"]:
                if col in block.columns:
                    block[col] = pd.to_numeric(block[col], errors="coerce") * r_to_k

    return main_blocks


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


def compute_entropy_rise(
    workbook_path: Path | str | None = None,
    sheet_name: str = "ds Manual Calculations",
    fluid: "Solution | None" = None,
) -> Dict[str, pd.DataFrame]:
    """
    Compute entropy rise Cp*ln(T2/T1) - R*ln(P2/P1) for each component/block.

    If a Cantera ``fluid`` is provided, Cp and R are evaluated at each inlet
    streamline state using the sheet's static temperature/pressure columns. Sheet
    units are expected to be Rankine (temperature) and psi (pressure).
    Otherwise, Cp/R columns from the sheet (``Cps``, ``Rs``) are used.
    """
    if workbook_path is None:
        workbook_path = Path(__file__).resolve().parent / "E3_HPC_Overall_Data.xlsx"
    workbook_path = Path(workbook_path)

    df = pd.read_excel(workbook_path, sheet_name=sheet_name, header=0)
    if "SL" not in df.columns:
        raise ValueError(f"Sheet '{sheet_name}' does not contain an 'SL' column.")

    def _split_blocks(frame: pd.DataFrame) -> Dict[str, pd.DataFrame]:
        frame = frame.copy()
        frame["SL"] = pd.to_numeric(frame["SL"], errors="coerce")
        frame = frame[frame["SL"].notna()].reset_index(drop=True)

        sl_series = frame["SL"]
        reset_indices: List[int] = [i for i, val in enumerate(sl_series) if val == 1]
        labels = ["inlet"]
        for stage in range(1, 11):
            labels.extend([f"rotor{stage}", f"stator{stage}"])
        blocks: Dict[str, pd.DataFrame] = {}
        if not reset_indices:
            blocks["inlet"] = frame.copy().reset_index(drop=True)
            return blocks

        reset_indices.append(len(frame))
        for idx in range(len(reset_indices) - 1):
            start = reset_indices[idx]
            end = reset_indices[idx + 1]
            name = labels[idx] if idx < len(labels) else f"section{idx}"
            block = frame.iloc[start:end].copy().reset_index(drop=True)
            blocks[name] = block
        return blocks

    blocks = _split_blocks(df)
    for name, block in blocks.items():
        ts_ratio = pd.to_numeric(block.get("Ts Exit / Ts Inlet", np.nan), errors="coerce")
        ps_ratio = pd.to_numeric(block.get("Ps Exit / Ps Inlet", np.nan), errors="coerce")
        if fluid is not None and "Ts" in block.columns and "Ps" in block.columns:
            T_rankine = pd.to_numeric(block["Ts"], errors="coerce")
            P_psi = pd.to_numeric(block["Ps"], errors="coerce")
            T_kelvin = T_rankine * (5.0 / 9.0)
            P_pa = P_psi * 6894.757
            cp_vals = []
            r_vals = []
            for T_k, P in zip(T_kelvin, P_pa):
                try:
                    fluid.TP = float(T_k), float(P)
                    cp_vals.append(fluid.cp_mass)
                    r_vals.append(fluid.cp - fluid.cv)
                except Exception:
                    cp_vals.append(np.nan)
                    r_vals.append(np.nan)
            cps = np.array(cp_vals)
            rs = np.array(r_vals)
        else:
            cps = pd.to_numeric(block.get("Cps", np.nan), errors="coerce")
            rs = pd.to_numeric(block.get("Rs", np.nan), errors="coerce")
        block["ds_calc"] = cps * np.log(ts_ratio) - rs * np.log(ps_ratio)
        blocks[name] = block
    return blocks


if __name__ == "__main__":
    data = load_overall_data()
    counts = load_blade_counts()
    for name, block in data.items():
        print(f"{name.upper()} (rows={len(block)})  blades={counts.get(name, 'n/a')}")
        print(block.head(3))
        print()
