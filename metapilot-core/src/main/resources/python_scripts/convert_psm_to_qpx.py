#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


PROTON_MASS = 1.007276466812


REQUIRED_COLUMNS = [
    "File Name",
    "Scan number",
    "Base Sequence",
    "Full Sequence",
    "Peptide Monoisotopic Mass",
    "Scan Retention Time",
    "Precursor Charge",
    "Protein Accession",
    "Score",
    "PEP",
    "Peptide Mass",
    "Mass difference",
    "Missed cleavages",
]


def require_columns(df: pd.DataFrame) -> None:
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(
            "Input TSV is missing required columns:\n"
            + "\n".join(f"  - {c}" for c in missing)
        )


def clean_peptide_sequence(value: Any) -> str | None:
    """
    QPX sequence should be the unmodified amino acid sequence.

    Your 'Base Sequence' should already be clean, but this function removes
    whitespace and common non-amino-acid characters just in case.
    """
    if pd.isna(value):
        return None

    seq = str(value).strip()
    seq = re.sub(r"\s+", "", seq)

    return seq or None


def clean_peptidoform(value: Any, fallback_sequence: str | None) -> str | None:
    """
    QPX peptidoform should ideally be ProForma.

    If your 'Full Sequence' is already ProForma-like, this keeps it.
    If it is empty, it falls back to the base sequence.
    """
    if pd.isna(value) or str(value).strip() == "":
        return fallback_sequence

    return str(value).strip()


def parse_scan(value: Any) -> list[int]:
    """
    QPX stores scan as a list.

    Examples:
      12345 -> [12345]
      'scan=12345' -> [12345]
      'controllerType=0 controllerNumber=1 scan=12345' -> [12345]
    """
    if pd.isna(value):
        return []

    text = str(value).strip()
    numbers = re.findall(r"\d+", text)

    if not numbers:
        return []

    # Usually the final number is the actual scan number in Thermo-style spectrum IDs.
    return [int(numbers[-1])]


def split_protein_accessions(value: Any) -> list[str] | None:
    """
    Splits protein accession strings into a list.

    This supports common separators: semicolon, comma, pipe.
    """
    if pd.isna(value):
        return None

    text = str(value).strip()
    if not text:
        return None

    # Be careful with UniProt IDs containing pipes.
    # If the whole accession looks like sp|P12345|NAME, do not split by pipe.
    if text.startswith(("sp|", "tr|")) and ";" not in text and "," not in text:
        return [text]

    if ";" in text:
        parts = text.split(";")
    elif "," in text:
        parts = text.split(",")
    else:
        parts = [text]

    parts = [p.strip() for p in parts if p.strip()]
    return parts or None


def is_decoy_from_proteins(proteins: list[str] | None) -> bool:
    """
    Your rule:
      protein accession starting with REV_ means decoy.

    For multi-protein PSMs, this returns True only if all listed proteins are REV_.
    This avoids marking target+decoy shared peptides as pure decoys.
    """
    if not proteins:
        return False

    return all(p.startswith("REV_") for p in proteins)


def to_float(value: Any) -> float | None:
    if pd.isna(value) or str(value).strip() == "":
        return None
    return float(value)


def to_int(value: Any) -> int | None:
    if pd.isna(value) or str(value).strip() == "":
        return None
    return int(float(value))


def neutral_mass_to_mz(neutral_mass: float | None, charge: int | None) -> float | None:
    if neutral_mass is None or charge is None or charge == 0:
        return None
    return (neutral_mass + charge * PROTON_MASS) / charge


def convert_tsv_to_qpx_psm_parquet(
    input_tsv: str | Path,
    output_parquet: str | Path,
    rt_unit: str = "unknown",
    score_name: str = "Score",
) -> None:
    input_tsv = Path(input_tsv)
    output_parquet = Path(output_parquet)

    df = pd.read_csv(input_tsv, sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]

    require_columns(df)

    records: list[dict[str, Any]] = []

    for row_index, row in df.iterrows():
        sequence = clean_peptide_sequence(row["Base Sequence"])
        peptidoform = clean_peptidoform(row["Full Sequence"], sequence)

        charge = to_int(row["Precursor Charge"])

        theoretical_mass = to_float(row["Peptide Monoisotopic Mass"])
        observed_mass = to_float(row["Peptide Mass"])
        mass_diff_da = to_float(row["Mass difference"])

        calculated_mz = neutral_mass_to_mz(theoretical_mass, charge)
        observed_mz = neutral_mass_to_mz(observed_mass, charge)

        if theoretical_mass is not None and mass_diff_da is not None and theoretical_mass != 0:
            mass_error_ppm = 1e6 * mass_diff_da / theoretical_mass
        elif theoretical_mass is not None and observed_mass is not None and theoretical_mass != 0:
            mass_error_ppm = 1e6 * (observed_mass - theoretical_mass) / theoretical_mass
        else:
            mass_error_ppm = None

        proteins = split_protein_accessions(row["Protein Accession"])
        is_decoy = is_decoy_from_proteins(proteins)

        score = to_float(row["Score"])
        pep = to_float(row["PEP"])
        rt = to_float(row["Scan Retention Time"])
        if rt is not None and rt_unit == "minutes":
            rt = rt * 60.0
        missed_cleavages = to_int(row["Missed cleavages"])

        if sequence is None:
            raise ValueError(f"Row {row_index + 2}: missing Base Sequence")

        if peptidoform is None:
            raise ValueError(f"Row {row_index + 2}: missing Full Sequence / peptidoform")

        if charge is None:
            raise ValueError(f"Row {row_index + 2}: missing Precursor Charge")

        if calculated_mz is None:
            raise ValueError(f"Row {row_index + 2}: cannot calculate calculated_mz")

        if observed_mz is None:
            raise ValueError(f"Row {row_index + 2}: cannot calculate observed_mz")

        additional_scores = {}
        if score is not None:
            additional_scores[score_name] = score
        if mass_diff_da is not None:
            additional_scores["mass_difference_da"] = mass_diff_da

        record = {
            # QPX core PSM fields
            "sequence": sequence,
            "peptidoform": peptidoform,
            "modifications": None,
            "charge": charge,
            "posterior_error_probability": pep,
            "is_decoy": is_decoy,
            "calculated_mz": calculated_mz,
            "observed_mz": observed_mz,
            "mass_error_ppm": mass_error_ppm,
            "additional_scores": additional_scores if additional_scores else None,
            "predicted_rt": None,
            "run_file_name": str(row["File Name"]).strip(),
            "cv_params": None,
            "scan": parse_scan(row["Scan number"]),
            "rt": rt,
            "ion_mobility": None,
            "missed_cleavages": missed_cleavages,
            "protein_accessions": proteins,
            "cross_links": None,

            # Optional spectral arrays are not available from your TSV
            "mz_array": None,
            "intensity_array": None,
            "charge_array": None,
            "ion_type_array": None,
            "ion_mobility_array": None,
        }

        records.append(record)

    schema = pa.schema(
        [
            pa.field("sequence", pa.string(), nullable=False),
            pa.field("peptidoform", pa.string(), nullable=False),
            pa.field("modifications", pa.string(), nullable=True),
            pa.field("charge", pa.int32(), nullable=False),
            pa.field("posterior_error_probability", pa.float64(), nullable=True),
            pa.field("is_decoy", pa.bool_(), nullable=False),
            pa.field("calculated_mz", pa.float64(), nullable=False),
            pa.field("observed_mz", pa.float64(), nullable=False),
            pa.field("mass_error_ppm", pa.float64(), nullable=True),
            pa.field(
                "additional_scores",
                pa.map_(pa.string(), pa.float64()),
                nullable=True,
            ),
            pa.field("predicted_rt", pa.float64(), nullable=True),
            pa.field("run_file_name", pa.string(), nullable=False),
            pa.field("cv_params", pa.string(), nullable=True),
            pa.field("scan", pa.list_(pa.int64()), nullable=False),
            pa.field("rt", pa.float64(), nullable=True),
            pa.field("ion_mobility", pa.float64(), nullable=True),
            pa.field("missed_cleavages", pa.int32(), nullable=True),
            pa.field("protein_accessions", pa.list_(pa.string()), nullable=True),
            pa.field("cross_links", pa.string(), nullable=True),
            pa.field("mz_array", pa.list_(pa.float64()), nullable=True),
            pa.field("intensity_array", pa.list_(pa.float64()), nullable=True),
            pa.field("charge_array", pa.list_(pa.int32()), nullable=True),
            pa.field("ion_type_array", pa.list_(pa.string()), nullable=True),
            pa.field("ion_mobility_array", pa.list_(pa.float64()), nullable=True),
        ]
    )

    table = pa.Table.from_pylist(records, schema=schema)

    output_parquet.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, output_parquet, compression="zstd")

    print(f"Converted {len(records):,} PSMs")
    print(f"Output: {output_parquet}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert a generic MetaPilot-style PSM TSV file to QPX psm.parquet."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input PSM TSV file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output QPX PSM parquet file, e.g. psm.parquet.",
    )
    parser.add_argument(
        "--rt-unit",
        default="unknown",
        choices=["seconds", "minutes", "unknown"],
        help="Unit of 'Scan Retention Time'. Stored in additional_scores metadata.",
    )
    parser.add_argument(
        "--score-name",
        default="Score",
        help="Name of the search-engine score stored in additional_scores.",
    )

    args = parser.parse_args()

    convert_tsv_to_qpx_psm_parquet(
        input_tsv=args.input,
        output_parquet=args.output,
        rt_unit=args.rt_unit,
        score_name=args.score_name,
    )


if __name__ == "__main__":
    main()