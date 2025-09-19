#!/usr/bin/env python3
"""Compute read pair contamination and relatedness metrics from STAR gene counts."""

from __future__ import annotations

import itertools
import logging
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd


def setup_logger(log_path: str | None) -> logging.Logger:
    """Configure a logger that writes to the provided log file."""

    logger = logging.getLogger("relatedness")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    handler: logging.Handler
    if log_path:
        Path(log_path).parent.mkdir(parents=True, exist_ok=True)
        handler = logging.FileHandler(log_path)
    else:
        handler = logging.StreamHandler()

    handler.setFormatter(
        logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    )
    logger.addHandler(handler)
    logger.propagate = False
    return logger


def get_column(strandedness: str | None) -> int:
    """Map STAR strandedness descriptors to ReadsPerGene column indices."""

    if strandedness in {None, "", "none"}:
        return 1
    if strandedness == "yes":
        return 2
    if strandedness == "reverse":
        return 3
    raise ValueError(
        (
            "'strandedness' column should be empty or have the value "
            "'none', 'yes' or 'reverse', instead has the value {}"
        ).format(repr(strandedness))
    )


def load_counts(
    count_paths: Sequence[str],
    sample_names: Sequence[str],
    unit_names: Sequence[str],
    strands: Sequence[str],
    logger: logging.Logger,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load STAR ReadsPerGene outputs into a count matrix and metadata table."""

    if not (
        len(count_paths)
        == len(sample_names)
        == len(unit_names)
        == len(strands)
    ):
        raise ValueError(
            "Expected the number of read count files, samples, units and strands "
            f"to match, got {len(count_paths)} files, {len(sample_names)} samples, "
            f"{len(unit_names)} units and {len(strands)} strand entries."
        )

    series_list: List[pd.Series] = []
    metadata_records = []

    for path, sample, unit, strand in zip(
        count_paths, sample_names, unit_names, strands
    ):
        column = get_column(strand)
        table = pd.read_table(
            path, header=None, index_col=0, usecols=[0, column], skiprows=4
        )
        series = table.iloc[:, 0].astype(float)
        key = f"{sample}::{unit}"
        series.name = key
        series_list.append(series)
        metadata_records.append(
            {"key": key, "sample": sample, "unit": unit, "path": str(path)}
        )

    counts = pd.concat(series_list, axis=1).fillna(0.0)
    counts.columns.name = "sample_unit"

    metadata = pd.DataFrame.from_records(metadata_records).set_index("key")

    logger.info(
        "Loaded %d read-pair count tables spanning %d genes.",
        len(series_list),
        counts.shape[0],
    )
    return counts, metadata


def attach_patient_metadata(
    metadata: pd.DataFrame, samples_path: str, logger: logging.Logger
) -> pd.DataFrame:
    """Attach patient identifiers to each sample-unit entry."""

    samples_df = pd.read_csv(samples_path, sep="\t", dtype=str)

    if "patient" not in samples_df.columns:
        raise ValueError(
            "The samples sheet must contain a 'patient' column to assess relatedness."
        )

    sample_patients = samples_df.set_index("sample_name")["patient"].to_dict()

    metadata = metadata.copy()
    patients = metadata["sample"].map(sample_patients)
    patients = patients.fillna("")
    missing_mask = patients.eq("") | patients.isna()
    if missing_mask.any():
        logger.warning(
            "Assigning sample identifiers as patient labels for %d entries without explicit patient information.",
            int(missing_mask.sum()),
        )
        patients.loc[missing_mask] = metadata.loc[missing_mask, "sample"]

    metadata["patient"] = patients
    logger.info(
        "Detected %d unique patients across %d read pairs.",
        metadata["patient"].nunique(),
        metadata.shape[0],
    )
    return metadata


def compute_correlations(counts: pd.DataFrame, logger: logging.Logger) -> pd.DataFrame:
    """Compute a Pearson correlation matrix on log-transformed counts."""

    if counts.empty:
        logger.warning("No counts were provided; returning an empty correlation matrix.")
        return pd.DataFrame()

    if counts.shape[1] == 1:
        logger.info(
            "Only one read pair available; correlation matrix will contain a single entry."
        )
        return pd.DataFrame([[1.0]], index=counts.columns, columns=counts.columns)

    filtered = counts.loc[(counts > 0).any(axis=1)]
    if filtered.empty:
        filtered = counts

    log_counts = np.log1p(filtered)
    corr = log_counts.corr(method="pearson")
    for column in corr.columns:
        corr.loc[column, column] = 1.0

    logger.info(
        "Computed correlation matrix for %d read pairs (shape %s).",
        corr.shape[0],
        corr.shape,
    )
    return corr


def _best_match(
    correlations: pd.Series, valid_keys: Iterable[str]
) -> Tuple[str | None, float]:
    """Return the key with the highest correlation among *valid_keys*."""

    key_list = list(valid_keys)
    if not key_list:
        return None, float("nan")

    if correlations.empty:
        return None, float("nan")

    subset = correlations.loc[key_list]
    subset = subset.dropna()
    if subset.empty:
        return None, float("nan")

    best_key = subset.idxmax()
    return best_key, float(subset.loc[best_key])


def summarize_readpair_contamination(
    metadata: pd.DataFrame,
    corr: pd.DataFrame,
    contamination_threshold: float,
    same_patient_threshold: float,
) -> pd.DataFrame:
    """Summarize per-read-pair relationships and potential contamination."""

    rows = []
    for key, row in metadata.iterrows():
        correlations = corr.loc[key] if not corr.empty else pd.Series(dtype=float)
        correlations = correlations.drop(index=key, errors="ignore")

        best_key, best_corr = _best_match(correlations, correlations.index)

        same_patient_keys = metadata.index[metadata["patient"] == row["patient"]]
        same_patient_keys = [k for k in same_patient_keys if k != key]
        best_same_key, best_same_corr = _best_match(correlations, same_patient_keys)

        cross_patient_keys = metadata.index[metadata["patient"] != row["patient"]]
        best_cross_key, best_cross_corr = _best_match(correlations, cross_patient_keys)

        cross_flag = (
            not np.isnan(best_cross_corr)
            and best_cross_corr >= contamination_threshold
        )

        if best_same_key is None:
            same_flag = False
            same_note = "no_additional_samples_for_patient"
        else:
            same_flag = bool(
                np.isnan(best_same_corr)
                or best_same_corr < same_patient_threshold
            )
            same_note = "below_threshold" if same_flag else ""

        rows.append(
            {
                "sample": row["sample"],
                "unit": row["unit"],
                "patient": row["patient"],
                "best_match_sample": metadata.loc[best_key, "sample"]
                if best_key
                else "",
                "best_match_unit": metadata.loc[best_key, "unit"]
                if best_key
                else "",
                "best_match_patient": metadata.loc[best_key, "patient"]
                if best_key
                else "",
                "best_match_correlation": best_corr,
                "best_same_patient_sample": metadata.loc[best_same_key, "sample"]
                if best_same_key
                else "",
                "best_same_patient_unit": metadata.loc[best_same_key, "unit"]
                if best_same_key
                else "",
                "best_same_patient_correlation": best_same_corr,
                "best_cross_patient_sample": metadata.loc[best_cross_key, "sample"]
                if best_cross_key
                else "",
                "best_cross_patient_unit": metadata.loc[best_cross_key, "unit"]
                if best_cross_key
                else "",
                "best_cross_patient_patient": metadata.loc[best_cross_key, "patient"]
                if best_cross_key
                else "",
                "best_cross_patient_correlation": best_cross_corr,
                "cross_patient_flag": bool(cross_flag),
                "low_same_patient_flag": bool(same_flag),
                "same_patient_note": same_note,
            }
        )

    return pd.DataFrame(rows)


def summarize_patient_relationships(
    metadata: pd.DataFrame,
    corr: pd.DataFrame,
    same_patient_threshold: float,
) -> pd.DataFrame:
    """Summarize correlations among samples that share the same patient."""

    records = []
    for patient, indices in metadata.groupby("patient").groups.items():
        key_list = list(indices)
        if len(key_list) < 2:
            key = key_list[0]
            records.append(
                {
                    "patient": patient,
                    "patient_sample_count": len(key_list),
                    "sample_a": metadata.loc[key, "sample"],
                    "unit_a": metadata.loc[key, "unit"],
                    "sample_b": "",
                    "unit_b": "",
                    "correlation": float("nan"),
                    "flag_low": False,
                    "note": "only_one_sample_for_patient",
                }
            )
            continue

        for key_a, key_b in itertools.combinations(key_list, 2):
            correlation = (
                corr.loc[key_a, key_b] if not corr.empty else float("nan")
            )
            flag_low = bool(
                np.isnan(correlation) or correlation < same_patient_threshold
            )
            note = "correlation_nan" if np.isnan(correlation) else ""
            if flag_low and note == "":
                note = "below_threshold"

            records.append(
                {
                    "patient": patient,
                    "patient_sample_count": len(key_list),
                    "sample_a": metadata.loc[key_a, "sample"],
                    "unit_a": metadata.loc[key_a, "unit"],
                    "sample_b": metadata.loc[key_b, "sample"],
                    "unit_b": metadata.loc[key_b, "unit"],
                    "correlation": correlation,
                    "flag_low": flag_low,
                    "note": note,
                }
            )

    return pd.DataFrame.from_records(records)


def summarize_cross_patient_matches(
    metadata: pd.DataFrame,
    corr: pd.DataFrame,
    contamination_threshold: float,
) -> pd.DataFrame:
    """Report cross-patient pairs whose correlation exceeds the threshold."""

    rows = []
    for key_a, key_b in itertools.combinations(metadata.index, 2):
        if metadata.loc[key_a, "patient"] == metadata.loc[key_b, "patient"]:
            continue

        correlation = (
            corr.loc[key_a, key_b] if not corr.empty else float("nan")
        )
        if np.isnan(correlation) or correlation < contamination_threshold:
            continue

        rows.append(
            {
                "sample_a": metadata.loc[key_a, "sample"],
                "unit_a": metadata.loc[key_a, "unit"],
                "patient_a": metadata.loc[key_a, "patient"],
                "sample_b": metadata.loc[key_b, "sample"],
                "unit_b": metadata.loc[key_b, "unit"],
                "patient_b": metadata.loc[key_b, "patient"],
                "correlation": correlation,
                "flag": True,
            }
        )

    if rows:
        df = pd.DataFrame(rows).sort_values(
            by="correlation", ascending=False
        ).reset_index(drop=True)
    else:
        df = pd.DataFrame(
            columns=
            [
                "sample_a",
                "unit_a",
                "patient_a",
                "sample_b",
                "unit_b",
                "patient_b",
                "correlation",
                "flag",
            ]
        )

    return df


def main() -> None:
    logger = setup_logger(snakemake.log[0] if snakemake.log else None)

    count_paths = [str(path) for path in snakemake.input.counts]
    sample_names = list(snakemake.params.sample_names)
    unit_names = list(snakemake.params.unit_names)
    strands = list(snakemake.params.strands)

    counts, metadata = load_counts(
        count_paths=count_paths,
        sample_names=sample_names,
        unit_names=unit_names,
        strands=strands,
        logger=logger,
    )

    metadata = attach_patient_metadata(
        metadata=metadata,
        samples_path=str(snakemake.input.samples),
        logger=logger,
    )

    corr = compute_correlations(counts=counts, logger=logger)

    mode = snakemake.params.mode
    contamination_threshold = float(snakemake.params.contamination_threshold)
    same_patient_threshold = float(snakemake.params.same_patient_threshold)

    logger.info(
        "Running relatedness analysis in '%s' mode (contamination threshold %.3f, same patient threshold %.3f).",
        mode,
        contamination_threshold,
        same_patient_threshold,
    )

    if mode == "readpair":
        result = summarize_readpair_contamination(
            metadata=metadata,
            corr=corr,
            contamination_threshold=contamination_threshold,
            same_patient_threshold=same_patient_threshold,
        )
    elif mode == "patient":
        result = summarize_patient_relationships(
            metadata=metadata,
            corr=corr,
            same_patient_threshold=same_patient_threshold,
        )
    elif mode == "cross":
        result = summarize_cross_patient_matches(
            metadata=metadata,
            corr=corr,
            contamination_threshold=contamination_threshold,
        )
    else:
        raise ValueError(f"Unsupported mode '{mode}'.")

    Path(snakemake.output[0]).parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(snakemake.output[0], sep="\t", index=False)

    logger.info(
        "Wrote %d rows to %s.", len(result.index), snakemake.output[0]
    )


if __name__ == "__main__":
    main()
