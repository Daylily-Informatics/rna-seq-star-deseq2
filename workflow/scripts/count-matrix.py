#!/usr/bin/env python3
"""Build a count matrix from STAR gene counts outputs."""

import argparse
from typing import List

import pandas as pd


def parse_comma_separated(value: str, name: str) -> List[str]:
    items = [item.strip() for item in value.split(",") if item.strip()]
    if not items:
        raise argparse.ArgumentTypeError(f"{name} list cannot be empty")
    return items


def get_column(strandedness: str) -> int:
    if strandedness in {"", "none", None}:
        return 1  # non stranded protocol
    if strandedness == "yes":
        return 2  # 3rd column
    if strandedness == "reverse":
        return 3  # 4th column, usually for Illumina truseq
    raise ValueError(
        (
            "'strandedness' column should be empty or have the value "
            "'none', 'yes' or 'reverse', instead has the value {}"
        ).format(repr(strandedness))
    )


def build_matrix(inputs: List[str], strands: List[str], samples: List[str]) -> pd.DataFrame:
    counts = []
    for path, strand, sample in zip(inputs, strands, samples):
        table = pd.read_table(
            path, index_col=0, usecols=[0, get_column(strand)], header=None, skiprows=4
        )
        table.columns = [sample]
        counts.append(table)

    matrix = pd.concat(counts, axis=1)
    matrix.index.name = "gene"
    # collapse technical replicates
    matrix = matrix.groupby(matrix.columns, axis=1, sort=False).sum()
    return matrix


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--samples",
        required=True,
        help="Comma-separated list of sample names matching the input files.",
    )
    parser.add_argument(
        "--strands",
        required=True,
        help="Comma-separated list of strandedness values matching the input files.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to write the combined count matrix to.",
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="STAR ReadsPerGene.out.tab files to combine.",
    )

    args = parser.parse_args()

    samples = parse_comma_separated(args.samples, "Sample")
    strands = parse_comma_separated(args.strands, "Strand")
    inputs = args.inputs

    if not inputs:
        parser.error("At least one input file must be provided.")

    if not (len(inputs) == len(samples) == len(strands)):
        parser.error(
            "The number of inputs, samples, and strands must match: "
            f"got {len(inputs)} inputs, {len(samples)} samples and {len(strands)} strands."
        )

    matrix = build_matrix(inputs, strands, samples)
    matrix.to_csv(args.output, sep="\t")


if __name__ == "__main__":
    main()
