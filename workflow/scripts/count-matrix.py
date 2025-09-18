import json
import sys
from pathlib import Path

import pandas as pd


sys.stderr = open(snakemake.log[0], "w")


STAR_COLUMNS = {"none": 1, "yes": 2, "reverse": 3}


def infer_star_column(infer_path: Path, fallback: str, sample: str) -> int:
    fallback = fallback if isinstance(fallback, str) else "none"
    column = STAR_COLUMNS.get(fallback, 1)

    if not infer_path.exists():
        print(
            f"[count-matrix] sample={sample} infer_experiment missing; fallback to {fallback} (col {column})"
        )
        return column

    sf = sr = None
    with infer_path.open() as handle:
        for line in handle:
            if line.startswith("1++,1--,2+-,2-+"):
                try:
                    sf = float(line.split(":", 1)[1].strip())
                except ValueError:
                    pass
            elif line.startswith("1+-,1-+,2++,2--"):
                try:
                    sr = float(line.split(":", 1)[1].strip())
                except ValueError:
                    pass

    if sf is not None and sf > 0.6:
        column = 2  # STAR column 3
    elif sr is not None and sr > 0.6:
        column = 3  # STAR column 4
    else:
        column = 1  # STAR column 2

    print(
        json.dumps(
            {
                "sample": sample,
                "infer": str(infer_path),
                "sf": sf,
                "sr": sr,
                "fallback": fallback,
                "column": column,
            }
        )
    )
    return column


reads_files = list(snakemake.input.get("reads", []))
infer_files = list(snakemake.input.get("infer", []))
samples = list(snakemake.params.samples)
fallbacks = list(snakemake.params.fallback_strand)

if not (len(reads_files) == len(samples) == len(fallbacks)):
    raise ValueError("Input/parameter lengths are inconsistent for count matrix generation")

if len(infer_files) not in (0, len(reads_files)):
    raise ValueError("Mismatch between gene count files and infer_experiment outputs")

counts = []
for read_file, infer_file, sample, fallback in zip(
    reads_files,
    infer_files,
    samples,
    fallbacks,
):
    column = infer_star_column(Path(infer_file), fallback, sample)
    table = pd.read_table(
        read_file,
        index_col=0,
        usecols=[0, column],
        header=None,
        skiprows=4,
    )
    table.columns = [sample]
    counts.append(table)

if not counts:
    raise ValueError("No STAR gene count files detected")

matrix = pd.concat(counts, axis=1)
matrix.index.name = "gene"
matrix = matrix.groupby(matrix.columns, axis=1, sort=False).sum()
matrix.to_csv(snakemake.output[0], sep="\t")
