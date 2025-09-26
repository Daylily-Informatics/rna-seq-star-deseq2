import glob
import re

import pandas as pd
from snakemake.utils import validate

validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)


def _slugify_build(build: str) -> str:
    slug = re.sub(r"[^0-9a-zA-Z]+", "", str(build)).lower()
    return slug or "build"


GENOME_BUILD_SLUG = _slugify_build(config["ref"].get("build", ""))
ANALYSIS_ROOT = f"results/rna/{GENOME_BUILD_SLUG}"


def sample_path(sample: str, *parts: str) -> str:
    return "/".join([ANALYSIS_ROOT, sample, *parts])


def cohort_path(*parts: str) -> str:
    return "/".join([ANALYSIS_ROOT, "cohort", *parts])


def trimmed_fastq_path(sample: str, unit: str, read: str) -> str:
    labels = {"fq1": "R1", "fq2": "R2", "single": "single"}
    if read not in labels:
        raise ValueError(f"Unsupported read label '{read}' for trimmed fastq path")
    return sample_path(sample, "fastq", "trimmed", f"{sample}.{unit}.trimmed.{labels[read]}.fastq.gz")


def trimmed_qc_path(sample: str, unit: str, suffix: str) -> str:
    return sample_path(sample, "fastq", "trimmed", f"{sample}.{unit}.trimmed.{suffix}")


def star_bam_path(sample: str, unit: str) -> str:
    return sample_path(sample, "align", "star", f"{sample}.{unit}.star.sortedByCoord.bam")


def star_counts_path(sample: str, unit: str) -> str:
    return sample_path(sample, "align", "star", f"{sample}.{unit}.star.ReadsPerGene.out.tab")


def fastqc_output_path(sample: str, unit: str, read: str, ext: str) -> str:
    return sample_path(sample, "qc", "fastqc", f"{sample}.{unit}.{read}.fastqc.{ext}")


def rseqc_output_path(sample: str, unit: str, suffix: str) -> str:
    return sample_path(sample, "qc", "rseqc", f"{sample}.{unit}.{suffix}")


def get_final_output():
    final_output = expand(
        cohort_path("deseq2", "diffexp", "{contrast}.diffexp.symbol.tsv"),
        contrast=config["diffexp"]["contrasts"],
    )
    final_output.append(cohort_path("counts", "rnaseq.counts.symbol.tsv"))
    final_output.append(cohort_path("reports", "multiqc", "multiqc_report.html"))

    if config["pca"]["activate"]:
        pca_variables = list(config["diffexp"]["variables_of_interest"])
        if config["diffexp"]["batch_effects"]:
            pca_variables.extend(config["diffexp"]["batch_effects"])
        if config["pca"]["labels"]:
            pca_variables.extend(config["pca"]["labels"])
        final_output.extend(
            expand(
                cohort_path("deseq2", "pca", "{variable}.pca.svg"),
                variable=pca_variables,
            )
        )
    return final_output


validate(samples, schema="../schemas/samples.schema.yaml")

#units = (
#    pd.read_csv(config["units"], sep="\t", dtype={"sample_name": str, "unit_name": str})
#    .set_index(["sample_name", "unit_name"], drop=False)
#    .sort_index()
#)
units = (
    pd.read_csv(
        config["units"],
        sep="\t",
        dtype=str,               # ← all cols strings
        keep_default_na=False,   # ← "" stays "", no NaN
        na_filter=False,
    )
    .applymap(lambda x: x.strip() if isinstance(x, str) else x)
    .set_index(["sample_name", "unit_name"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")


wildcard_constraints:
    sample="|".join(samples["sample_name"]),
    unit="|".join(units["unit_name"]),


def get_cutadapt_input(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]

    if pd.isna(unit["fq1"]):
        # SRA sample (always paired-end for now)
        accession = unit["sra"]
        return expand("sra/{accession}_{read}.fastq", accession=accession, read=[1, 2])

    if unit["fq1"].endswith("gz"):
        ending = ".gz"
    else:
        ending = ""

    if pd.isna(unit["fq2"]):
        # single end local sample
        return "pipe/cutadapt/{S}/{U}.fq1.fastq{E}".format(
            S=unit.sample_name, U=unit.unit_name, E=ending
        )
    else:
        # paired end local sample
        return expand(
            "pipe/cutadapt/{S}/{U}.{{read}}.fastq{E}".format(
                S=unit.sample_name, U=unit.unit_name, E=ending
            ),
            read=["fq1", "fq2"],
        )


def get_cutadapt_pipe_input(wildcards):
    files = list(
        sorted(glob.glob(units.loc[wildcards.sample].loc[wildcards.unit, wildcards.fq]))
    )
    assert len(files) > 0
    return files


def is_paired_end(sample):
    sample_units = units.loc[sample]
    fq2_null = sample_units["fq2"].isnull()
    sra_null = sample_units["sra"].isnull()
    paired = ~fq2_null | ~sra_null
    all_paired = paired.all()
    all_single = (~paired).all()
    assert (
        all_single or all_paired
    ), "invalid units for sample {}, must be all paired end or all single end".format(
        sample
    )
    return all_paired


def get_fq(wildcards):
    if config["trimming"]["activate"]:
        # activated trimming, use trimmed data
        if is_paired_end(wildcards.sample):
            # paired-end sample
            return dict(
                zip(
                    ["fq1", "fq2"],
                    [
                        trimmed_fastq_path(wildcards.sample, wildcards.unit, "fq1"),
                        trimmed_fastq_path(wildcards.sample, wildcards.unit, "fq2"),
                    ],
                )
            )
        # single end sample
        return {
            "fq1": trimmed_fastq_path(wildcards.sample, wildcards.unit, "single")
        }
    else:
        # no trimming, use raw reads
        u = units.loc[(wildcards.sample, wildcards.unit)]
        if pd.isna(u["fq1"]):
            # SRA sample (always paired-end for now)
            accession = u["sra"]
            return dict(
                zip(
                    ["fq1", "fq2"],
                    expand(
                        "sra/{accession}_{group}.fastq",
                        accession=accession,
                        group=["R1", "R2"],
                    ),
                )
            )
        if not is_paired_end(wildcards.sample):
            return {"fq1": f"{u.fq1}"}
        else:
            return {"fq1": f"{u.fq1}", "fq2": f"{u.fq2}"}


def get_fastqc_fastq(wildcards):
    fastqs = get_fq(wildcards)
    read_map = {"R1": "fq1", "R2": "fq2"}
    try:
        key = read_map[wildcards.read]
    except KeyError:
        raise ValueError(
            "Invalid read value '{read}' for sample {sample} unit {unit}".format(
                read=wildcards.read, sample=wildcards.sample, unit=wildcards.unit
            )
        )

    fastq = fastqs.get(key)
    if fastq is None:
        raise ValueError(
            "Read {read} not available for sample {sample} unit {unit}".format(
                read=wildcards.read, sample=wildcards.sample, unit=wildcards.unit
            )
        )

    return fastq


def get_strandedness(units):
    if "strandedness" in units.columns:
        return units["strandedness"].tolist()
    else:
        strand_list = ["none"]
        return strand_list * units.shape[0]


def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6


def is_activated(xpath):
    c = config
    for entry in xpath.split("/"):
        c = c.get(entry, {})
    return bool(c.get("activate", False))


def get_bioc_species_name():
    first_letter = config["ref"]["species"][0]
    subspecies = config["ref"]["species"].split("_")[1]
    return first_letter + subspecies


def get_contrast(wildcards):
    return config["diffexp"]["contrasts"][wildcards.contrast]
