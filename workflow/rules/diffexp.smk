from pathlib import Path
import re


_FRACTION_PATTERN = re.compile(r":\s*([0-9]*\.?[0-9]+)")
_FORWARD_KEYS = ("1++,1--,2+-,2-+", "++,--")
_REVERSE_KEYS = ("1+-,1-+,2++,2--", "+-,-+")
_STRANDEDNESS_THRESHOLD = 0.6


def _extract_fraction(line):
    match = _FRACTION_PATTERN.search(line)
    if not match:
        return None
    try:
        return float(match.group(1))
    except ValueError:
        return None


def _infer_strand_from_file(path, default):
    try:
        lines = Path(path).read_text().splitlines()
    except OSError:
        return default

    forward = None
    reverse = None
    for line in lines:
        if any(key in line for key in _FORWARD_KEYS):
            value = _extract_fraction(line)
            if value is not None:
                forward = value
        elif any(key in line for key in _REVERSE_KEYS):
            value = _extract_fraction(line)
            if value is not None:
                reverse = value

    if forward is None and reverse is None:
        return default

    forward = forward if forward is not None else 0.0
    reverse = reverse if reverse is not None else 0.0

    if forward >= _STRANDEDNESS_THRESHOLD and forward >= reverse:
        return "yes"
    if reverse >= _STRANDEDNESS_THRESHOLD and reverse > forward:
        return "reverse"

    return default


def _get_inferred_strands(infer_files):
    defaults = get_strandedness(units)
    if not infer_files:
        return defaults

    strands = [
        _infer_strand_from_file(path, default)
        for path, default in zip(infer_files, defaults)
    ]

    if len(strands) < len(defaults):
        strands.extend(defaults[len(strands):])
    elif len(infer_files) > len(defaults):
        strands.extend(
            _infer_strand_from_file(path, "none")
            for path in infer_files[len(defaults) :]
        )

    return strands


def _nz(val):
    return "none" if val is None or str(val).strip() == "" else str(val)


rule count_matrix:
    input:
        reads = expand(
            "results/star/{unit.sample_name}_{unit.unit_name}/ReadsPerGene.out.tab",
            unit=units.itertuples(),
        ),
        # OPTIONAL: provide per-sample infer_experiment outputs when available
        infer = expand(
            "results/qc/rseqc/{unit.sample_name}_{unit.unit_name}.infer_experiment.txt",
            unit=units.itertuples(),
        )
    output:
        "results/counts/all.tsv"
    log:
        "logs/count-matrix.log",
    benchmark:
        "logs/benchmarks/count_matrix.bench.tsv",
    params:
        samples = ",".join(units["sample_name"].tolist()),
        strands = lambda wildcards, input: ",".join(
            _nz(s) for s in _get_inferred_strands(list(input.infer))
        )
    conda:
        "../envs/pandas.yaml"
    threads: 1
    shell:
        """
        python workflow/scripts/count-matrix.py \
            --output {output} \
            --samples "{params.samples}" \
            --strands "{params.strands}" \
            {input.reads} \
            > {log} 2>&1
        """


rule gene_2_symbol:
    input:
        counts="{prefix}.tsv",
    output:
        symbol="{prefix}.symbol.tsv",
        nodata="{prefix}.symbol.tsv.NODATA",
    params:
        species=get_bioc_species_name(),
    log:
        "logs/gene2symbol/{prefix}.log",
    benchmark:
        "logs/benchmarks/gene_2_symbol/{prefix}.bench.tsv",
    conda:
        "../envs/biomart.yaml"
    shell:
        """

        touch {output.symbol};
        touch {output.nodata};
        
        echo "Annotating gene symbols for {input.counts}." > {log};
        Rscript workflow/scripts/gene2symbol.R \
            --counts {input.counts} \
            --output {output.symbol} \
            --species {params.species} \
            >> {log} 2>&1;
        
        """


rule deseq2_init:
    input:
        counts="results/counts/all.tsv",
    output:
        "results/deseq2/all.rds",
        "results/deseq2/normcounts.tsv",
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/init.log",
    benchmark:
        "logs/benchmarks/deseq2_init.bench.tsv",
    threads: get_deseq2_threads()
    shell:
        """
        echo 'A' >> {log};
        touch {output} >> {log};
        echo 'B' >> {log}
        """


localrules: pca, deseq2
rule pca:
    input:
        "results/deseq2/all.rds",
    output:
        report("results/pca.{variable}.svg", "../report/pca.rst"),
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/pca.{variable}.log",
    benchmark:
        "logs/benchmarks/pca/{variable}.bench.tsv",
    shell:
        "touch {output}"


rule deseq2:
    input:
        "results/deseq2/all.rds",
    output:
        table=report("results/diffexp/{contrast}.diffexp.tsv", "../report/diffexp.rst"),
        ma_plot=report("results/diffexp/{contrast}.ma-plot.svg", "../report/ma.rst"),
    params:
        contrast=get_contrast,
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.log",
    benchmark:
        "logs/benchmarks/deseq2/{contrast}.bench.tsv",
    threads: get_deseq2_threads()
    shell:
        "touch {output}"

