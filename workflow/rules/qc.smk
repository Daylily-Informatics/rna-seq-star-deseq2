## RSEQC

FASTQC_WILDCARDS = []
for unit in units.itertuples():
    reads = ("R1", "R2") if is_paired_end(unit.sample_name) else ("R1",)
    for read in reads:
        FASTQC_WILDCARDS.append(
            {"sample": unit.sample_name, "unit": unit.unit_name, "read": read}
        )

FASTQC_ZIP_OUTPUTS = [
    fastqc_output_path(entry["sample"], entry["unit"], entry["read"], "zip")
    for entry in FASTQC_WILDCARDS
]

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

def get_multiqc_inputs(_wildcards=None):
    inputs = list(FASTQC_ZIP_OUTPUTS)
    for unit in units.itertuples():
        sample = unit.sample_name
        unit_name = unit.unit_name
        inputs.append(star_bam_path(sample, unit_name))
        inputs.append(rseqc_output_path(sample, unit_name, "junctionanno.junction.bed"))
        inputs.append(
            rseqc_output_path(
                sample, unit_name, "junctionsat.junctionSaturation_plot.pdf"
            )
        )
        inputs.append(rseqc_output_path(sample, unit_name, "infer_experiment.txt"))
        inputs.append(rseqc_output_path(sample, unit_name, "stats.txt"))
        inputs.append(
            rseqc_output_path(
                sample, unit_name, "inner_distance_freq.inner_distance.txt"
            )
        )
        inputs.append(rseqc_output_path(sample, unit_name, "readdistribution.txt"))
        inputs.append(
            rseqc_output_path(sample, unit_name, "readdup.DupRate_plot.pdf")
        )
        inputs.append(rseqc_output_path(sample, unit_name, "readgc.GC_plot.pdf"))
        inputs.append(
            f"logs/rseqc/rseqc_junction_annotation/{sample}_{unit_name}.log"
        )
    return inputs


rule fastqc:
    input:
        fastq=lambda wildcards: get_fastqc_fastq(wildcards),
    output:
        html=fastqc_output_path("{sample}", "{unit}", "{read}", "html"),
        zip=fastqc_output_path("{sample}", "{unit}", "{read}", "zip"),
    threads: 4
    log:
        "logs/fastqc/{sample}_{unit}_{read}.log",
    benchmark:
        "logs/benchmarks/fastqc/{sample}_{unit}_{read}.bench.tsv",
    params:
        extra="",
    wrapper:
        "v3.5.3/bio/fastqc"


rule rseqc_gtf2bed:
    input:
        "resources/genome.gtf",
    output:
        bed=cohort_path("qc", "rseqc", "annotation.bed"),
        db=temp(cohort_path("qc", "rseqc", "annotation.db")),
    threads: 32
    log:
        "logs/rseqc_gtf2bed.log",
    benchmark:
        "logs/benchmarks/rseqc_gtf2bed.bench.tsv",
    conda:
        "../envs/gffutils.yaml"
    shell:
        """
        set -euo pipefail
        python workflow/scripts/gtf2bed.py \
            --input {input} \
            --bed {output.bed} \
            --db {output.db} \
            > {log} 2>&1
        """


rule rseqc_junction_annotation:
    input:
        bam=star_bam_path("{sample}", "{unit}"),
        bed=cohort_path("qc", "rseqc", "annotation.bed"),
    output:
        bed=rseqc_output_path("{sample}", "{unit}", "junctionanno.junction.bed"),
    priority: 1
    log:
        "logs/rseqc/rseqc_junction_annotation/{sample}_{unit}.log",
    benchmark:
        "logs/rseqc/rseqc_junction_annotation/{sample}_{unit}.bench.tsv"
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix=rseqc_output_path("{sample}", "{unit}", "junctionanno.junction"),
    threads: 32
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} > {log} 2>&1;
        """


rule rseqc_junction_saturation:
    input:
        bam=star_bam_path("{sample}", "{unit}"),
        bed=cohort_path("qc", "rseqc", "annotation.bed"),
    output:
        pdf=rseqc_output_path(
            "{sample}", "{unit}", "junctionsat.junctionSaturation_plot.pdf"
        ),
    priority: 1
    threads: 32
    log:
        "logs/rseqc/rseqc_junction_saturation/{sample}_{unit}.log",
    benchmark:
        "logs/benchmarks/rseqc_junction_saturation/{sample}_{unit}.bench.tsv",
    params:
        extra=r"-q 255",
        prefix=rseqc_output_path(
            "{sample}", "{unit}", "junctionsat.junctionSaturation_plot"
        ),
    conda:
        "../envs/rseqc.yaml"
    shell:
        """junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} > {log} 2>&1;
	"""


rule rseqc_stat:
    input:
        star_bam_path("{sample}", "{unit}"),
    output:
        stats=rseqc_output_path("{sample}", "{unit}", "stats.txt"),
    priority: 1
    log:
        "logs/rseqc/rseqc_stat/{sample}_{unit}.log",
    benchmark:
        "logs/benchmarks/rseqc_stat/{sample}_{unit}.bench.tsv",
    threads: 32
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        bam_stat.py -i {input} > {output.stats} 2> {log};
        """


rule rseqc_infer:
    input:
        bam=star_bam_path("{sample}", "{unit}"),
        bed=cohort_path("qc", "rseqc", "annotation.bed"),
    output:
        txt=rseqc_output_path("{sample}", "{unit}", "infer_experiment.txt"),
    priority: 1
    log:
        "logs/rseqc/rseqc_infer/{sample}_{unit}.log",
    benchmark:
        "logs/benchmarks/rseqc_infer/{sample}_{unit}.bench.tsv",
    threads: 32
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        infer_experiment.py -r {input.bed} -i {input.bam} > {output.txt} 2> {log};
        """


rule rseqc_innerdis:
    input:
        bam=star_bam_path("{sample}", "{unit}"),
        bed=cohort_path("qc", "rseqc", "annotation.bed"),
    output:
        txt=rseqc_output_path(
            "{sample}", "{unit}", "inner_distance_freq.inner_distance.txt"
        ),
    priority: 1
    log:
        "logs/rseqc/rseqc_innerdis/{sample}_{unit}.log",
    benchmark:
        "logs/benchmarks/rseqc_innerdis/{sample}_{unit}.bench.tsv",
    params:
        prefix=rseqc_output_path(
            "{sample}", "{unit}", "inner_distance_freq.inner_distance"
        ),
    threads: 32
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1;
        """


rule rseqc_readdis:
    input:
        bam=star_bam_path("{sample}", "{unit}"),
        bed=cohort_path("qc", "rseqc", "annotation.bed"),
    threads: 32
    output:
        txt=rseqc_output_path("{sample}", "{unit}", "readdistribution.txt"),
    priority: 1
    log:
        "logs/rseqc/rseqc_readdis/{sample}_{unit}.log",
    benchmark:
        "logs/benchmarks/rseqc_readdis/{sample}_{unit}.bench.tsv",
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        read_distribution.py -r {input.bed} -i {input.bam} > {output.txt} 2> {log};
        """


rule rseqc_readdup:
    input:
        star_bam_path("{sample}", "{unit}"),
    output:
        pdf=rseqc_output_path("{sample}", "{unit}", "readdup.DupRate_plot.pdf"),
    threads: 32
    priority: 1
    log:
        "logs/rseqc/rseqc_readdup/{sample}_{unit}.log",
    benchmark:
        "logs/benchmarks/rseqc_readdup/{sample}_{unit}.bench.tsv",
    params:
        prefix=rseqc_output_path("{sample}", "{unit}", "readdup.DupRate_plot"),
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1;
        """

rule rseqc_readgc:
    input:
        star_bam_path("{sample}", "{unit}"),
    output:
        pdf=rseqc_output_path("{sample}", "{unit}", "readgc.GC_plot.pdf"),
    threads: 32
    priority: 1
    log:
        "logs/rseqc/rseqc_readgc/{sample}_{unit}.log",
    benchmark:
        "logs/benchmarks/rseqc_readgc/{sample}_{unit}.bench.tsv",
    params:
        prefix=rseqc_output_path("{sample}", "{unit}", "readgc.GC_plot"),
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        read_GC.py -i {input} -o {params.prefix} > {log} 2>&1;
        """

rule multiqc:
    input:
        lambda wildcards: get_multiqc_inputs(wildcards),
    threads: 32
    output:
        html=cohort_path("reports", "multiqc", "multiqc_report.html"),
    log:
        "logs/multiqc.log",
    benchmark:
        "logs/benchmarks/multiqc.bench.tsv",
    container:
        "docker://daylilyinformatics/daylily_multiqc:0.2"
    shell:
        """
        outdir={output.html.rpartition('/')[0]}; \
        mkdir -p "$outdir"; \
        multiqc -f \
         --template default \
        --filename {output.html.rpartition('/')[-1]} \
        --outdir "$outdir" \
        -i 'deseq2 RNASEQ QC' {ANALYSIS_ROOT} >> {log} 2>&1;
        """
