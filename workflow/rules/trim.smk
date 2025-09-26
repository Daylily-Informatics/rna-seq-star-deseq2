rule get_sra:
    output:
        "sra/{accession}_1.fastq",
        "sra/{accession}_2.fastq",
    threads: 31
    log:
        "logs/get-sra/{accession}.log",
    benchmark:
        "logs/benchmarks/get_sra/{accession}.bench.tsv",
    wrapper:
        "v3.5.3/bio/sra-tools/fasterq-dump"


rule cutadapt_pipe:
    input:
        get_cutadapt_pipe_input,
    output:
        pipe("pipe/cutadapt/{sample}/{unit}.{fq}.{ext}"),
    log:
        "logs/pipe-fastqs/catadapt/{sample}_{unit}.{fq}.{ext}.log",
    benchmark:
        "logs/pipe-fastqs/catadapt/{sample}_{unit}.{fq}.{ext}.bench.tsv",
    wildcard_constraints:
        ext=r"fastq|fastq\.gz",
    threads: 1  ## this does something special when running using pipe() output directives
    shell:
        "cat {input} > {output} 2> {log}"


rule cutadapt_pe:
    input:
        get_cutadapt_input,
    output:
        fastq1=trimmed_fastq_path("{sample}", "{unit}", "fq1"),
        fastq2=trimmed_fastq_path("{sample}", "{unit}", "fq2"),
        qc=trimmed_qc_path("{sample}", "{unit}", "paired.qc.txt"),
    log:
        "logs/cutadapt/{sample}_{unit}.log",
    benchmark:
        "logs/cutadapt/{sample}_{unit}.pe.bench.tsv",
    params:
        extra=config["params"]["cutadapt-pe"],
        adapters=lambda w: str(units.loc[w.sample].loc[w.unit, "adapters"]),
    threads: 64
    wrapper:
        "v3.5.3/bio/cutadapt/pe"


rule cutadapt_se:
    input:
        get_cutadapt_input,
    output:
        fastq=trimmed_fastq_path("{sample}", "{unit}", "single"),
        qc=trimmed_qc_path("{sample}", "{unit}", "single.qc.txt"),
    log:
        "logs/cutadapt/{sample}_{unit}.log",
    benchmark:
        "logs/cutadapt/{sample}_{unit}.se.bench.tsv",
    params:
        extra=config["params"]["cutadapt-se"],
        adapters=lambda w: str(units.loc[w.sample].loc[w.unit, "adapters"]),
    threads: 64
    wrapper:
        "v3.5.3/bio/cutadapt/se"
