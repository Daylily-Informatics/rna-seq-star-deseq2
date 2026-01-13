rule align:
    input:
        fq1=lambda wildcards: get_fq(wildcards)["fq1"],
        fq2=lambda wildcards: get_fq(wildcards).get("fq2", []),
        index="resources/star_genome",
        gtf="resources/genome.gtf",
    output:
        aln=star_bam_path("{sample}", "{unit}"),
        reads_per_gene=star_counts_path("{sample}", "{unit}"),
    log:
        "logs/star/{sample}_{unit}.log",
    benchmark:
        "logs/star/{sample}_{unit}.bench.tsv",
    threads: 127
    resources:
        tmpdir="/dev/shm"
    params:
        idx=lambda wc, input: input.index,
        extra=lambda wc, input: f'--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile {input.gtf} {config["params"]["star"]} ',
    wrapper:
        "v3.5.3/bio/star/align"
