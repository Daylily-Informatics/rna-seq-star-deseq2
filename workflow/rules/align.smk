rule align:
    input:
        unpack(get_fq),
        index="resources/star_genome",
        gtf="resources/genome.gtf",
    output:
        aln=lambda wc: star_bam_path(wc.sample, wc.unit),
        reads_per_gene=lambda wc: star_counts_path(wc.sample, wc.unit),
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
