def get_align_inputs(wildcards):
    files = dict(get_fq(wildcards))
    files.update({
        "index": "resources/star_genome",
        "gtf": "resources/genome.gtf",
    })
    return files


rule align:
    input:
        lambda wildcards: get_align_inputs(wildcards),
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
