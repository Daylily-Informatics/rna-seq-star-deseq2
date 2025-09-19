from snakemake.exceptions import WorkflowError


def _get_rsem_aligner():
    aligner = config.get("rsem_aligner", "bowtie2")
    if isinstance(aligner, str):
        aligner = aligner.lower()
    else:
        raise WorkflowError(
            f"Expected 'rsem_aligner' to be a string but received {type(aligner)}."
        )
    if aligner not in {"star", "bowtie2"}:
        raise WorkflowError(
            "Unsupported RSEM aligner '{aligner}'. Choose 'star' or 'bowtie2'.".format(
                aligner=aligner
            )
        )
    return aligner


RSEM_ALIGNER = _get_rsem_aligner()


rule rsem_index_star:
    input:
        fasta="resources/genome.fasta",
        gtf="resources/genome.gtf",
    output:
        transcripts="resources/star_rsem.transcripts.fa",
    params:
        prefix=lambda wc, output: output.transcripts.replace(".transcripts.fa", ""),
        extra=config["params"]["rsem"],
    log:
        "logs/rsem/prepare_reference.log",
    threads: 190
    cache: True
    container: "docker://daylilyinformatics/rsem:1.3.3.5"
    shell:
        """
        rsem-prepare-reference -p {threads} --star --gtf {input.gtf} {input.fasta} {params.extra} {params.prefix} &> {log}
        """

rule rsem_index_bowtie:
    input:
        fasta="resources/genome.fasta",
        gtf="resources/genome.gtf",
    output:
        transcripts="resources/bowtie_rsem.transcripts.fa",
    params:
        prefix=lambda wc, output: output.transcripts.replace(".transcripts.fa", ""),
        extra=config["params"]["rsem"],
    log:
        "logs/rsem/prepare_reference.log",
    threads: 190
    cache: True
    container: "docker://daylilyinformatics/rsem:1.3.3.5"
    shell:
        """
        rsem-prepare-reference -p {threads} --bowtie --gtf {input.gtf} {input.fasta} {params.extra} {params.prefix} &> {log}
        """

rule rsem_index_bowtie2:
    input:
        fasta="resources/genome.fasta",
        gtf="resources/genome.gtf",
    output:
        transcripts="resources/bowtie2_rsem.transcripts.fa",
    params:
        prefix=lambda wc, output: output.transcripts.replace(".transcripts.fa", ""),
        extra=config["params"]["rsem"],
    log:
        "logs/rsem/prepare_reference.log",
    threads: 190
    cache: True
    container: "docker://daylilyinformatics/rsem:1.3.3.5"
    shell:
        """
        rsem-prepare-reference -p {threads} --bowtie2 --gtf {input.gtf} {input.fasta} {params.extra} {params.prefix} &> {log}
        """



rule rsem_star:
    input:
        unpack(get_fq),
        ref="resources/star_rsem.transcripts.fa",
    output:
        genes="results/rsem/{sample}_{unit}.genes.results",
        isoforms="results/rsem/{sample}_{unit}.isoforms.results",
        bam="results/rsem/{sample}_{unit}.genome.sorted.bam",
        bai="results/rsem/{sample}_{unit}.genome.sorted.bam.bai",
        wig="results/rsem/{sample}_{unit}.genome.sorted.wig",
    log:
        "logs/rsem/{sample}_{unit}.star.log",
    benchmark:
        "logs/rsem/{sample}_{unit}.star.bench.tsv",
    threads: 190
    params:
        prefix=lambda wc, input: input.ref.replace(".transcripts.fa", ""),
        extra=config["params"]["rsem"],
        paired=lambda wc: "--paired-end" if is_paired_end(wc.sample) else "",
        fq_inputs=lambda wc, input: " ".join([input.fq1] + ([input.fq2] if is_paired_end(wc.sample) else [])),
        star_gzipped=lambda wc, input: (
            "--star-gzipped-read-file"
                if any(
                str(f).endswith(".gz")
                for f in [input.fq1]
                + ([input.fq2] if is_paired_end(wc.sample) else [])
            )
            else ""
        ),
    container: "docker://daylilyinformatics/rsem:1.3.3.5"
    shell:
        """
        (
        rsem-calculate-expression --star {params.star_gzipped} {params.paired} --sort-bam-by-coordinate --output-genome-bam -p {threads} {params.extra} {params.fq_inputs} {params.prefix} results/rsem/{wildcards.sample}_{wildcards.unit} &&
        rsem-bam2wig {output.bam} {output.wig} {wildcards.sample}_{wildcards.unit}
        ) &> {log}
        """



#    rule rsem_bowtie2_quant:
#        input:
#            bam="results/rsem/{sample}_{unit}.genome.sorted.bam",
#            bai="results/rsem/{sample}_{unit}.genome.sorted.bam.bai",
#            ref="resources/bowtie2_rsem.transcripts.fa",
#        output:
#            genes="results/rsem/x{sample}_{unit}.genes.results",
#            isoforms="results/rsem/x{sample}_{unit}.isoforms.results",
#        log:
#            "logs/rsem/{sample}_{unit}.log",
#        benchmark:
#            "logs/rsem/{sample}_{unit}.bench.tsv",
#        threads: 190
#        params:
#            prefix=lambda wc, input: input.ref.replace(".transcripts.fa", ""),
#            extra=config["params"]["rsem"],
#            paired=lambda wc: "--paired-end" if is_paired_end(wc.sample) else "",
#        container: "docker://daylilyinformatics/rsem:1.3.3.5"
#        shell:
#            """
#            rsem-calculate-expression --bowtie2 --alignments {params.paired} -p {threads} {params.extra} {input.bam} {params.prefix} results/rsem/{wildcards.sample}_{wildcards.unit} &> {log}
#            """
#
#    rule rsem_bowtie2:
#        input:
#            unpack(get_fq),
#            ref="resources/bowtie2_rsem.transcripts.fa",
#        output:
#            genes="results/rsem/{sample}_{unit}.genes.results",
#            isoforms="results/rsem/{sample}_{unit}.isoforms.results",
#            bam="results/rsem/{sample}_{unit}.genome.sorted.bam",
#            bai="results/rsem/{sample}_{unit}.genome.sorted.bam.bai",
#            wig="results/rsem/{sample}_{unit}.genome.sorted.wig",
#        log:
#            "logs/rsem/{sample}_{unit}.bowtie.log",
#        benchmark:
#            "logs/rsem/{sample}_{unit}.bowtie.bench.tsv",
#        threads: 190
#        params:
#            prefix=lambda wc, input: input.ref.replace(".transcripts.fa", ""),
#            extra=config["params"]["rsem"],
#            paired=lambda wc: "--paired-end" if is_paired_end(wc.sample) else "",
#            fq_inputs=lambda wc, input: " ".join([input.fq1] + ([input.fq2] if is_paired_end(wc.sample) else [])),
#        container: "docker://daylilyinformatics/rsem:1.3.3.5"
#        shell:
#            """
#            (
#            rsem-calculate-expression --bowtie2 {params.paired} --sort-bam-by-coordinate --output-genome-bam -p {threads} {params.extra} {params.fq_inputs} {params.prefix} results/rsem/{wildcards.sample}_{wildcards.unit} &&
#            rsem-bam2wig {output.bam} {output.wig} {wildcards.sample}_{wildcards.unit}
#            ) &> {log}
#            """

#   ruleorder:
#        rsem_bowtie2 > rsem_bowtie2_quant

