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
    threads: 8
    cache: True
    container: "docker://daylilyinformatics/rsem:1.3.3.4"
    shell:
        """
        rsem-prepare-reference --star --gtf {input.gtf} {input.fasta} {params.extra} {params.prefix} &> {log}
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
    threads: 8
    cache: True
    container: "docker://daylilyinformatics/rsem:1.3.3.4"
    shell:
        """
        rsem-prepare-reference --bowtie --gtf {input.gtf} {input.fasta} {params.extra} {params.prefix} &> {log}
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
    threads: 8
    cache: True
    container: "docker://daylilyinformatics/rsem:1.3.3.4"
    shell:
        """
        rsem-prepare-reference --bowtie2 --gtf {input.gtf} {input.fasta} {params.extra} {params.prefix} &> {log}
        """

rule rsem_bowtie2_quant:
    input:
        bam="results/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam",
        ref="resources/bowtie2_rsem.transcripts.fa",
    output:
        genes="results/rsem/{sample}_{unit}.genes.results",
        isoforms="results/rsem/{sample}_{unit}.isoforms.results",
    log:
        "logs/rsem/{sample}_{unit}.log",
    benchmark:
        "logs/rsem/{sample}_{unit}.bench.tsv",
    threads: 8
    params:
        prefix=lambda wc, input: input.ref.replace("/rsem.transcripts.fa", ""),
        extra=config["params"]["rsem"],
        paired=lambda wc: "--paired-end" if is_paired_end(wc.sample) else "",
    container: "docker://daylilyinformatics/rsem:1.3.3.4"
    shell:
        """
        rsem-calculate-expression --bowtie2 --alignments {params.paired} -p {threads} {params.extra} {input.bam} {params.prefix} results/rsem/{wildcards.sample}_{wildcards.unit} &> {log}
        """


rule rsem_bowtie2:
    input:
        unpack(get_fq),
        ref="resources/bowtie2_rsem.transcripts.fa",
    output:
        genes="results/rsem/{sample}_{unit}.genes.results",
        isoforms="results/rsem/{sample}_{unit}.isoforms.results",
        bam="results/rsem/{sample}_{unit}.genome.sorted.bam",
        bai="results/rsem/{sample}_{unit}.genome.sorted.bam.bai",
        wig="results/rsem/{sample}_{unit}.genome.sorted.wig",
    log:
        "logs/rsem/{sample}_{unit}.bowtie.log",
    benchmark:
        "logs/rsem/{sample}_{unit}.bowtie.bench.tsv",
    threads: 8
    params:
        prefix=lambda wc, input: input.ref.replace("/rsem.transcripts.fa", ""),
        extra=config["params"]["rsem"],
        paired=lambda wc: "--paired-end" if is_paired_end(wc.sample) else "",
        fq_inputs=lambda wc, input: " ".join([input.fq1] + ([input.fq2] if is_paired_end(wc.sample) else [])),
    container: "docker://daylilyinformatics/rsem:1.3.3.4"
    shell:
        """
        (
        rsem-calculate-expression --bowtie2 {params.paired} --sort-bam-by-coordinate --output-genome-bam -p {threads} {params.extra} {params.fq_inputs} {params.prefix} results/rsem/{wildcards.sample}_{wildcards.unit} &&
        rsem-bam2wig {output.bam} {output.wig} {wildcards.sample}_{wildcards.unit}
        ) &> {log}
        """

ruleorder:
    rsem_bowtie2 > rsem_bowtie2_quant

