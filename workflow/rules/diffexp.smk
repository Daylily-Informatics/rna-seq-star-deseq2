rule count_matrix:
    input:
        expand(
            "results/star/{unit.sample_name}_{unit.unit_name}/ReadsPerGene.out.tab",
            unit=units.itertuples(),
        ),
    output:
        "results/counts/all.tsv",
    log:
        "logs/count-matrix.log",
    params:
        samples=lambda wildcards: ",".join(units["sample_name"].tolist()),
        strands=lambda wildcards: ",".join(get_strandedness(units)),
    conda:
        "../envs/pandas.yaml"
    shell:
        """
        python workflow/scripts/count-matrix.py \
            --output {output} \
            --samples "{params.samples}" \
            --strands "{params.strands}" \
            {input} \
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
    conda:
        "../envs/biomart.yaml"
    shell:
        """
        line_count=$(awk 'END {{print NR}}' {input.counts})
        if [ "$line_count" -le 1 ]; then
            echo "Input {input.counts} has $line_count line(s); skipping gene symbol annotation." > {log}
            touch {output.symbol}
            touch {output.nodata}
        else
            rm -f {output.nodata}
            echo "Annotating gene symbols for {input.counts}." > {log}
            Rscript workflow/scripts/gene2symbol.R \
                --counts {input.counts} \
                --output {output.symbol} \
                --species {params.species} \
                >> {log} 2>&1
        fi
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
    threads: get_deseq2_threads()
    shell:
        "touch {output}"

