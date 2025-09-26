# Snakemake workflow: rna-seq-star-deseq2 on AWS ParallelCluster

This repository packages a production-oriented RNA-seq differential expression workflow that is tuned for Daylily Informatics' AWS ParallelCluster environments. It layers AWS-aware execution, caching, and reporting conventions on top of the familiar STAR ➜ feature counting ➜ DESeq2 analysis stack so teams can iterate on transcriptomics projects without reinventing infrastructure.

## What this repository provides

- **Turnkey Snakemake workflow** for paired-end RNA-seq, orchestrated by the `pcluster-slurm` Snakemake executor plugin to seamlessly submit jobs to AWS ParallelCluster Slurm schedulers.
- **Reference and resource management helpers** that cache STAR indices, Singularity images, and Conda environments on shared FSx/NFS mounts for rapid re-use across projects.
- **Quality control dashboards** generated with FastQC, MultiQC, and RSeQC, with example outputs captured under [`docs/results_tree.log`](docs/results_tree.log) for quick browsing of the report structure.
- **Differential expression reporting** via DESeq2, producing count matrices and annotated contrasts suitable for downstream visualization or knowledge bases.
- **Reproducible environment bootstrapping** with helper scripts (for example, `bin/install_miniconda`) and configuration templates that reduce the friction of onboarding new analysts to Daylily's RNA-seq stack.

### Example artifacts

The `docs/` directory contains reference directory trees from a representative run:

- [`docs/resources_tree.log`](docs/resources_tree.log) shows the curated genome resources bundle (FASTA, GTF, and STAR genome directory) expected by the workflow.
- [`docs/results_tree.log`](docs/results_tree.log) enumerates the outputs produced during a small treated-vs-untreated comparison, including STAR alignment BAMs, read count tables, DESeq2 normalized counts, MA plots, and QC summaries.

These examples illustrate the layout teams can rely on when integrating results with downstream analytics or long-term storage policies.

## Common use cases

- **Rapid validation of new reference builds** by rebuilding STAR indices on the cluster and verifying the generated QC reports against expected metrics.
- **Budget-aware large cohort processing** where the workflow's support for Slurm comments (e.g., `SMK_SLURM_COMMENT`) and job partition targeting simplifies cost tracking across Daylily's organizational projects.
- **Iterative methods development** for Daylily scientists who need to test alternative quantification or normalization strategies while preserving a stable baseline pipeline for comparison.

## Prerequisites

### `daylily-ephemeral-cluster` (using AWS Parallel Cluster)
- This has been developed to run on an AWS Parallel Cluster slurm headnode, specifically one created using [https://github.com/Daylily-Informatics/daylily-ephemeral-cluster](https://github.com/Daylily-Informatics/daylily-ephemeral-cluster).

### Other
1. An active AWS ParallelCluster deployment with Slurm (either a self-managed cluster or the [daylily-ephemeral-cluster](https://github.com/Daylily-Informatics/daylily-ephemeral-cluster)).
2. Conda or Miniconda available on the head node. The provided `bin/install_miniconda` script can be used if Conda is not already present.
3. User access to shared FSx (or comparable) storage for caching environments, containers, and references.

## Getting started

Clone the repository (it includes small example data sets):

```bash
git clone git@github.com:Daylily-Informatics/rna-seq-star-deseq2.git
cd rna-seq-star-deseq2
```

### Build the Snakemake environment (Snakemake v9.11.4.1)

Install the Daylily Informatics fork of Snakemake that bundles AWS ParallelCluster integration alongside the executor plugin dependencies.

```bash
conda create -n snakemake -c conda-forge tabulate yaml
conda activate snakemake
pip install git+https://github.com/Daylily-Informatics/snakemake-aws@v9.11.4.3
pip install snakemake-executor-plugin-pcluster-slurm==0.0.31
pip install snakedeploy
conda activate srrda
snakemake --version  # 9.11.4.1
```

### Prime caches (recommended for large-scale work)

```bash
conda activate srrda
mkdir -p /fsx/resources/environments/containers/ubuntu/rnaseq_cache/
export SNAKEMAKE_OUTPUT_CACHE=/fsx/resources/environments/containers/ubuntu/rnaseq_cache/
export TMPDIR=/fsx/scratch/
```

Prepare run descriptors:

```bash
cp config/units.tsv.template config/units.tsv
[[ "$(uname)" == "Darwin" ]] && sed -i "" "s|REGSUB_PWD|$PWD|g" config/units.tsv || sed -i "s|REGSUB_PWD|$PWD|g" config/units.tsv
```

Build Conda environments ahead of time to minimize surprises during production submissions (can take ~1 hour):

```bash
snakemake --use-conda --use-singularity \
  --singularity-prefix /fsx/resources/environments/containers/ubuntu/ \
  --singularity-args "-B /tmp:/tmp -B /fsx:/fsx -B /home/$USER:/home/$USER -B $PWD/:$PWD" \
  --conda-prefix /fsx/resources/environments/containers/ubuntu/ \
  --executor pcluster-slurm \
  --default-resources slurm_partition=i192,i128 runtime=86400 mem_mb=36900 tmpdir=/fsx/scratch \
  --cache -p \
  -k \
  --max-threads 20000 \
  --restart-times 2 \
  --cores 20000 -j 14 -n \
  --conda-create-envs-only
```

> Note: Run once with `--conda-create-envs-only` to populate environments, then rerun without `-n` to execute the workflow. Setting `--max-threads` and `--cores` above your head-node CPU count works around Slurm thread detection quirks on AWS ParallelCluster.

### Submit the pipeline

Update `config/units.tsv`, `config/samples.tsv`, and `config/config.yaml` with your project-specific metadata. Then launch the workflow:

```bash
snakemake --use-conda --use-singularity \
  --singularity-prefix /fsx/resources/environments/containers/ubuntu/ \
  --singularity-args "-B /tmp:/tmp -B /fsx:/fsx -B /home/$USER:/home/$USER -B $PWD/:$PWD" \
  --conda-prefix /fsx/resources/environments/containers/ubuntu/ \
  --executor pcluster-slurm \
  --default-resources slurm_partition=i192,i128 runtime=86400 mem_mb=36900 tmpdir=/fsx/scratch \
  --cache -p \
  -k \
  --restart-times 2 \
  --max-threads 20000 \
  --cores 20000 -j 14 \
  --include-aws-benchmark-metrics
```

Monitor job states with `watch squeue` and adjust partitions (`slurm_partition=...`) to match the compute queues defined for your cluster.

## Running with your data

1. Place sample FASTQ paths and associated metadata in `config/units.tsv` and `config/samples.tsv`.
2. Review `config/config.yaml` for alignment, quantification, and contrast options.
3. Perform a dry run with `-n` to validate DAG construction and resource requests before launching full scale analyses.

## Divergence notice

This repository has evolved beyond the original public RNA-seq workflow. Previous references to that project have been removed to reduce confusion; the documentation above reflects the Daylily-specific tooling now maintained here.
