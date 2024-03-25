# Getting Started

## Tutorial

A tutorial is available with a small(ish) dataset where biologically meaningful
results can be produced. This can help get an understanding of a good workflow
to use different modules. You can also follow along with your own data and just
skip analyses you don't want. If you prefer to just jump in instead, below
describes how to quickly get a new project up and running.

## Requirements

This pipeline can be run on Linux systems with Conda and Apptainer/Singularity
installed. All other dependencies will be handled with the workflow, and thus,
sufficient storage space is needed for these installations (~10GB, but this
needs verification). It can be run on a local workstation with sufficient
resources and storage space (dataset dependent), but is aimed at execution on
high performance computing systems with job queuing systems.

Data-wise, you'll need a reference genome (uncompressed) and some sequencing
data for your samples. The latter can be either raw fastq files, bam alignments
to the reference, or accession numbers for already published fastq files.

## Deploying the workflow

The pipeline can be deployed in two ways: (1) using Snakedeploy which will
deploy the pipeline as a module (recommended); (2) clone the repository at the
version/branch you prefer (recommended if you will change any workflow code).

Both methods require a Snakemake environment to run the pipeline in.

### Preparing the environment

First, create an environment for Snakemake, including Snakedeploy if you intend
to deploy that way:

```bash
mamba create -c conda-forge -c bioconda --name snakemake snakemake snakedeploy
```

If you already have a Snakemake environment, you can use that, so long as
it has `snakemake` (not just `snakemake-minimal`) installed. Snakemake
versions >=7.25 are likely to work, but most testing is on 7.32.4. You can
change the name of the environment if you'd like to isolate it any other
Snakemake environments you might have. The only main dependency of this
pipeline is `snakemake`.

Activate the Snakemake environment:

```bash
conda activate snakemake
```

### Deploying with Snakedeploy

Make your working directory:

```bash
mkdir -p /path/to/work-dir
cd /path/to/work-dir
```

And deploy the workflow, using the tag for the version you want to deploy:

```bash
snakedeploy deploy-workflow https://github.com/zjnolen/angsd-snakemake-pipeline . --tag v0.2.0
```

This will generate a simple Snakefile in a `workflow` folder that loads the
pipeline as a module. It will also download the template `config.yaml`,
`samples.tsv`, and `units.tsv` in the `config` folder.

### Cloning from GitHub

Go to the folder you would like you working directory to be created in and
clone the GitHub repo:

```bash
git clone https://github.com/zjnolen/angsd-snakemake-pipeline.git
```

If you would like, you can change the name of the directory:

```bash
mv angsd-snakemake-pipeline work-dir-name
```

Move into the working directory (`angsd-snakemake-pipeline` or `work-dir-name`
if you changed it) and checkout the version you would like to use:

```bash
git checkout v0.2.0
```

This can also be used to checkout specific branches or commits.

## Configuring the workflow

Now you are ready to configure the workflow, see the documentation for that
[here](config.md).
