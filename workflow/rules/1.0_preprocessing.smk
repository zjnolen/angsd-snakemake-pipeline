# Rules for adapter removal, quality trimming, and pair collapsing (for historical
# samples only) with fastp


rule fastp_mergedout:
    """Process historical reads with fastp, collapsing overlapping read pairs"""
    input:
        unpack(get_raw_fastq),
    output:
        trimmed=temp(
            expand(
                "results/preprocessing/fastp/{{sample}}_{{unit}}_{{lib}}.{read}.discard.fastq.gz",
                read=["R1", "R2"],
            )
        ),
        merged=temp("results/preprocessing/fastp/{sample}_{unit}_{lib}.merged.fastq.gz"),
        html=report(
            "results/preprocessing/qc/fastp/{sample}_{unit}_{lib}_paired.html",
            category="Quality Control",
            subcategory="Trimming Reports",
            labels={
                "Sample": "{sample}",
                "Unit": "{unit}",
                "Lib": "{lib}",
                "Type": "fastp Report",
            },
        ),
        json="results/preprocessing/qc/fastp/{sample}_{unit}_{lib}_paired.json",
    log:
        "logs/preprocessing/fastp/{sample}_{unit}_{lib}.merged.log",
    benchmark:
        "benchmarks/preprocessing/fastp/{sample}_{unit}_{lib}.merged.log"
    params:
        extra=config["params"]["fastp"]["extra"] + " --merge",
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        runtime=lambda wildcards, attempt: attempt * 240,
    wrapper:
        "v2.5.0/bio/fastp"


rule fastp_pairedout:
    """Process modern reads with fastp, trimming adapters and low quality bases"""
    input:
        unpack(get_raw_fastq),
    output:
        trimmed=temp(
            expand(
                "results/preprocessing/fastp/{{sample}}_{{unit}}_{{lib}}.{read}.fastq.gz",
                read=["R1", "R2"],
            )
        ),
        html=report(
            "results/preprocessing/qc/fastp/{sample}_{unit}_{lib}_paired.html",
            category="Quality Control",
            subcategory="Trimming Reports",
            labels={
                "Sample": "{sample}",
                "Unit": "{unit}",
                "Lib": "{lib}",
                "Type": "fastp Report",
            },
        ),
        json="results/preprocessing/qc/fastp/{sample}_{unit}_{lib}_paired.json",
    log:
        "logs/preprocessing/fastp/{sample}_{unit}_{lib}.paired.log",
    benchmark:
        "benchmarks/preprocessing/fastp/{sample}_{unit}_{lib}.paired.log"
    params:
        extra=config["params"]["fastp"]["extra"],
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        runtime=lambda wildcards, attempt: attempt * 240,
    wrapper:
        "v2.5.0/bio/fastp"
