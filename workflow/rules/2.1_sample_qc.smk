# Rules for various assessments of samples, mainly regarding alignment quality


localrules:
    combine_sample_qc,


rule samtools_flagstat:
    """
    Estimate statistics for samtools flags. Needed to assess mapping percentage
    """
    input:
        "results/mapping/mapped/{prefix}.bam",
    output:
        "results/mapping/mapped/{prefix}.flagstat",
    log:
        "logs/mapping/samtools/flagstat/{prefix}.log",
    benchmark:
        "benchmarks/mapping/samtools/flagstat/{prefix}.log"
    wrapper:
        "v2.6.0/bio/samtools/flagstat"


rule qualimap:
    """
    Estimate general mapping statistics for each final sample bam file.
    """
    input:
        unpack(get_final_bam),
    output:
        directory("results/mapping/qc/qualimap/{sample}.{ref}"),
        pdf=report(
            "results/mapping/qc/qualimap/{sample}.{ref}/report.pdf",
            category="Quality Control",
            subcategory="Mapping Reports",
            labels={"Sample": "{sample}", "Ref": "{ref}", "Type": "Qualimap Report"},
        ),
        txt="results/mapping/qc/qualimap/{sample}.{ref}/genome_results.txt",
    params:
        extra="-outformat pdf",
    log:
        "logs/mapping/qualimap/{sample}.{ref}.log",
    benchmark:
        "benchmarks/mapping/qualimap/{sample}.{ref}.log"
    resources:
        runtime=360,
    wrapper:
        "v2.6.0/bio/qualimap/bamqc"


rule endo_cont:
    """
    Estimate the proportion of reads mapping to the reference as a proxy for endogenous
    content.
    """
    input:
        get_endo_cont_stat,
    output:
        "results/mapping/qc/endogenous_content/{sample}.{ref}.endo",
    conda:
        "../envs/shell.yaml"
    log:
        "logs/mapping/endogenous_content/{sample}.{ref}.log",
    benchmark:
        "benchmarks/mapping/endogenous_content/{sample}.{ref}.log"
    shell:
        r"""
        (total=$(grep -E "^[0-9]+ \+ [0-9]+ in total" {input} \
            | awk '{{print $1}}')
        mapped=$(grep -E "^[0-9]+ \+ [0-9]+ mapped" {input} \
            | awk '{{print $1}}')
        primary=$(grep -E "^[0-9]+ \+ [0-9]+ primary mapped" {input} \
            | awk '{{print $1}}')

        echo $total $mapped $primary {wildcards.sample} | \
            awk '{{printf "%s\t%.3f\t%.3f\n",$4,$2/$1*100,$3/$1*100}}' \
            > {output}) 2> {log}
        """


rule compile_endo_cont:
    """
    Merge per sample endogenous content estimates into a single table.
    """
    input:
        lambda w: expand(
            "results/mapping/qc/endogenous_content/{sample}.{{ref}}.endo",
            sample=get_samples_from_pop("all"),
        ),
    output:
        "results/datasets/{dataset}/qc/endogenous_content/{dataset}.{ref}_all.endo.tsv",
    log:
        "logs/datasets/{dataset}/qc/endogenous_content/{dataset}.{ref}_compile-endocont.log",
    benchmark:
        "benchmarks/datasets/{dataset}/qc/endogenous_content/{dataset}.{ref}_compile-endocont.log"
    conda:
        "../envs/shell.yaml"
    resources:
        runtime=lambda wildcards, attempt: attempt * 15,
    shell:
        """
        (printf "sample\tperc.map\tperc.prim.map\n" > {output}
        cat {input} >> {output}) 2> {log}
        """


rule ind_unfiltered_depth:
    """
    Estimate unfiltered sample depth, only removing reads using default ANGSD filters
    """
    input:
        bamlist="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
        bams=get_bamlist_bams,
        bais=get_bamlist_bais,
        ref="results/ref/{ref}/{ref}.fa",
    output:
        sample_hist="results/mapping/qc/ind_depth/unfiltered/{dataset}.{ref}_{population}{dp}_allsites-unfilt.depthSample",
        global_hist=temp(
            "results/mapping/qc/ind_depth/unfiltered/{dataset}.{ref}_{population}{dp}_allsites-unfilt.depthGlobal"
        ),
        arg="results/mapping/qc/ind_depth/unfiltered/{dataset}.{ref}_{population}{dp}_allsites-unfilt.arg",
    log:
        "logs/mapping/ind_depth/unfiltered/{dataset}.{ref}_{population}{dp}_allsites-unfilt.log",
    benchmark:
        "benchmarks/mapping/ind_depth/unfiltered/{dataset}.{ref}_{population}{dp}_allsites-unfilt.log"
    container:
        angsd_container
    params:
        out=lambda w, output: os.path.splitext(output.arg)[0],
        maxdepth=config["params"]["angsd"]["maxdepth"],
    threads: lambda wildcards, attempt: attempt
    resources:
        runtime=lambda wildcards, attempt: attempt * 120,
    shell:
        """
        angsd -doDepth 1 -doCounts 1 -maxDepth {params.maxdepth} \
            -bam {input.bamlist} -nThreads {threads} -out {params.out} &> {log}
        """


rule ind_filtered_depth:
    """
    Estimate depth at positions using filters for main workflow. This describes the
    depth of positions that will go into likelihood calculations
    """
    input:
        bamlist="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
        bams=get_bamlist_bams,
        bais=get_bamlist_bais,
        ref="results/ref/{ref}/{ref}.fa",
        sites="results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_{sites}-filts.sites",
        idx="results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_{sites}-filts.sites.idx",
    output:
        sample_hist="results/datasets/{dataset}/qc/ind_depth/filtered/{dataset}.{ref}_{population}{dp}_{sites}-filts.depthSample",
        global_hist=temp(
            "results/datasets/{dataset}/qc/ind_depth/filtered/{dataset}.{ref}_{population}{dp}_{sites}-filts.depthGlobal"
        ),
        arg="results/datasets/{dataset}/qc/ind_depth/filtered/{dataset}.{ref}_{population}{dp}_{sites}-filts.arg",
    log:
        "logs/{dataset}/ind_depth/filtered/{dataset}.{ref}_{population}{dp}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/ind_depth/filtered/{dataset}.{ref}_{population}{dp}_{sites}-filts.log"
    container:
        angsd_container
    params:
        extra=config["params"]["angsd"]["extra"],
        mapQ=config["mapQ"],
        baseQ=config["baseQ"],
        out=lambda w, output: os.path.splitext(output.arg)[0],
        maxdepth=config["params"]["angsd"]["maxdepth"],
    threads: lambda wildcards, attempt: attempt
    resources:
        runtime=lambda wildcards, attempt: attempt * 60,
    shell:
        """
        angsd -doDepth 1 -doCounts 1 -maxDepth {params.maxdepth} \
            -bam {input.bamlist} -ref {input.ref} -nThreads {threads} \
            {params.extra} -minMapQ {params.mapQ} -minQ {params.baseQ} \
            -sites {input.sites} -out {params.out} &> {log}
        """


rule summarize_ind_depth:
    """
    Get mean and standard deviation of depth from individual depth distributions
    """
    input:
        sample_hist="{prefix}{dataset}.{ref}_{sample}{dp}_{group}.depthSample",
        bed=get_total_bed,
    output:
        sample_summ="{prefix}{dataset}.{ref}_{sample}{dp}_{group}.depth.sum",
    log:
        "logs/summarize_ind_depth/{prefix}{dataset}.{ref}_{sample}{dp}_{group}.log",
    benchmark:
        "benchmarks/summarize_ind_depth/{prefix}{dataset}.{ref}_{sample}{dp}_{group}.log"
    conda:
        "../envs/r.yaml"
    threads: lambda wildcards, attempt: attempt
    script:
        "../scripts/calc_depth.R"


rule merge_ind_depth:
    """
    Combine depth summaries for all individuals
    """
    input:
        depth=lambda w: expand(
            "{{prefix}}{{dataset}}.{{ref}}_{sample}{{dp}}_{{group}}.depthSample",
            sample=get_samples_from_pop("all"),
        ),
        summary=lambda w: expand(
            "{{prefix}}{{dataset}}.{{ref}}_{sample}{{dp}}_{{group}}.depth.sum",
            sample=get_samples_from_pop("all"),
        ),
    output:
        dep="{prefix}{dataset}.{ref}_all{dp}_{group}.depth",
        sum="{prefix}{dataset}.{ref}_all{dp}_{group}.depth.sum",
    log:
        "logs/merge_depth/{prefix}{dataset}.{ref}_all{dp}_{group}.log",
    benchmark:
        "benchmarks/merge_depth/{prefix}{dataset}.{ref}_all{dp}_{group}.log"
    conda:
        "../envs/shell.yaml"
    shell:
        """
        (cat {input.depth} > {output.dep}
        printf "sample\t{wildcards.group}.depth.mean\t{wildcards.group}.depth.stdev\n" \
            > {output.sum}
        cat {input.summary} >> {output.sum}) 2> {log}
        """


rule combine_sample_qc:
    """
    Compile summary table of all sample QC measures
    """
    input:
        unpack(get_sample_qcs),
    output:
        "results/datasets/{dataset}/qc/{dataset}.{ref}_all{dp}.sampleqc.tsv",
    log:
        "logs/datasets/{dataset}/combine_sample_qc/{dataset}.{ref}{dp}.log",
    benchmark:
        "benchmarks/datasets/{dataset}/combine_sample_qc/{dataset}.{ref}{dp}.log"
    conda:
        "../envs/shell.yaml"
    shadow:
        "minimal"
    shell:
        """
        (for i in {input}; do
            head -n 1 $i > ${{i}}.tmp
            tail -n +2 $i | sort -k1 >> ${{i}}.tmp
        done

        cut -d '\t' -f 1 {input.inds}.tmp > {output}

        for i in {input}; do
            cut -d '\t' -f 2- ${{i}}.tmp | paste {output} - > {output}.tmp
            mv {output}.tmp {output}
        done) 2> {log}
        """


rule sample_qc_summary:
    """
    Convert sample QC summary table to html for report
    """
    input:
        "results/datasets/{dataset}/qc/{dataset}.{ref}_all{dp}.sampleqc.tsv",
    output:
        report(
            "results/datasets/{dataset}/qc/{dataset}.{ref}_all{dp}.sampleqc.html",
            category="Quality Control",
            subcategory="Sample coverage and endogenous content",
            labels={"Type": "Table"},
        ),
    log:
        "logs/{dataset}/combine_sample_qc/{dataset}.{ref}{dp}_tsv2html.log",
    benchmark:
        "benchmarks/{dataset}/combine_sample_qc/{dataset}.{ref}{dp}_tsv2html.log"
    conda:
        "../envs/r-rectable.yaml"
    script:
        "../scripts/tsv2html.Rmd"
