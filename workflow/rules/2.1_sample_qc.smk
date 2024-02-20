# Rules for various assessments of samples, mainly regarding alignment quality


localrules:
    combine_sample_qc,


rule samtools_flagstat:
    """
    Estimate statistics for samtools flags. Needed to assess mapping percentage
    """
    input:
        "results/{prefix}.bam",
    output:
        "results/{prefix}.flagstat",
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


rule qualimap_userprovided:
    """
    Estimate general mapping statistics for each user-provided sample bam file.
    """
    input:
        unpack(get_final_bam),
    output:
        directory(
            "results/datasets/{dataset}/qc/user-provided-bams/qualimap/{sample}.{ref}"
        ),
        pdf=report(
            "results/datasets/{dataset}/qc/user-provided-bams/qualimap/{sample}.{ref}/report.pdf",
            category="Quality Control",
            subcategory="Mapping Reports",
            labels={"Sample": "{sample}", "Ref": "{ref}", "Type": "Qualimap Report"},
        ),
        txt="results/datasets/{dataset}/qc/user-provided-bams/qualimap/{sample}.{ref}/genome_results.txt",
    params:
        extra="-outformat pdf",
    log:
        "logs/mapping/qualimap/{dataset}.{sample}.{ref}.log",
    benchmark:
        "benchmarks/mapping/qualimap/{dataset}.{sample}.{ref}.log"
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
        unpack(get_endo_cont_stat),
    output:
        endo="results/datasets/{dataset}/qc/endogenous_content/{dataset}.{sample}.{ref}.endo",
    conda:
        "../envs/shell.yaml"
    log:
        "logs/datasets/{dataset}/qc/endogenous_content/{dataset}.{sample}.{ref}.log",
    benchmark:
        "benchmarks/datasets/{dataset}/qc/endogenous_content/{dataset}.{sample}.{ref}.log"
    script:
        "../scripts/calc_endocont.sh"


rule compile_endo_cont:
    """
    Merge per sample endogenous content estimates into a single table.
    """
    input:
        lambda w: expand(
            "results/datasets/{{dataset}}/qc/endogenous_content/{{dataset}}.{sample}.{{ref}}.endo",
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
        (printf "sample\tperc.collapsed.map\tperc.uncollapsed.map\tperc.total.map\n" > {output}
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
        unpack(filt_depth),
        bamlist="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
        bams=get_bamlist_bams,
        bais=get_bamlist_bais,
        ref="results/ref/{ref}/{ref}.fa",
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


rule doFasta:
    """
    Call consensus fasta of individual in filtered regions.
    """
    input:
        bam="results/datasets/{dataset}/bams/{sample}.{ref}.bam",
        bai="results/datasets/{dataset}/bams/{sample}.{ref}.bam.bai",
        #sites="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.sites",
        #idx="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.sites.idx",
    output:
        fa="results/datasets/{dataset}/fastas/{sample}.{ref}.consensus.fa.gz",
        arg="results/datasets/{dataset}/fastas/{sample}.{ref}.consensus.arg"
    log:
        "logs/{dataset}/angsd/doFasta/{sample}.{ref}.consensus.log",
    benchmark:
        "benchmarks/{dataset}/angsd/doFasta/{sample}.{ref}.consensus.log"
    container:
        angsd_container
    resources:
        runtime=lambda w, attempt: attempt * 360,
    params:
        extra=config["params"]["angsd"]["extra"],
        mapQ=config["mapQ"],
        baseQ=config["baseQ"],
        out=lambda w, output: os.path.splitext(output.arg)[0],
    shell:
        """
        angsd -doFasta 1 -i {input.bam} -nThreads {threads} {params.extra} \
            -minMapQ {params.mapQ} -minQ {params.baseQ} -doCounts 1 -out {params.out} \
            2> {log}
        """


rule samtools_faidx_sample_cons:
    """Index reference genome using samtools (fai index used by several tools)"""
    input:
        "results/datasets/{dataset}/fastas/{sample}.{ref}.consensus.fa.gz",
    output:
        "results/datasets/{dataset}/fastas/{sample}.{ref}.consensus.fa.gz.fai",
    log:
        "logs/{dataset}/samtools/faidx/{sample}.{ref}.consensus.log",
    benchmark:
        "benchmarks/{dataset}/samtools/faidx/{sample}.{ref}.consensus.log"
    wrapper:
        "v2.4.0/bio/samtools/faidx"


rule doAncError:
    """
    Calculates error rates per individual using an 'error free' individual and
    an outgroup. Currently, the outgroup is forced to be the reference, so this
    is only suitable if that is truly the case.
    """
    input:
        ref="results/ref/{ref}/{ref}.fa",
        fai="results/ref/{ref}/{ref}.fa.fai",
        errfree=expand(
            "results/datasets/{{dataset}}/fastas/{sample}.{{ref}}.consensus.fa.gz",
            sample=config["params"]["angsd"]["error_free_ind"],
        ),
        errfreefai=expand(
            "results/datasets/{{dataset}}/fastas/{sample}.{{ref}}.consensus.fa.gz.fai",
            sample=config["params"]["angsd"]["error_free_ind"],
        ),
        bam="results/datasets/{dataset}/bams/{sample}.{ref}.bam",
        bai="results/datasets/{dataset}/bams/{sample}.{ref}.bam.bai",
        sites="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.sites",
        idx="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.sites.idx",
    output:
        err="results/datasets/{dataset}/qc/doAncError/{sample}/{sample}.{ref}_{sites}-filts.ancError",
    log:
        "logs/{dataset}/angsd/doAncError/{sample}.{ref}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/angsd/doAncError/{sample}.{ref}_{sites}-filts.log"
    container:
        angsd_container
    resources:
        runtime=lambda w, attempt: attempt * 360,
    params:
        extra=config["params"]["angsd"]["extra"],
        mapQ=config["mapQ"],
        baseQ=config["baseQ"],
        out=lambda w, output: os.path.splitext(output.err)[0],
    shell:
        """
        angsd -doAncError 2 -i {input.bam} -anc {input.ref} -ref {input.errfree} \
            -nThreads {threads} {params.extra} -minMapQ {params.mapQ} \
            -minQ {params.baseQ} -sites {input.sites} -out {params.out} 2> {log}
        """


rule estError:
    """
    Estimates error rate for individuals from ANGSD outputs
    """
    input:
        "results/datasets/{dataset}/qc/doAncError/{sample}/{sample}.{ref}_{sites}-filts.ancError",
    output:
        table="results/datasets/{dataset}/qc/doAncError/{sample}/{sample}.{ref}_{sites}-filts.errorEst.txt",
        pdf="results/datasets/{dataset}/qc/doAncError/{sample}/{sample}.{ref}_{sites}-filts.errorEst.pdf",
    log:
        "logs/{dataset}/angsd/estError/{sample}.{ref}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/angsd/estError/{sample}.{ref}_{sites}-filts.log"
    conda:
        "../envs/r.yaml"
    shadow:
        "minimal"
    shell:
        """
        (cd results/datasets/{wildcards.dataset}/qc/doAncError/{wildcards.sample}
        Rscript \
            <(curl https://raw.githubusercontent.com/ANGSD/angsd/66a5961fcbf3b691cf39f96f4bd90868efa002ea/R/estError.R) \
            file={wildcards.sample}.{wildcards.ref}_{wildcards.sites}-filts.ancError
        
        for i in errorEst*; do
            mv $i {wildcards.sample}.{wildcards.ref}_{wildcards.sites}-filts.$i
        done

        sed -i 's/ind/{wildcards.sample}/g' \
            {wildcards.sample}.{wildcards.ref}_{wildcards.sites}-filts.errorEst.txt) \
            2> {log}
        """


rule cat_error:
    input:
        expand(
            "results/datasets/{{dataset}}/qc/doAncError/{sample}/{sample}.{{ref}}_{{sites}}-filts.errorEst.txt",
            sample=samples.index,
        ),
    output:
        est="results/datasets/{dataset}/qc/doAncError/{dataset}.{ref}_all_{sites}-filts.errorEst.tsv",
        overallest="results/datasets/{dataset}/qc/doAncError/{dataset}.{ref}_all_{sites}-filts.errorEstOverall.tsv",
    log:
        "logs/{dataset}/angsd/estError/{dataset}.{ref}_all_{sites}-filts_combine.log",
    benchmark:
        "benchmarks/{dataset}/angsd/estError/{dataset}.{ref}_all_{sites}-filts_combine.log"
    conda:
        "../envs/shell.yaml"
    shell:
        """
        (echo 'C->A    G->A    T->A    A->C    G->C    T->C    A->G    C->G    T->G    A->T    C->T    G->T' > {output.est}
        echo 'sample    error%' > {output.overallest}
        for i in {input}; do
            head -n2 $i | tail -n1 | tr -d '\"' >> {output.est}
            tail -n1 $i | tr ' ' '\t' | tr -d '\"' | cut -f1,2 >> {output.overallest}
        done) 2> {log}
        """
