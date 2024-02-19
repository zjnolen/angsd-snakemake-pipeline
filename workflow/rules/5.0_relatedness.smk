# Pairwise individual relatedness with NGSrelate and  R0, R1, KING-robust kinship
# method from Waples et al. 2019, MolEcol


rule est_kinship_stats_sfs:
    """
    Uses the equations from Waples et al. 2019, MolEcol to estimate R0, R1, and KING-
    robust kinship between all sample pairings.
    """
    input:
        sfs="results/datasets/{dataset}/analyses/sfs/{dataset}.{ref}_{ind1}-{ind2}{dp}_{sites}-filts.sfs",
    output:
        "results/datasets/{dataset}/analyses/kinship/ibsrelate_sfs/{dataset}.{ref}_{ind1}-{ind2}{dp}_{sites}-filts.kinship",
    log:
        "logs/{dataset}/kinship/ibsrelate_sfs/{dataset}.{ref}_{ind1}-{ind2}{dp}_{sites}-filts_kinship.log",
    benchmark:
        "benchmarks/{dataset}/kinship/ibsrelate_sfs/{dataset}.{ref}_{ind1}-{ind2}{dp}_{sites}-filts_kinship.log"
    wildcard_constraints:
        ind1="|".join([i for i in samples.index.tolist()]),
        ind2="|".join([i for i in samples.index.tolist()]),
    conda:
        "../envs/r.yaml"
    resources:
        runtime=lambda wildcards, attempt: attempt * 15,
    script:
        "../scripts/kinship_sfs.R"


rule compile_kinship_stats_sfs:
    """
    Compiles kinship stats for all pairs into a single table.
    """
    input:
        get_kinship,
    output:
        "results/datasets/{dataset}/analyses/kinship/ibsrelate_sfs/{dataset}.{ref}_all{dp}_{sites}-filts.kinship",
    log:
        "logs/{dataset}/kinship/ibsrelate_sfs/{dataset}.{ref}_all{dp}_{sites}-filts_compile-stats.log",
    benchmark:
        "benchmarks/{dataset}/kinship/ibsrelate_sfs/{dataset}.{ref}_all{dp}_{sites}-filts_compile-stats.log"
    conda:
        "../envs/shell.yaml"
    resources:
        runtime=lambda wildcards, attempt: attempt * 15,
    shell:
        """
        (printf "ind1\tind2\tR0\tR1\tKING\n" > {output}
        cat {input} >> {output}) 2> {log}
        """


rule doGlf1_ibsrelate:
    input:
        unpack(filt_depth),
        bam="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
        bams=get_bamlist_bams,
        bais=get_bamlist_bais,
        ref="results/ref/{ref}/{ref}.fa",
        regions="results/datasets/{dataset}/filters/chunks/{ref}_chunk{chunk}.rf",
    output:
        glf=temp(
            "results/datasets/{dataset}/glfs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.glf.gz"
        ),
        arg="results/datasets/{dataset}/glfs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.arg",
    log:
        "logs/{dataset}/angsd/doGlf1/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/angsd/doGlf1/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.log"
    container:
        angsd_container
    params:
        gl_model=config["params"]["angsd"]["gl_model"],
        extra=config["params"]["angsd"]["extra"],
        mapQ=config["mapQ"],
        baseQ=config["baseQ"],
        out=lambda w, output: os.path.splitext(output.arg)[0],
    resources:
        runtime=lambda wildcards, attempt: attempt * 360,
    threads: lambda wildcards, attempt: attempt * 2
    shell:
        """
        angsd -doGlf 1 -bam {input.bam} -GL {params.gl_model} -ref {input.ref} \
            -nThreads {threads} {params.extra} -minMapQ {params.mapQ} \
            -minQ {params.baseQ} -sites {input.sites} -rf {input.regions} \
            -out {params.out} &> {log}
        """


rule ibsrelate:
    input:
        "results/datasets/{dataset}/glfs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.glf.gz",
    output:
        "results/datasets/{dataset}/analyses/kinship/ibsrelate_ibs/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.ibspair",
    log:
        "logs/{dataset}/angsd/ibsrelate_ibs/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/angsd/ibsrelate_ibs/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.log"
    container:
        angsd_container
    params:
        nind=lambda w: len(get_samples_from_pop(w.population)),
        maxsites=config["chunk_size"],
        out=lambda w, output: os.path.splitext(output[0])[0],
    resources:
        runtime="7d",
    threads: lambda wildcards, attempt: attempt * 10
    shell:
        """
        ibs -glf {input} -model 0 -nInd {params.nind} -allpairs 1 \
            -m {params.maxsites} -outFileName {params.out}
        """


rule est_kinship_stats_ibs:
    """
    Uses the equations from Waples et al. 2019, MolEcol to estimate R0, R1, and KING-
    robust kinship between all sample pairings.
    """
    input:
        ibs=expand(
            "results/datasets/{{dataset}}/analyses/kinship/ibsrelate_ibs/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}_{{sites}}-filts.ibspair",
            chunk=chunklist,
        ),
        inds="results/datasets/{dataset}/poplists/{dataset}_all.indiv.list",
    output:
        "results/datasets/{dataset}/analyses/kinship/ibsrelate_ibs/{dataset}.{ref}_{population}{dp}_{sites}-filts.kinship",
    log:
        "logs/{dataset}/kinship/ibsrelate_ibs/{dataset}.{ref}_{population}{dp}_{sites}-filts_kinship.log",
    benchmark:
        "benchmarks/{dataset}/kinship/ibsrelate_ibs/{dataset}.{ref}_{population}{dp}_{sites}-filts_kinship.log"
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/kinship_ibs.R"


rule kinship_table_html:
    """
    Converts kinship table to html for report.
    """
    input:
        "results/datasets/{dataset}/analyses/kinship/ibsrelate_{type}/{dataset}.{ref}_all{dp}_{sites}-filts.kinship",
    output:
        report(
            "results/datasets/{dataset}/analyses/kinship/ibsrelate_{type}/{dataset}.{ref}_all{dp}_{sites}-filts.kinship.html",
            category="Relatedness",
            subcategory="IBSrelate - {type}",
            labels=lambda w: {"Filter": "{sites}", **dp_report(w), "Type": "Table"},
        ),
    log:
        "logs/{dataset}/kinship/ibsrelate_{type}/{dataset}.{ref}_all{dp}_{sites}-filts_tsv2html.log",
    benchmark:
        "benchmarks/{dataset}/kinship/ibsrelate_{type}/{dataset}.{ref}_all{dp}_{sites}-filts_tsv2html.log"
    conda:
        "../envs/r-rectable.yaml"
    script:
        "../scripts/tsv2html.Rmd"


rule ngsrelate:
    """
    Estimates inbreeding and relatedness measures using NGSrelate.
    """
    input:
        beagle=expand(
            "results/datasets/{{dataset}}/beagles/pruned/{{dataset}}.{{ref}}_all{{dp}}_{{sites}}-filts.pruned_maxkbdist-{maxkb}_minr2-{r2}.beagle.gz",
            maxkb=config["params"]["ngsld"]["max_kb_dist_pruning_dataset"],
            r2=config["params"]["ngsld"]["pruning_min-weight_dataset"],
        ),
        inds="results/datasets/{dataset}/poplists/{dataset}_all.indiv.list",
    output:
        relate="results/datasets/{dataset}/analyses/kinship/ngsrelate/{dataset}.{ref}_all{dp}_{sites}-filts_relate.tsv",
        samples="results/datasets/{dataset}/analyses/kinship/ngsrelate/{dataset}.{ref}_all{dp}_{sites}-filts_samples.list",
    log:
        "logs/{dataset}/kinship/ngsrelate/{dataset}.{ref}_all{dp}_{sites}-filts.log",
    container:
        ngsrelate_container
    threads: lambda wildcards, attempt: attempt * 4
    params:
        nind=lambda w: len(get_samples_from_pop("all")),
    resources:
        runtime=lambda wildcards, attempt: attempt * 360,
    shell:
        r"""
        (nsites=$(zcat {input.beagle} | tail -n +2 | wc -l)
        echo "nsites nind"
        echo $nsites {params.nind}
        cut -f1 {input.inds} | tail -n +2 > {output.samples}
        ngsRelate -G {input.beagle} -n {params.nind} -L $nsites -O {output.relate} \
            -z {output.samples}) &> {log}
        """


rule ngsrelate_summary:
    """
    Converts NGSrelate table to html.
    """
    input:
        "results/datasets/{dataset}/analyses/kinship/ngsrelate/{dataset}.{ref}_all{dp}_{sites}-filts_relate.tsv",
    output:
        report(
            "results/datasets/{dataset}/analyses/kinship/ngsrelate/{dataset}.{ref}_all{dp}_{sites}-filts_relate.html",
            category="Relatedness",
            subcategory="NgsRelate",
            labels=lambda w: {"Filter": "{sites}", **dp_report(w), "Type": "Table"},
        ),
    log:
        "logs/{dataset}/kinship/ngsrelate/{dataset}.{ref}_all{dp}_{sites}-filts_tsv2html.log",
    conda:
        "../envs/r-rectable.yaml"
    script:
        "../scripts/tsv2html.Rmd"
