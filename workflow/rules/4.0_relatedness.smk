# Pairwise individual relatedness with NGSrelate and  R0, R1, KING-robust kinship
# method from Waples et al. 2019, MolEcol


rule est_kinship_stats:
    """
    Uses the equations from Waples et al. 2019, MolEcol to estimate R0, R1, and KING-
    robust kinship between all sample pairings.
    """
    input:
        sfs="results/datasets/{dataset}/analyses/sfs/{dataset}.{ref}_{ind1}-{ind2}{dp}.sfs",
    output:
        "results/datasets/{dataset}/analyses/kinship/waples2019/{dataset}.{ref}_{ind1}-{ind2}{dp}.kinship",
    log:
        "logs/{dataset}/kinship/waples2019/{dataset}.{ref}_{ind1}-{ind2}{dp}_kinship.log",
    benchmark:
        "benchmarks/{dataset}/kinship/waples2019/{dataset}.{ref}_{ind1}-{ind2}{dp}_kinship.log"
    wildcard_constraints:
        ind1="|".join([i for i in samples.index.tolist()]),
        ind2="|".join([i for i in samples.index.tolist()]),
    conda:
        "../envs/r.yaml"
    resources:
        time=lambda wildcards, attempt: attempt * 15,
    script:
        "../scripts/kinship.R"


rule compile_kinship_stats:
    """
    Compiles kinship stats for all pairs into a single table.
    """
    input:
        get_kinship,
    output:
        "results/datasets/{dataset}/analyses/kinship/waples2019/{dataset}.{ref}_all{dp}.kinship",
    log:
        "logs/{dataset}/kinship/waples2019/{dataset}.{ref}_all{dp}_compile-stats.log",
    benchmark:
        "benchmarks/{dataset}/kinship/waples2019/{dataset}.{ref}_all{dp}_compile-stats.log"
    conda:
        "../envs/shell.yaml"
    resources:
        time=lambda wildcards, attempt: attempt * 15,
    shell:
        """
        (echo "ind1    ind2    R0    R1    KING" > {output}
        cat {input} >> {output}) 2> {log}
        """


rule kinship_table_html:
    """
    Converts kinship table to html for report.
    """
    input:
        "results/datasets/{dataset}/analyses/kinship/waples2019/{dataset}.{ref}_all{dp}.kinship",
    output:
        report(
            "results/datasets/{dataset}/analyses/kinship/waples2019/{dataset}.{ref}_all{dp}.kinship.html",
            category="Kinship",
            labels={"Topic": "Waples et al. 2019 R0,R1,KING", "Type": "Table"},
        ),
    log:
        "logs/{dataset}/kinship/waples2019/{dataset}.{ref}_all{dp}_tsv2html.log",
    benchmark:
        "benchmarks/{dataset}/kinship/waples2019/{dataset}.{ref}_all{dp}_tsv2html.log"
    conda:
        "../envs/r-rectable.yaml"
    script:
        "../scripts/tsv2html.Rmd"


rule ngsrelate:
    """
    Estimates inbreeding and relatedness measures using NGSrelate.
    """
    input:
        beagle="results/datasets/{dataset}/beagles/pruned/{dataset}.{ref}_all{dp}_pruned.beagle.gz",
        bamlist="results/datasets/{dataset}/bamlists/{dataset}.{ref}_all{dp}.bamlist",
        inds="results/datasets/{dataset}/poplists/{dataset}_all.indiv.list",
    output:
        relate="results/datasets/{dataset}/analyses/kinship/ngsrelate/{dataset}.{ref}_all{dp}_relate.tsv",
        samples="results/datasets/{dataset}/analyses/kinship/ngsrelate/{dataset}.{ref}_all{dp}_samples.list",
    log:
        "logs/{dataset}/kinship/ngsrelate/{dataset}.{ref}_all{dp}.log",
    container:
        ngsrelate_container
    threads: lambda wildcards, attempt: attempt * 4
    params:
        nind=lambda w: len(get_samples_from_pop("all")),
    resources:
        time=lambda wildcards, attempt: attempt * 360,
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
        "results/datasets/{dataset}/analyses/kinship/ngsrelate/{dataset}.{ref}_all{dp}_relate.tsv",
    output:
        report(
            "results/datasets/{dataset}/analyses/kinship/ngsrelate/{dataset}.{ref}_all{dp}_relate.html",
            category="Kinship",
            labels={"Topic": "NgsRelate", "Type": "Table"},
        ),
    log:
        "logs/{dataset}/kinship/ngsrelate/{dataset}.{ref}_all{dp}_tsv2html.log",
    conda:
        "../envs/r-rectable.yaml"
    script:
        "../scripts/tsv2html.Rmd"