# Rules for estimating diversity and neutrality statistics


rule realSFS_saf2theta:
    """
    Generates a thetas index for each population from SAF and SFS.
    """
    input:
        safidx="results/datasets/{dataset}/safs/{dataset}.{ref}_{population}{dp}_{sites}-filts.saf.idx",
        saf_others=multiext(
            "results/datasets/{dataset}/safs/{dataset}.{ref}_{population}{dp}_{sites}-filts.saf",
            ".pos.gz",
            ".gz",
        ),
        sfs="results/datasets/{dataset}/analyses/sfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.sfs",
    output:
        thetasidx=temp(
            "results/datasets/{dataset}/analyses/thetas/{dataset}.{ref}_{population}{dp}_{sites}-filts.thetas.idx"
        ),
        thetas=temp(
            "results/datasets/{dataset}/analyses/thetas/{dataset}.{ref}_{population}{dp}_{sites}-filts.thetas.gz"
        ),
    container:
        angsd_container
    log:
        "logs/{dataset}/realSFS/saf2theta/{dataset}.{ref}_{population}{dp}_{sites}-filts.log",
    params:
        out=lambda w, output: output.thetas.removesuffix(".thetas.gz"),
        fold=config["params"]["realsfs"]["fold"],
    resources:
        runtime=lambda wildcards, attempt: attempt * 120,
    shell:
        """
        realSFS saf2theta {input.safidx} -sfs {input.sfs} -fold {params.fold} \
            -outname {params.out} &> {log}
        """


rule thetaStat:
    """
    Estimates diversity and neutrality statistics in sliding windows of a set size.
    """
    input:
        thetasidx="results/datasets/{dataset}/analyses/thetas/{dataset}.{ref}_{population}{dp}_{sites}-filts.thetas.idx",
        thetas="results/datasets/{dataset}/analyses/thetas/{dataset}.{ref}_{population}{dp}_{sites}-filts.thetas.gz",
    output:
        thetas="results/datasets/{dataset}/analyses/thetas/{dataset}.{ref}_{population}{dp}_{sites}-filts.thetaWindows.{win}_{step}.pestPG",
    container:
        angsd_container
    log:
        "logs/{dataset}/thetaStat/{dataset}.{ref}_{population}{dp}_{sites}-filts.{win}_{step}.log",
    params:
        out=lambda w, output: os.path.splitext(output.thetas)[0],
    resources:
        runtime=lambda wildcards, attempt: attempt * 120,
    shell:
        """
        thetaStat do_stat {input.thetasidx} -win {wildcards.win} -type 0 \
            -step {wildcards.step} -outnames {params.out} &> {log}
        """


rule plot_thetas:
    """
    Plots per population distributions of theta, pi, and Tajima's D across the genome.
    """
    input:
        expand(
            "results/datasets/{{dataset}}/analyses/thetas/{{dataset}}.{{ref}}_{population}{{dp}}_{{sites}}-filts.thetaWindows.{{win}}_{{step}}.pestPG",
            population=pop_list,
        ),
    output:
        watterson=report(
            "results/datasets/{dataset}/plots/thetas/{dataset}.{ref}_all{dp}_{sites}-filts.window_{win}_{step}.density.watterson.pdf",
            category="04.1 Watterson's Theta",
            labels=lambda w: {
                "Filter": "{sites}",
                **dp_report(w),
                "Win size": "{win}bp",
                "Win step": "{step}bp",
                "Type": "Violin Plot",
            },
        ),
        pi=report(
            "results/datasets/{dataset}/plots/thetas/{dataset}.{ref}_all{dp}_{sites}-filts.window_{win}_{step}.density.pi.pdf",
            category="04.2 Nucleotide Diversity (Pi)",
            labels=lambda w: {
                "Filter": "{sites}",
                **dp_report(w),
                "Win size": "{win}bp",
                "Win step": "{step}bp",
                "Type": "Violin Plot",
            },
        ),
        tajima=report(
            "results/datasets/{dataset}/plots/thetas/{dataset}.{ref}_all{dp}_{sites}-filts.window_{win}_{step}.density.tajima.pdf",
            category="04.3 Tajima's D",
            labels=lambda w: {
                "Filter": "{sites}",
                **dp_report(w),
                "Win size": "{win}bp",
                "Win step": "{step}bp",
                "Type": "Violin Plot",
            },
        ),
        tables=expand(
            "results/datasets/{{dataset}}/analyses/thetas/{{dataset}}.{{ref}}_all{{dp}}_{{sites}}-filts.window_{{win}}_{{step}}.{stat}.mean.tsv",
            stat=["watterson", "tajima", "pi"],
        ),
    log:
        "logs/{dataset}/thetaStat/{dataset}.{ref}_all{dp}_{sites}-filts.{win}_{step}.plot.log",
    benchmark:
        "benchmarks/{dataset}/thetaStat/{dataset}.{ref}_all{dp}_{sites}-filts.{win}_{step}.plot.log"
    conda:
        "../envs/r.yaml"
    params:
        popnames=pop_list,
        plotpre=lambda w, output: output["watterson"].removesuffix(".watterson.pdf"),
        tabpre=lambda w, output: output["tables"][0].removesuffix(".watterson.mean.tsv"),
        minsites=config["params"]["thetas"]["minsites"],
    script:
        "../scripts/plot_thetas.R"


rule theta_tables:
    """
    Converts Theta table to html.
    """
    input:
        "results/datasets/{dataset}/analyses/thetas/{dataset}.{ref}_all{dp}_{sites}-filts.window_{win}_{step}.{stat}.mean.tsv",
    output:
        report(
            "results/datasets/{dataset}/analyses/thetas/{dataset}.{ref}_all{dp}_{sites}-filts.window_{win}_{step}.{stat}.mean.html",
            category=lambda w: theta_report(w),
            labels=lambda w: {
                "Filter": "{sites}",
                **dp_report(w),
                "Win size": "{win}bp",
                "Win step": "{step}bp",
                "Type": "Table (Means)",
            },
        ),
    log:
        "logs/{dataset}/thetaStat/{dataset}.{ref}_all{dp}_{sites}-filts.{win}_{step}.{stat}_tsv2html.log",
    conda:
        "../envs/r-rectable.yaml"
    script:
        "../scripts/tsv2html.Rmd"
