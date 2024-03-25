# Rules for estimating pairwise population Fst


rule realSFS_fst_index:
    """
    Generates Fst index file for producing Fst estimates.
    """
    input:
        saf1="results/datasets/{dataset}/safs/{dataset}.{ref}_{population1}{dp}_{sites}-filts.saf.idx",
        saf1_others=multiext(
            "results/datasets/{dataset}/safs/{dataset}.{ref}_{population1}{dp}_{sites}-filts.saf",
            ".pos.gz",
            ".gz",
        ),
        saf2="results/datasets/{dataset}/safs/{dataset}.{ref}_{population2}{dp}_{sites}-filts.saf.idx",
        saf2_others=multiext(
            "results/datasets/{dataset}/safs/{dataset}.{ref}_{population2}{dp}_{sites}-filts.saf",
            ".pos.gz",
            ".gz",
        ),
        sfs="results/datasets/{dataset}/analyses/sfs/{dataset}.{ref}_{population1}-{population2}{dp}_{sites}-filts.sfs",
    output:
        fstidx=temp(
            "results/datasets/{dataset}/analyses/fst/{dataset}.{ref}_{population1}-{population2}{dp}_{sites}-filts.fst.idx"
        ),
        fstgz=temp(
            "results/datasets/{dataset}/analyses/fst/{dataset}.{ref}_{population1}-{population2}{dp}_{sites}-filts.fst.gz"
        ),
    container:
        angsd_container
    log:
        "logs/{dataset}/realSFS/fst/index/{dataset}.{ref}_{population1}-{population2}{dp}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/realSFS/fst/index/{dataset}.{ref}_{population1}-{population2}{dp}_{sites}-filts.log"
    params:
        out=lambda w, output: output.fstidx.removesuffix(".fst.idx"),
        fst=config["params"]["fst"]["whichFst"],
    resources:
        runtime=lambda wildcards, attempt: attempt * 120,
    shell:
        """
        realSFS fst index -whichFst {params.fst} \
            {input.saf1} {input.saf2} -sfs {input.sfs} \
            -fstout {params.out} &> {log}
        """


rule realSFS_fst_stats:
    """
    Estimates global weighted and unweighted Fst for each population pair.
    """
    input:
        fstidx="results/datasets/{dataset}/analyses/fst/{dataset}.{ref}_{population1}-{population2}{dp}_{sites}-filts.fst.idx",
        fstgz="results/datasets/{dataset}/analyses/fst/{dataset}.{ref}_{population1}-{population2}{dp}_{sites}-filts.fst.gz",
    output:
        fstglob="results/datasets/{dataset}/analyses/fst/{dataset}.{ref}_{population1}-{population2}{dp}_{sites}-filts.fst.global.tsv",
    container:
        angsd_container
    log:
        "logs/{dataset}/realSFS/fst/stats/{dataset}.{ref}_{population1}-{population2}{dp}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/realSFS/fst/stats/{dataset}.{ref}_{population1}-{population2}{dp}_{sites}-filts.log"
    resources:
        runtime=lambda wildcards, attempt: attempt * 60,
    shell:
        r"""
        realSFS fst stats {input.fstidx} | \
            awk '{{print "{wildcards.population1}\t{wildcards.population2}\t"\
            $1"\t"$2}}' > {output.fstglob} 2> {log}
        """


rule realSFS_fst_stats2:
    """
    Estimates fst in sliding windows across the genome.
    """
    input:
        fstidx="results/datasets/{dataset}/analyses/fst/{dataset}.{ref}_{population1}-{population2}{dp}_{sites}-filts.fst.idx",
        fstgz="results/datasets/{dataset}/analyses/fst/{dataset}.{ref}_{population1}-{population2}{dp}_{sites}-filts.fst.gz",
    output:
        fstwin="results/datasets/{dataset}/analyses/fst/{dataset}.{ref}_{population1}-{population2}{dp}_{sites}-filts.fst.window_{win}_{step}.tsv",
    container:
        angsd_container
    log:
        "logs/{dataset}/realSFS/fst/stats2/{dataset}.{ref}_{population1}-{population2}{dp}_{sites}-filts.window_{win}_{step}.log",
    benchmark:
        "benchmarks/{dataset}/realSFS/fst/stats2/{dataset}.{ref}_{population1}-{population2}{dp}_{sites}-filts.window_{win}_{step}.log"
    shell:
        r"""
        realSFS fst stats2 {input.fstidx} -win {wildcards.win} -step {wildcards.step} \
            -type 0 | awk '{{print "{wildcards.population1}\t{wildcards.population2}\t"\
            $0}}' > {output.fstwin} 2> {log}
        """


rule aggregate_fst_global:
    """
    Aggregates global fst estimates for all population pairs in one table.
    """
    input:
        unpack(get_fst),
    output:
        glob="results/datasets/{dataset}/analyses/fst/{dataset}.{ref}_{unit}pairs{dp}_{sites}-filts.fst.{scale}.tsv",
    log:
        "logs/{dataset}/realSFS/fst/aggregate/{dataset}.{ref}_{unit}pairs{dp}_{sites}-filts.{scale}.log",
    benchmark:
        "benchmarks/{dataset}/realSFS/fst/aggregate/{dataset}.{ref}_{unit}pairs{dp}_{sites}-filts.{scale}.log"
    conda:
        "../envs/shell.yaml"
    wildcard_constraints:
        unit="ind|pop",
        scale="global",
    shell:
        """
        (printf "pop1\tpop2\tunweight.fst\tweight.fst\n" > {output.glob}
        cat {input} >> {output.glob}) 2> {log}
        """


rule aggregate_fst_window:
    """
    Aggregates fst windows for all population pairs in one table.
    """
    input:
        unpack(get_fst),
    output:
        window="results/datasets/{dataset}/analyses/fst/{dataset}.{ref}_{unit}pairs{dp}_{sites}-filts.fst.{scale}_{win}_{step}.tsv",
    log:
        "logs/{dataset}/realSFS/fst/aggregate/{dataset}.{ref}_{unit}pairs{dp}_{sites}-filts.{scale}_{win}_{step}.log",
    benchmark:
        "benchmarks/{dataset}/realSFS/fst/aggregate/{dataset}.{ref}_{unit}pairs{dp}_{sites}-filts.{scale}_{win}_{step}.log"
    conda:
        "../envs/shell.yaml"
    wildcard_constraints:
        unit="ind|pop",
        scale="window",
    shell:
        """
        (printf "pop1\tpop2\twindow\tchr\twindow_center\tNsites\tweight.fst\n" \
            > {output.window}
        for i in {input}; do
            tail -n +2 $i >> {output.window}
        done) 2> {log}
        """


rule plot_fst:
    """
    Plots fst heatmap of global fst estimates for all population pairs.
    """
    input:
        "results/datasets/{dataset}/analyses/fst/{dataset}.{ref}_{unit}pairs{dp}_{sites}-filts.fst.global.tsv",
    output:
        report(
            "results/datasets/{dataset}/plots/fst/{dataset}.{ref}_{unit}pairs{dp}_{sites}-filts.fst.global.pdf",
            category="05 Fst",
            subcategory="Global",
            labels=lambda w: {
                "Filter": "{sites}",
                **dp_report(w),
                **unit_report(w),
                "Type": "Heatmap",
            },
        ),
    log:
        "logs/{dataset}/realSFS/fst/plot/{dataset}.{ref}_{unit}pairs{dp}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/realSFS/fst/plot/{dataset}.{ref}_{unit}pairs{dp}_{sites}-filts.log"
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/plot_fst.R"
