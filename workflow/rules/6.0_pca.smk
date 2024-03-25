# Rules for performing principal component analysis


rule remove_excl_pca_admix:
    """
    Excludes requested individuals from beagle file only for PCA and Admix. This is 
    largely to exclude close relatives from these analyses (where they can impact 
    results), while allowing them in all other analyses.
    """
    input:
        "results/datasets/{dataset}/beagles/pruned/{dataset}.{ref}_all{dp}_{sites}-filts.pruned_maxkbdist-{maxkb}_minr2-{r2}.beagle.gz",
    output:
        "results/datasets/{dataset}/beagles/pruned/{dataset}.{ref}_all_excl_pca-admix{dp}_{sites}-filts.pruned_maxkbdist-{maxkb}_minr2-{r2}.beagle.gz",
    log:
        "logs/{dataset}/ngsLD/excl_pca_admix_beagle/{dataset}.{ref}_all_excl_pca-admix{dp}_{sites}-filts.pruned_maxkbdist-{maxkb}_minr2-{r2}.log",
    benchmark:
        "benchmarks/{dataset}/ngsLD/excl_pca_admix_beagle/{dataset}.{ref}_all_excl_pca-admix{dp}_{sites}-filts.pruned_maxkbdist-{maxkb}_minr2-{r2}.log"
    conda:
        "../envs/shell.yaml"
    params:
        remcols=get_excl_ind_cols,
    shell:
        """
        zcat {input} | cut -f{params.remcols} --complement | gzip > {output} 2> {log}
        """


rule pca_pcangsd:
    """
    Produces covariance matrix from SNP genotype likelihood data with PCAngsd.
    """
    input:
        beagle=expand(
            "results/datasets/{{dataset}}/beagles/pruned/{{dataset}}.{{ref}}_{{population}}{{dp}}_{{sites}}-filts.pruned_maxkbdist-{maxkb}_minr2-{r2}.beagle.gz",
            maxkb=config["params"]["ngsld"]["max_kb_dist_pruning_dataset"],
            r2=config["params"]["ngsld"]["pruning_min-weight_dataset"],
        ),
    output:
        cov="results/datasets/{dataset}/analyses/pcangsd/{dataset}.{ref}_{population}{dp}_{sites}-filts.cov",
    log:
        "logs/{dataset}/pcangsd/{dataset}.{ref}_{population}{dp}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/pcangsd/{dataset}.{ref}_{population}{dp}_{sites}-filts.log"
    container:
        pcangsd_container
    params:
        prefix=lambda w, output: os.path.splitext(output.cov)[0],
    threads: lambda wildcards, attempt: attempt
    resources:
        runtime=lambda wildcards, attempt: attempt * 60,
    shell:
        """
        pcangsd -b {input.beagle} -o {params.prefix} &> {log}
        """


rule plot_pca:
    """
    Plots PCA for various pairs of PCs.
    """
    input:
        "results/datasets/{dataset}/analyses/pcangsd/{dataset}.{ref}_{population}{dp}_{sites}-filts.cov",
        "results/datasets/{dataset}/poplists/{dataset}_{population}.indiv.list",
    output:
        report(
            "results/datasets/{dataset}/plots/pca/{dataset}.{ref}_{population}{dp}_{sites}-filts_pc{xpc}-{ypc}.svg",
            category="03.1 PCA",
            labels=lambda w: {
                "Filter": "{sites}",
                **dp_report(w),
                "PCs": "PC{xpc}-PC{ypc}",
                "Type": "Scatterplot",
            },
        ),
    log:
        "logs/{dataset}/pcangsd/{dataset}.{ref}_{population}{dp}_{sites}-filts_pc{xpc}-{ypc}_plot.log",
    benchmark:
        "benchmarks/{dataset}/pcangsd/{dataset}.{ref}_{population}{dp}_{sites}-filts_pc{xpc}-{ypc}_plot.log"
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/plot_PCA.R"
