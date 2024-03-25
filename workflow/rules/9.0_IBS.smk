# Estimation of pairwise identity by state genetic distance for all individuals


rule angsd_doIBS:
    """
    Estimates IBS matrix as well as position specific states for all individuals.
    """
    input:
        bamlist="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
        bams=get_bamlist_bams,
        bais=get_bamlist_bais,
        sites="results/datasets/{dataset}/filters/snps/{dataset}.{ref}_{population}{dp}_{sites}-filts_snps.sites",
        idx="results/datasets/{dataset}/filters/snps/{dataset}.{ref}_{population}{dp}_{sites}-filts_snps.sites.idx",
    output:
        ibs="results/datasets/{dataset}/analyses/IBS/{dataset}.{ref}_{population}{dp}_{sites}-filts.ibs.gz",
        ibsmat="results/datasets/{dataset}/analyses/IBS/{dataset}.{ref}_{population}{dp}_{sites}-filts.ibsMat",
        arg="results/datasets/{dataset}/analyses/IBS/{dataset}.{ref}_{population}{dp}_{sites}-filts.arg",
    log:
        "logs/{dataset}/angsd/doIBS/{dataset}.{ref}_{population}{dp}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/angsd/doIBS/{dataset}.{ref}_{population}{dp}_{sites}-filts.log"
    container:
        angsd_container
    params:
        mapQ=config["mapQ"],
        baseQ=config["baseQ"],
        ibs=config["params"]["ibs"]["doibs"],
        out=lambda w, output: os.path.splitext(output.arg)[0],
    threads: 8
    resources:
        runtime=lambda wildcards, attempt: attempt * 2880,
    shell:
        """
        angsd -doIBS {params.ibs} -bam {input.bamlist} -nThreads {threads} -doCounts 1 \
            -minMapQ {params.mapQ} -minQ {params.baseQ} -sites {input.sites} \
            -doMajorMinor 3 -makeMatrix 1 -out {params.out} &> {log}
        """
