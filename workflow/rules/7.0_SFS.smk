# Rules for generating site frequency spectra


rule realSFS_1dSFS:
    """
    Generate a 1D site frequency spectrum.
    """
    input:
        saf="results/datasets/{dataset}/safs/{dataset}.{ref}_{population}{dp}_{sites}-filts.saf.idx",
        others=multiext(
            "results/datasets/{dataset}/safs/{dataset}.{ref}_{population}{dp}_{sites}-filts.saf",
            ".pos.gz",
            ".gz",
        ),
    output:
        sfs=ensure(
            "results/datasets/{dataset}/analyses/sfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.sfs",
            non_empty=True,
        ),
    container:
        angsd_container
    log:
        "logs/{dataset}/realSFS/1dSFS/{dataset}.{ref}_{population}{dp}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/realSFS/1dSFS/{dataset}.{ref}_{population}{dp}_{sites}-filts.log"
    params:
        fold=config["params"]["realsfs"]["fold"],
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        runtime=lambda wildcards, attempt: attempt * 720,
    shell:
        """
        realSFS {input.saf} -fold {params.fold} -P {threads} \
            > {output.sfs} 2> {log}
        """


rule realSFS_1dSFS_bootstrap:
    """
    Generate a 1D site frequency spectrum.
    """
    input:
        saf="results/datasets/{dataset}/safs/{dataset}.{ref}_{population}{dp}_{sites}-filts.saf.idx",
        others=multiext(
            "results/datasets/{dataset}/safs/{dataset}.{ref}_{population}{dp}_{sites}-filts.saf",
            ".pos.gz",
            ".gz",
        ),
    output:
        sfs=ensure(
            "results/datasets/{dataset}/analyses/sfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.boot.sfs",
            non_empty=True,
        ),
    container:
        angsd_container
    log:
        "logs/{dataset}/realSFS/1dSFS/{dataset}.{ref}_{population}{dp}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/realSFS/1dSFS/{dataset}.{ref}_{population}{dp}_{sites}-filts.log"
    params:
        fold=config["params"]["realsfs"]["fold"],
        boot=config["params"]["realsfs"]["sfsboot"],
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        runtime="7d",
    shell:
        """
        realSFS {input.saf} -fold {params.fold} -P {threads} \
            -bootstrap {params.boot} > {output.sfs} 2> {log}
        """


rule realSFS_2dSFS:
    """
    Generate a 2D site frequency spectrum.
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
    output:
        sfs=ensure(
            "results/datasets/{dataset}/analyses/sfs/{dataset}.{ref}_{population1}-{population2}{dp}_{sites}-filts.sfs",
            non_empty=True,
        ),
    container:
        angsd_container
    log:
        "logs/{dataset}/realSFS/2dSFS/{dataset}.{ref}_{population1}-{population2}{dp}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/realSFS/2dSFS/{dataset}.{ref}_{population1}-{population2}{dp}_{sites}-filts.log"
    wildcard_constraints:
        population1="|".join(
            [i for i in samples.index.tolist()]
            + [i for i in samples.population.values.tolist()]
        ),
        population2="|".join(
            [i for i in samples.index.tolist()]
            + [i for i in samples.population.values.tolist()]
        ),
    params:
        fold=config["params"]["realsfs"]["fold"],
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        runtime=lambda wildcards, attempt: attempt * 720,
    shell:
        """
        realSFS {input.saf1} {input.saf2} -fold {params.fold} \
            -P {threads} > {output.sfs} 2> {log}
        """


rule realSFS_2dSFS_bootstrap:
    """
    Generate a 2D site frequency spectrum.
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
    output:
        sfs=ensure(
            "results/datasets/{dataset}/analyses/sfs/{dataset}.{ref}_{population1}-{population2}{dp}_{sites}-filts.boot.sfs",
            non_empty=True,
        ),
    container:
        angsd_container
    log:
        "logs/{dataset}/realSFS/2dSFS/{dataset}.{ref}_{population1}-{population2}{dp}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/realSFS/2dSFS/{dataset}.{ref}_{population1}-{population2}{dp}_{sites}-filts.log"
    wildcard_constraints:
        population1="|".join(
            [i for i in samples.index.tolist()]
            + [i for i in samples.population.values.tolist()]
        ),
        population2="|".join(
            [i for i in samples.index.tolist()]
            + [i for i in samples.population.values.tolist()]
        ),
    params:
        fold=config["params"]["realsfs"]["fold"],
        boot=config["params"]["realsfs"]["sfsboot"],
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        runtime=lambda wildcards, attempt: attempt * 720,
    shell:
        """
        realSFS {input.saf1} {input.saf2} -fold {params.fold} \
            -bootstrap {params.boot} -P {threads} > {output.sfs} 2> {log}
        """
