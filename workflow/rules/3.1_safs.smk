# Rules for creating site allele frequency files used by ANGSD


rule angsd_doSaf_pop:
    """
    Generate a site allele frequency file for a given population and genome chunk.
    """
    input:
        unpack(filt_depth),
        unpack(get_anc_ref),
        bam="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
        bams=get_bamlist_bams,
        bais=get_bamlist_bais,
        ref="results/ref/{ref}/{ref}.fa",
        reffai="results/ref/{ref}/{ref}.fa.fai",
        regions="results/datasets/{dataset}/filters/chunks/{ref}_chunk{chunk}.rf",
    output:
        saf=temp(
            "results/datasets/{dataset}/safs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.saf.gz"
        ),
        safidx=temp(
            "results/datasets/{dataset}/safs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.saf.idx"
        ),
        safpos=temp(
            "results/datasets/{dataset}/safs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.saf.pos.gz"
        ),
        arg="results/datasets/{dataset}/safs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.arg",
    log:
        "logs/{dataset}/angsd/doSaf/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/angsd/doSaf/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.log"
    container:
        angsd_container
    wildcard_constraints:
        population="|".join(
            ["all"]
            + ["all_excl_pca-admix"]
            + [i for i in samples.population.values.tolist()]
            + [i for i in samples.depth.values.tolist()]
        ),
    params:
        gl_model=config["params"]["angsd"]["gl_model"],
        extra=config["params"]["angsd"]["extra"],
        extra_saf=config["params"]["angsd"]["extra_saf"],
        mapQ=config["mapQ"],
        baseQ=config["baseQ"],
        out=lambda w, output: os.path.splitext(output.arg)[0],
    resources:
        runtime=lambda wildcards, attempt: attempt * 180,
    threads: lambda wildcards, attempt: attempt * 2
    shell:
        """
        angsd -doSaf 1 -bam {input.bam} -GL {params.gl_model} -ref {input.ref} \
            -nThreads {threads} {params.extra} -minMapQ {params.mapQ} \
            -minQ {params.baseQ} -sites {input.sites} -anc {input.anc} \
            {params.extra_saf} -rf {input.regions} -out {params.out} &> {log}
        """


rule realSFS_catsaf:
    """
    Merge genome chunks into one file per population.
    """
    input:
        safs=lambda w: expand(
            "results/datasets/{{dataset}}/safs/chunks/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}_{{sites}}-filts.saf.idx",
            chunk=chunklist,
        ),
        safgz=lambda w: expand(
            "results/datasets/{{dataset}}/safs/chunks/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}_{{sites}}-filts.saf.gz",
            chunk=chunklist,
        ),
        safposgz=lambda w: expand(
            "results/datasets/{{dataset}}/safs/chunks/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}_{{sites}}-filts.saf.pos.gz",
            chunk=chunklist,
        ),
    output:
        multiext(
            "results/datasets/{dataset}/safs/{dataset}.{ref}_{population}{dp}_{sites}-filts.saf",
            ".idx",
            ".pos.gz",
            ".gz",
        ),
    log:
        "logs/{dataset}/realSFS/cat/{dataset}.{ref}_{population}{dp}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/realSFS/cat/{dataset}.{ref}_{population}{dp}_{sites}-filts.log"
    container:
        angsd_container
    wildcard_constraints:
        population="|".join(
            ["all"]
            + ["all_excl_pca-admix"]
            + [i for i in samples.population.values.tolist()]
            + [i for i in samples.depth.values.tolist()]
        ),
    params:
        out=lambda w, output: output[0].removesuffix(".saf.idx"),
    resources:
        runtime=lambda wildcards, attempt: attempt * 60,
    shell:
        """
        realSFS cat {input.safs} -P 1 -outnames {params.out} 2> {log}
        """


rule angsd_doSaf_sample:
    """
    Generate a site allele frequency file for a given subsampled population and genome
    chunk.
    """
    input:
        unpack(filt_depth),
        unpack(get_anc_ref),
        bam="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
        bams=get_bamlist_bams,
        bais=get_bamlist_bais,
        ref="results/ref/{ref}/{ref}.fa",
        reffai="results/ref/{ref}/{ref}.fa.fai",
    output:
        saf="results/datasets/{dataset}/safs/{dataset}.{ref}_{population}{dp}_{sites}-filts.saf.gz",
        safidx="results/datasets/{dataset}/safs/{dataset}.{ref}_{population}{dp}_{sites}-filts.saf.idx",
        safpos="results/datasets/{dataset}/safs/{dataset}.{ref}_{population}{dp}_{sites}-filts.saf.pos.gz",
        arg="results/datasets/{dataset}/safs/{dataset}.{ref}_{population}{dp}_{sites}-filts.arg",
    log:
        "logs/{dataset}/angsd/doSaf/{dataset}.{ref}_{population}{dp}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/angsd/doSaf/{dataset}.{ref}_{population}{dp}_{sites}-filts.log"
    container:
        angsd_container
    wildcard_constraints:
        population="|".join([i for i in samples.index.tolist()]),
    params:
        gl_model=config["params"]["angsd"]["gl_model"],
        extra=config["params"]["angsd"]["extra"],
        extra_saf=config["params"]["angsd"]["extra_saf"],
        mindepthind=config["params"]["angsd"]["mindepthind_heterozygosity"],
        mapQ=config["mapQ"],
        baseQ=config["baseQ"],
        out=lambda w, output: os.path.splitext(output.arg)[0],
    resources:
        runtime="120m",
    shell:
        """
        (angsd -doSaf 1 -bam {input.bam} -GL {params.gl_model} -ref {input.ref} \
            -nThreads {threads} {params.extra} -minMapQ {params.mapQ} \
            -minQ {params.baseQ} -sites {input.sites} -anc {input.anc} \
            -setMinDepthInd {params.mindepthind} {params.extra_saf} \
            -out {params.out}) &> {log}
        """
