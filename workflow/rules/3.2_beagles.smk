# Rules for generating beagle (SNP genotype likelihoods) and maf (minor allele
# frequency) files per population


rule angsd_doGlf2:
    """
    Generates beagle and minor allele frequency files for a given population and genome 
    chunk. Calls SNPs from the whole dataset, and uses these same sites across all 
    population beagle files, even if a population is fixed for a certain allele.
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
        beagle=temp(
            "results/datasets/{dataset}/beagles/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.beagle.gz"
        ),
        maf=temp(
            "results/datasets/{dataset}/beagles/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.mafs.gz"
        ),
        arg="results/datasets/{dataset}/beagles/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.arg",
    log:
        "logs/{dataset}/angsd/doGlf2/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/angsd/doGlf2/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.log"
    container:
        angsd_container
    params:
        gl_model=config["params"]["angsd"]["gl_model"],
        extra=config["params"]["angsd"]["extra"],
        extra_beagle=config["params"]["angsd"]["extra_beagle"],
        mapQ=config["mapQ"],
        baseQ=config["baseQ"],
        snp_pval=config["params"]["angsd"]["snp_pval"],
        maf=config["params"]["angsd"]["domaf"],
        minmaf=config["params"]["angsd"]["min_maf"],
        majmin=config["params"]["angsd"]["domajorminor"],
        counts=get_docounts,
        nind=lambda w: len(get_samples_from_pop(w.population)),
        out=lambda w, output: os.path.splitext(output.arg)[0],
    threads: lambda wildcards, attempt: attempt
    resources:
        runtime=lambda wildcards, attempt: attempt * 720,
    shell:
        """
        angsd -doGlf 2 -bam {input.bam} -GL {params.gl_model} -ref {input.ref} \
            -doMajorMinor {params.majmin} -doMaf {params.maf} -minMaf {params.minmaf} \
            -SNP_pval {params.snp_pval} -nThreads {threads} {params.extra} \
            -minMapQ {params.mapQ} -minQ {params.baseQ} -sites {input.sites} \
            -anc {input.anc} {params.extra_beagle} -rf {input.regions} \
            {params.counts} -out {params.out} &> {log}
        """


rule merge_beagle:
    """
    Merge beagle files across genome chunks per population.
    """
    input:
        lambda w: expand(
            "results/datasets/{{dataset}}/beagles/chunks/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}_{{sites}}-filts.beagle.gz",
            chunk=chunklist,
        ),
    output:
        beagle="results/datasets/{dataset}/beagles/{dataset}.{ref}_{population}{dp}_{sites}-filts.beagle.gz",
    log:
        "logs/{dataset}/angsd/doGlf2/{dataset}.{ref}_{population}{dp}_{sites}-filts_merge-beagle.log",
    benchmark:
        "benchmarks/{dataset}/angsd/doGlf2/{dataset}.{ref}_{population}{dp}_{sites}-filts_merge-beagle.log"
    conda:
        "../envs/shell.yaml"
    resources:
        runtime=lambda wildcards, attempt: attempt * 60,
    shell:
        r"""
        (set +o pipefail;
        zcat {input} | head -n 1 | gzip > {output}

        for f in {input}; do
            zcat $f | tail -n +2 | gzip | cat >> {output.beagle}
        done) 2> {log}
        """


rule merge_maf:
    """
    Merges maf files across genome chunks per population.
    """
    input:
        lambda w: expand(
            "results/datasets/{{dataset}}/beagles/chunks/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}_{{sites}}-filts.mafs.gz",
            chunk=chunklist,
        ),
    output:
        maf="results/datasets/{dataset}/mafs/{dataset}.{ref}_{population}{dp}_{sites}-filts.mafs.gz",
    log:
        "logs/{dataset}/angsd/doGlf2/{dataset}.{ref}_{population}{dp}_{sites}-filts_merge-mafs.log",
    benchmark:
        "benchmarks/{dataset}/angsd/doGlf2/{dataset}.{ref}_{population}{dp}_{sites}-filts_merge-mafs.log"
    conda:
        "../envs/shell.yaml"
    resources:
        runtime=lambda wildcards, attempt: attempt * 60,
    shell:
        r"""
        (set +o pipefail;
        zcat {input} | head -n 1 | gzip > {output.maf}

        for f in {input}; do
            zcat $f | tail -n +2 | gzip | cat >> {output.maf}
        done) 2> {log}
        """


rule snpset:
    """
    Extracts SNP coordinates from maf file to get a list of variable sites in the 
    dataset.
    """
    input:
        "results/datasets/{dataset}/mafs/{dataset}.{ref}_all{dp}_{sites}-filts.mafs.gz",
    output:
        "results/datasets/{dataset}/filters/snps/{dataset}.{ref}{dp}_{sites}-filts_snps.sites",
    log:
        "logs/{dataset}/filters/snps/{dataset}.{ref}{dp}_{sites}-filts_snps.log",
    benchmark:
        "benchmarks/{dataset}/filters/snps/{dataset}.{ref}{dp}_{sites}-filts_snps.log"
    conda:
        "../envs/shell.yaml"
    shell:
        """
        zcat {input} | tail -n +2 | cut -f1-4 > {output} 2> {log}
        """
