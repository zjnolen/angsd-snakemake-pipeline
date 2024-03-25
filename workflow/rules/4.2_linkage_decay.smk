rule combine_LDdecay_files:
    input:
        ldgz=expand(
            "results/datasets/{{dataset}}/analyses/ngsLD/chunks/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}_{{sites}}-filts.ld_maxkbdist-{maxkb}_rndsample-{rndsmp}.gz",
            chunk=chunklist,
            maxkb=config["params"]["ngsld"]["max_kb_dist_decay"],
            rndsmp=config["params"]["ngsld"]["rnd_sample_decay"],
        ),
    output:
        ldgz=temp(
            "results/datasets/{dataset}/analyses/ngsLD/decay/{dataset}.{ref}_{population}{dp}_{sites}-filts.LDdecay.gz"
        ),
    log:
        "logs/{dataset}/ngsLD/combine_LDdecay_files/{dataset}.{ref}_{population}{dp}_{sites}-filts.ld_decay.log",
    benchmark:
        "benchmarks/{dataset}/ngsLD/combine_LDdecay_files/{dataset}.{ref}_{population}{dp}_{sites}-filts.ld_decay.log"
    conda:
        "../envs/shell.yaml"
    shell:
        """
        cat {input.ldgz} > {output.ldgz} 2> {log}
        """


rule fit_LD_decay:
    input:
        "results/datasets/{dataset}/analyses/ngsLD/decay/{dataset}.{ref}_{population}{dp}_{sites}-filts.LDdecay.gz",
    output:
        plot=report(
            "results/datasets/{dataset}/plots/LD_decay/{dataset}.{ref}_{population}{dp}_{sites}-filts.LDdecay.svg",
            category="01 Linkage Disequilibrium Decay",
            subcategory="{sites}",
            labels=lambda w: {
                "Population": "{population}",
                **dp_report(w),
                "Type": "Regression Plot",
            },
        ),
        fit="results/datasets/{dataset}/analyses/ngsLD/{dataset}.{ref}_{population}{dp}_{sites}-filts.LDdecay-fit.txt",
    log:
        "logs/{dataset}/ngsLD/fit_LD_decay/{dataset}.{ref}_{population}{dp}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/ngsLD/fit_LD_decay/{dataset}.{ref}_{population}{dp}_{sites}-filts.log"
    container:
        ngsld_container
    threads: lambda w, attempt: attempt
    resources:
        runtime=lambda w, attempt: attempt * 120,
    params:
        extra=config["params"]["ngsld"]["fit_LDdecay_extra"],
        nind=get_ngsld_n,
    shell:
        """
        echo {input} | fit_LDdecay.R {params.extra} {params.nind} --out {output.plot} \
                > {output.fit} 2> {log}
        """
