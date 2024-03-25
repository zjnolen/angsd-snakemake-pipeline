# Rules for visualizing and correcting signatures of post-mortem DNA damage


rule damageprofiler:
    """
    Estimates varous metrics related to post-mortem DNA damage. Informative 
    rather than corrective.
    """
    input:
        bam="results/mapping/bams/{sample}.{ref}.rmdup.realn.clip.bam",
        ref="results/ref/{ref}/{ref}.fa",
    output:
        multiext(
            "results/mapping/qc/damageprofiler/{sample}.{ref}/",
            "5pCtoT_freq.txt",
            "3pGtoA_freq.txt",
            "Length_plot.pdf",
            "DamagePlot_five_prime.svg",
            "DamagePlot.pdf",
            "DamagePlot_three_prime.svg",
            "DamageProfiler.log",
            "lgdistribution.txt",
            "edit_distance.svg",
            "edit_distance.pdf",
            "editDistance.txt",
            "Length_plot_combined_data.svg",
            "Length_plot_forward_reverse_separated.svg",
            "misincorporation.txt",
            "5p_freq_misincorporations.txt",
            "3p_freq_misincorporations.txt",
            "DNA_comp_genome.txt",
            "DNA_composition_sample.txt",
            "dmgprof.json",
        ),
    log:
        "logs/mapping/damageprofiler/{sample}.{ref}.log",
    benchmark:
        "benchmarks/mapping/damageprofiler/{sample}.{ref}.log"
    conda:
        "../envs/damageprofiler.yaml"
    params:
        out=lambda w, output: os.path.dirname(output[0]),
    threads: lambda wildcards, attempt: attempt
    resources:
        runtime=lambda wildcards, attempt: attempt * 60,
    shell:
        """
        damageprofiler -Xmx{resources.mem_mb}m -i {input.bam} -r {input.ref} \
            -o {params.out} &> {log}
        """


rule mapDamage2_rescaling:
    """
    Estimates various metrics related to post-mortem DNA damage and rescales base 
    quality scores to correct for damage.
    """
    input:
        bam="results/mapping/bams/{sample}.{ref}.rmdup.realn.clip.bam",
        ref="results/ref/{ref}/{ref}.fa",
    output:
        log="results/mapping/qc/mapdamage/{sample}.{ref}/Runtime_log.txt",
        GtoA3p="results/mapping/qc/mapdamage/{sample}.{ref}/3pGtoA_freq.txt",
        CtoT5p="results/mapping/qc/mapdamage/{sample}.{ref}/5pCtoT_freq.txt",
        dnacomp="results/mapping/qc/mapdamage/{sample}.{ref}/dnacomp.txt",
        frag_misincorp="results/mapping/qc/mapdamage/{sample}.{ref}/Fragmisincorporation_plot.pdf",
        len="results/mapping/qc/mapdamage/{sample}.{ref}/Length_plot.pdf",
        lg_dist="results/mapping/qc/mapdamage/{sample}.{ref}/lgdistribution.txt",
        misincorp="results/mapping/qc/mapdamage/{sample}.{ref}/misincorporation.txt",
        rescaled_bam="results/mapping/bams/{sample}.{ref}.rmdup.realn.clip.rescaled.bam",
    log:
        "logs/mapping/mapdamage/{sample}.{ref}.log",
    benchmark:
        "benchmarks/mapping/mapdamage/{sample}.{ref}.log"
    params:
        extra="--rescale",
    resources:
        runtime=1440,
        mem_mb=lambda w, attempt: attempt * 6400,
    wrapper:
        "v2.6.0/bio/mapdamage2"


rule dna_damage_multiqc:
    input:
        multiqc_input_dnadmg,
    output:
        report(
            "results/datasets/{dataset}/qc/dna-damage-mqc/dna-damage_all.{ref}_mqc.html",
            category="00 Quality Control",
            subcategory="4 DNA Damage",
            labels={"Type": "MultiQC Report"},
        ),
    log:
        "logs/mapping/dnadamage/{dataset}.{ref}_dnadmg-mqc.log",
    params:
        extra="--cl-config \"extra_fn_clean_exts: ['.rmdup']\" ",
    wrapper:
        "v3.5.0/bio/multiqc"
