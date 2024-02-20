# Creates a filtered list of sites for analyses to restrict the analyses to. This
# filtering limits analyses to autosomal scaffolds and removes sites with low
# mappability, low complexity, and excessively high or low depth. This filtering
# regime was adapted from Pečnerová et al. 2021 (Current Biology).


localrules:
    genome_bed,
    smallscaff_bed,
    sexlink_bed,
    genmap_filt_bed,
    repeat_sum,


rule genome_bed:
    """Create a bed file containing the entirety of the genome"""
    input:
        fai="results/ref/{ref}/{ref}.fa.fai",
    output:
        bed="results/ref/{ref}/beds/genome.bed",
        sum="results/ref/{ref}/beds/genome.bed.sum",
    log:
        "logs/ref/genome_bed/{ref}.log",
    benchmark:
        "benchmarks/ref/genome_bed/{ref}.log"
    conda:
        "../envs/shell.yaml"
    shell:
        r"""
        (# generate bed
        awk -v OFS='\t' '{{print $1, "0", $2}}' {input.fai} > {output.bed}

        # summarize bed
        len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' {output.bed})
        echo $len | awk '{{print "Total genome\t"$1"\t"$1/$1*100}}' > \
            {output.sum}) 2> {log}
        """


rule smallscaff_bed:
    """Create a bed file of all scaffolds under a specified size"""
    input:
        genbed="results/ref/{ref}/beds/genome.bed",
        gensum="results/ref/{ref}/beds/genome.bed.sum",
    output:
        bed="results/datasets/{dataset}/filters/small_scaffs/{ref}_scaff{size}bp.bed",
        sum="results/datasets/{dataset}/filters/small_scaffs/{ref}_scaff{size}bp.bed.sum",
    log:
        "logs/{dataset}/filters/smallscaff/{ref}_scaff{size}bp.log",
    benchmark:
        "benchmarks/{dataset}/filters/smallscaff/{ref}_scaff{size}bp.log"
    conda:
        "../envs/shell.yaml"
    params:
        minsize="{size}",
    shell:
        r"""
        (# generate bed
        awk '$3 < {params.minsize}' {input.genbed} > {output.bed}

        # summarize bed
        len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' {output.bed})
        echo $len $(awk -F "\t" '{{print $2}}' {input.gensum}) | \
            awk '{{print "Scaffolds<{params.minsize}bp\t"$2-$1"\t"($2-$1)/$2*100}}' \
            > {output.sum}) 2> {log}
        """


rule sexlink_bed:
    """Create bed files of specified sex-linked contigs and other excluded contigs"""
    input:
        genbed="results/ref/{ref}/beds/genome.bed",
        gensum="results/ref/{ref}/beds/genome.bed.sum",
    output:
        exclbed="results/datasets/{dataset}/filters/sex-link_mito_excl/{ref}_excl.bed",
        sum="results/datasets/{dataset}/filters/sex-link_mito_excl/{ref}_excl.bed.sum",
        sexbed="results/datasets/{dataset}/filters/sex-link_mito_excl/{ref}_sex-linked.bed",
    log:
        "logs/{dataset}/filters/sex-link_mito_excl/{ref}.log",
    benchmark:
        "benchmarks/{dataset}/filters/sex-link_mito_excl/{ref}.log"
    conda:
        "../envs/shell.yaml"
    params:
        sex=config["reference"]["sex-linked"],
        excl=config["reference"]["exclude"],
        mito=config["reference"]["mito"],
    shell:
        r"""
        (# generate beds
        printf '%s\n' {params.sex} | grep -f - {input.genbed} > {output.sexbed}
        printf '%s\n' {params.sex} {params.excl} {params.mito} | \
            grep -f - {input.genbed} > {output.exclbed}
        
        # summarize exclbed
        len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' {output.exclbed})
        echo $len $(awk -F "\t" '{{print $2}}' {input.gensum}) | \
            awk '{{print "Autosomes\t"$2-$1"\t"($2-$1)/$2*100}}' > {output.sum}
        ) 2> {log}
        """


rule genmap_index:
    """Index reference for Genmap"""
    input:
        ref="results/ref/{ref}/{ref}.fa",
    output:
        fold=directory("results/ref/{ref}/genmap/index"),
        files=multiext(
            "results/ref/{ref}/genmap/index/index",
            ".ids.concat",
            ".ids.limits",
            ".info.concat",
            ".info.limits",
            ".lf.drp",
            ".lf.drp.sbl",
            ".lf.drs",
            ".lf.drv",
            ".lf.drv.sbl",
            ".lf.pst",
            ".rev.lf.drp",
            ".rev.lf.drp.sbl",
            ".rev.lf.drs",
            ".rev.lf.drv",
            ".rev.lf.drv.sbl",
            ".rev.lf.pst",
            ".sa.ind",
            ".sa.len",
            ".sa.val",
            ".txt.concat",
            ".txt.limits",
        ),
    log:
        "logs/ref/genmap/index/{ref}.log",
    benchmark:
        "benchmarks/ref/genmap/index/{ref}.log"
    conda:
        "../envs/genmap.yaml"
    threads: lambda wildcards, attempt: attempt
    shell:
        """
        # genmap index annoyingly fails if directory already exists,
        # delete it to keep it happy
        rm -r {output.fold} 2> {log}

        genmap index -F {input.ref} -I {output.fold} &>> {log}
        """


rule genmap_map:
    """Estimate mappability of each site in the genome"""
    input:
        fold=rules.genmap_index.output.fold,
        files=rules.genmap_index.output.files,
    output:
        bed="results/ref/{ref}/genmap/map/{ref}_k{k}_e{e}.bedgraph",
    log:
        "logs/ref/genmap/map/{ref}_k{k}_e{e}.log",
    benchmark:
        "benchmarks/ref/genmap/map/{ref}_k{k}_e{e}.log"
    conda:
        "../envs/genmap.yaml"
    params:
        out=lambda w, output: os.path.splitext(output.bed)[0],
    threads: lambda wildcards, attempt: attempt
    resources:
        runtime=lambda wildcards, attempt: attempt * 360,
    shell:
        """
        genmap map -K {wildcards.k} -E {wildcards.e} -I {input.fold} \
            -O {params.out} -bg &> {log}
        """


rule windowgen:
    """
    Generate sliding windows across the genome to average mappability over.
    """
    input:
        bed="results/ref/{ref}/beds/genome.bed",
    output:
        bed="results/ref/{ref}/genmap/pileup/{ref}_genome_windows_k{k}.bed",
    log:
        "logs/ref/genmap/windowgen/{ref}_k{k}.log",
    benchmark:
        "benchmarks/ref/genmap/windowgen/{ref}_k{k}.log"
    conda:
        "../envs/bedtools.yaml"
    resources:
        runtime=720,
    shell:
        """
        bedtools makewindows -b {input.bed} -w {wildcards.k} -s 1 > {output.bed} 2> {log}
        """


rule pileup_mappability:
    """
    Generates a bed file containing the pileup mappability of a site, i.e. the
    mean mappability of all possible kmers mapping to it.
    """
    input:
        bed="results/ref/{ref}/genmap/pileup/{ref}_genome_windows_k{k}.bed",
        bgr="results/ref/{ref}/genmap/map/{ref}_k{k}_e{e}.bedgraph",
    output:
        bed="results/ref/{ref}/genmap/pileup/{ref}_pileup_mappability_k{k}_e{e}.bed",
        tmp=temp("results/ref/{ref}/genmap/map/{ref}_k{k}_e{e}.bedgraph.tmp"),
    log:
        "logs/ref/genmap/pileupmap/{ref}_k{k}_e{e}.log",
    benchmark:
        "benchmarks/ref/genmap/pileupmap/{ref}_k{k}_e{e}.log"
    conda:
        "../envs/bedops.yaml"
    resources:
        runtime=720,
    shell:
        r"""
        awk '{{print $1"\t"$2"\t"$3"\t"$2"-"$3"\t"$4}}' {input.bgr} > {output.tmp}
        bedmap --echo --wmean {input.bed} {output.tmp} | tr "|" "\t" | \
            awk '{{print $1"\t"$3-1"\t"$3"\t"$4}}' > {output.bed} 2> {log}
        """


rule genmap_filt_bed:
    """Create a bed containing all sites with a mappability below a set threshold"""
    input:
        genbed="results/ref/{ref}/genmap/pileup/{ref}_pileup_mappability_k{k}_e{e}.bed",
        gensum="results/ref/{ref}/beds/genome.bed.sum",
    output:
        bed="results/datasets/{dataset}/filters/pileupmap/{ref}_k{k}_e{e}_{thresh}.bed",
        tmp="results/datasets/{dataset}/filters/pileupmap/{ref}_k{k}_e{e}_{thresh}.bed.tmp",
        sum="results/datasets/{dataset}/filters/pileupmap/{ref}_k{k}_e{e}_{thresh}.bed.sum",
    log:
        "logs/{dataset}/filters/pileupmap/{ref}_k{k}_e{e}_{thresh}.log",
    benchmark:
        "benchmarks/{dataset}/filters/pileupmap/{ref}_k{k}_e{e}_{thresh}.log"
    conda:
        "../envs/bedtools.yaml"
    params:
        thresh=config["params"]["genmap"]["map_thresh"],
    shell:
        r"""
        # generate bed
        (awk '$4 < {params.thresh}' {input.genbed} > {output.tmp}
        bedtools merge -i {output.tmp} > {output.bed}

        #summarize bed
        len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' {output.bed})
        echo $len $(awk -F "\t" '{{print $2}}' {input.gensum}) | awk \
            '{{print "Pileup mappability K{wildcards.k}-E{wildcards.e} <{params.thresh}\t"$2-$1"\t"($2-$1)/$2*100}}' \
            > {output.sum}) 2> {log}
        """


rule repeat_builddatabase:
    """Build database for RepeatModeler"""
    input:
        ref="results/ref/{ref}/{ref}.fa",
    output:
        multiext(
            "results/ref/{ref}/repeatmodeler/{ref}.",
            "nhr",
            "nin",
            "nnd",
            "nni",
            "nog",
            "nsq",
            "translation",
        ),
    conda:
        "../envs/repeatmasker.yaml"
    log:
        "logs/ref/repeatmodeler/builddatabase/{ref}.log",
    benchmark:
        "benchmarks/ref/repeatmodeler/builddatabase/{ref}.log"
    params:
        db=lambda w, output: os.path.splitext(output[0])[0],
    shell:
        """
        BuildDatabase -name {params.db} {input.ref} &> {log}
        """


rule repeatmodeler:
    """Model repeats across the reference"""
    input:
        database=rules.repeat_builddatabase.output,
    output:
        fa="results/ref/{ref}/repeatmodeler/{ref}-families.fa",
        stk="results/ref/{ref}/repeatmodeler/{ref}-families.stk",
        log="results/ref/{ref}/repeatmodeler/{ref}-rmod.log",
    log:
        "logs/ref/repeatmodeler/repeatmodeler/{ref}.log",
    benchmark:
        "benchmarks/ref/repeatmodeler/repeatmodeler/{ref}.log"
    conda:
        "../envs/repeatmasker.yaml"
    params:
        db=lambda w, input: os.path.splitext(input[0])[0],
        ref="{ref}",
    threads: 10
    resources:
        runtime=10080,
    shadow:
        "minimal"
    shell:
        """
        RepeatModeler -database {params.db} -pa {threads} &> {log}
        """


rule repeatmasker:
    """Create gff file containing repeats within the genome"""
    input:
        unpack(get_repmaskin),
    output:
        gff="results/ref/{ref}/repeatmasker/{ref}.fa.out.gff",
    log:
        "logs/ref/repeatmasker/{ref}.log",
    benchmark:
        "benchmarks/ref/repeatmasker/{ref}.log"
    conda:
        "../envs/repeatmasker.yaml"
    params:
        out=lambda w, output: os.path.dirname(output.gff),
        libpre="-species" if config["analyses"]["repeatmasker"]["dfam_lib"] else "-lib",
        lib=lambda w, input: f"'{config['analyses']['repeatmasker']['dfam_lib']}'"
        if config["analyses"]["repeatmasker"]["dfam_lib"]
        else input.lib,
    threads: 5
    resources:
        runtime=720,
    shadow:
        "shallow"
    shell:
        """
        RepeatMasker -pa {threads} {params.libpre} {params.lib} -gff -x -no_is \
            -dir {params.out} {input.ref} &> {log}
        """


rule repeat_sum:
    """Summarize the proportion of the genome contained in the repeat gff"""
    input:
        unpack(get_rep_file),
        sum="results/ref/{ref}/beds/genome.bed.sum",
    output:
        sum="results/ref/{ref}/repeatmasker/{ref}.fa.out.sum",
        bed="results/ref/{ref}/repeatmasker/{ref}.fa.out.bed",
    log:
        "logs/ref/repeatmasker/summarize_gff/{ref}.log",
    benchmark:
        "benchmarks/ref/repeatmasker/summarize_gff/{ref}.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        r"""
        (bedtools merge -i {input.rep} > {output.bed}
        len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' {output.bed})
        echo $len $(awk -F "\t" '{{print $2}}' {input.sum}) | \
            awk '{{print "Repeats\t"$2-$1"\t"($2-$1)/$2*100}}' > {output.sum}) &> {log}
        """


if config["analyses"]["extreme_depth"]:

    rule angsd_depth:
        """Estimate global depth for different subsets of samples, performed in chunks"""
        input:
            bamlist="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
            regions="results/datasets/{dataset}/filters/chunks/{ref}_chunk{chunk}.rf",
            ref="results/ref/{ref}/{ref}.fa",
            bams=get_bamlist_bams,
            bais=get_bamlist_bais,
        output:
            posgz=temp(
                "results/datasets/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}_chunk{chunk}.pos.gz"
            ),
            hist=temp(
                "results/datasets/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}_chunk{chunk}.depthGlobal"
            ),
            samphist=temp(
                "results/datasets/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}_chunk{chunk}.depthSample"
            ),
            arg=temp(
                "results/datasets/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}_chunk{chunk}.arg"
            ),
        log:
            "logs/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}_chunk{chunk}.log",
        benchmark:
            "benchmarks/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}_chunk{chunk}.log"
        container:
            angsd_container
        params:
            extra=config["params"]["angsd"]["extra"],
            mapQ=config["mapQ"],
            baseQ=config["baseQ"],
            out=lambda w, output: os.path.splitext(output.arg)[0],
        threads: lambda wildcards, attempt: attempt * 2
        resources:
            runtime=lambda wildcards, attempt: attempt * 720,
        shell:
            """
            (nInd=$(cat {input.bamlist} | wc -l | awk '{{print $1+1}}')
            maxDP=$(echo 1000 $nInd | awk '{{print $1 * $2}}')

            angsd -bam {input.bamlist} -nThreads {threads} -rf {input.regions} \
                -ref {input.ref} -minMapQ {params.mapQ} -minQ {params.baseQ} -doCounts 1 \
                -dumpCounts 1 {params.extra} -doDepth 1 -maxDepth $maxDP \
                -out {params.out}) 2> {log}
            """

    rule combine_depths:
        """Merge global depth files across chunks"""
        input:
            lambda w: expand(
                "results/datasets/{{dataset}}/filters/depth/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}.depthGlobal",
                chunk=chunklist,
            ),
        output:
            "results/datasets/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}.depthGlobal",
        log:
            "logs/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}_combined.log",
        benchmark:
            "benchmarks/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}_combined.log"
        conda:
            "../envs/shell.yaml"
        shell:
            """
            cat {input} > {output} 2> {log}
            """

    rule summarize_depths:
        """Estimate mean and bounds of middle 95% of the global depth distribution"""
        input:
            "results/datasets/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}.depthGlobal",
        output:
            summ="results/datasets/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}_depth.summary",
            plot=report(
                "results/datasets/{dataset}/plots/depth_dist/{dataset}.{ref}_{population}{dp}_depth.svg",
                category="Quality Control",
                subcategory="Depth distributions and filters",
                labels={"Subset": "{population}", "Type": "Histogram"},
            ),
        log:
            "logs/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}_depth_extremes.log",
        benchmark:
            "benchmarks/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}_depth_extremes.log"
        conda:
            "../envs/r.yaml"
        params:
            lower=config["params"]["extreme_depth_filt"]["bounds"][0],
            upper=config["params"]["extreme_depth_filt"]["bounds"][1],
            method=config["params"]["extreme_depth_filt"]["method"],
        threads: lambda wildcards, attempt: attempt * 2
        script:
            "../scripts/depth_extremes.R"

    rule depth_bed:
        """Create a bed file containing regions of extreme depth in a subset"""
        input:
            genbed="results/ref/{ref}/beds/genome.bed",
            gensum="results/ref/{ref}/beds/genome.bed.sum",
            quants="results/datasets/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}_depth.summary",
            pos=lambda w: expand(
                "results/datasets/{{dataset}}/filters/depth/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}.pos.gz",
                chunk=chunklist,
            ),
        output:
            bed="results/datasets/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}_extreme-depth.bed",
            sum="results/datasets/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}_extreme-depth.bed.sum",
        log:
            "logs/{dataset}/filters/depth/bed/{dataset}.{ref}_{population}{dp}.log",
        benchmark:
            "benchmarks/{dataset}/filters/depth/bed/{dataset}.{ref}_{population}{dp}.log"
        conda:
            "../envs/bedtools.yaml"
        threads: lambda wildcards, attempt: attempt
        shell:
            r"""
            (lower=$(awk '{{print $2}}' {input.quants})
            upper=$(awk '{{print $3}}' {input.quants})
            for i in {input.pos}; do
                zcat $i | tail -n +2 | \
                awk -v lower=$lower -v upper=$upper '$3 > lower && $3 < upper'
            done | \
            awk '{{print $1"\t"$2-1"\t"$2}}' > {output.bed}
            bedtools merge -i {output.bed} > {output.bed}.tmp
            bedtools subtract -a {input.genbed} -b {output.bed}.tmp > {output.bed}
            rm {output.bed}.tmp
            
            # summarize bed
            len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' {output.bed})
            echo $len $(awk -F "\t" '{{print $2}}' {input.gensum}) | \
                awk '{{print "Depth ({wildcards.population})\t"$2-$1"\t"($2-$1)/$2*100}}' \
                > {output.sum}) 2> {log}
            """


# rule combine_depth_bed:
#     """
#     Combine beds for each subset to get a set of regions with extreme depth in any
#     subset
#     """
#     input:
#         beds=expand(
#             "results/datasets/{{dataset}}/filters/depth/{{dataset}}.{{ref}}_{population}{{dp}}_extreme-depth.bed",
#             population=["all"] + list(set(samples.depth.values)),
#         ),
#         gensum="results/ref/{ref}/beds/genome.bed.sum",
#     output:
#         bed="results/datasets/{dataset}/filters/depth/{dataset}.{ref}{dp}_extreme-depth.bed",
#         sum="results/datasets/{dataset}/filters/depth/{dataset}.{ref}{dp}_extreme-depth.bed.sum",
#     log:
#         "logs/{dataset}/filters/depth/bed/{dataset}.{ref}{dp}_combine-bed.log",
#     benchmark:
#         "benchmarks/{dataset}/filters/depth/bed/{dataset}.{ref}{dp}_combine-bed.log"
#     conda:
#         "../envs/bedtools.yaml"
#     shadow:
#         "minimal"
#     resources:
#         runtime=240,
#     shell:
#         """
#         # combine beds
#         (cat {input.beds} > {output.bed}.tmp
#         sort -k1,1 -k2,2n {output.bed}.tmp > {output.bed}.tmp.sort
#         rm {output.bed}.tmp

#         bedtools merge -i {output.bed}.tmp.sort > {output.bed}
#         rm {output.bed}.tmp.sort

#         # summarize bed
#         len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' {output.bed})
#         echo $len $(awk -F "\t" '{{print $2}}' {input.gensum}) | \
#             awk '{{print "Depth\t"$2-$1"\t"($2-$1)/$2*100}}' \
#             > {output.sum}) 2> {log}
#         """


rule angsd_missdata:
    """
    Print sites with data for more than a set proportion of individuals per population
    and across the whole dataset
    """
    input:
        bamlist="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
        regions="results/datasets/{dataset}/filters/chunks/{ref}_chunk{chunk}.rf",
        ref="results/ref/{ref}/{ref}.fa",
        bams=get_bamlist_bams,
        bais=get_bamlist_bais,
    output:
        posgz=temp(
            "results/datasets/{dataset}/filters/missdata/{dataset}.{ref}_{population}{dp}_chunk{chunk}_over{miss}.pos.gz"
        ),
        arg=temp(
            "results/datasets/{dataset}/filters/missdata/{dataset}.{ref}_{population}{dp}_chunk{chunk}_over{miss}.arg"
        ),
    log:
        "logs/{dataset}/filters/angsd_missdata/{dataset}.{ref}_{population}{dp}_chunk{chunk}_over{miss}.log",
    benchmark:
        "benchmarks/{dataset}/filters/angsd_missdata/{dataset}.{ref}_{population}{dp}_chunk{chunk}_over{miss}.log"
    container:
        angsd_container
    params:
        nind=lambda w: len(get_samples_from_pop(w.population)),
        extra=config["params"]["angsd"]["extra"],
        mapQ=config["mapQ"],
        baseQ=config["baseQ"],
        out=lambda w, output: os.path.splitext(output.arg)[0],
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        runtime=lambda wildcards, attempt: attempt * 360,
    shell:
        """
        (minInd=$(echo {params.nind} \
            | awk '{{print $1*{wildcards.miss}}}' \
            | awk '{{print int($1) + ( $1!=int($1) && $1>=0 )}}')
        
        angsd -bam {input.bamlist} -nThreads {threads} -rf {input.regions} \
            -ref {input.ref} -doCounts 1 -dumpCounts 1 -minInd $minInd \
            {params.extra} -minMapQ {params.mapQ} -minQ {params.baseQ} \
            -out {params.out}) 2> {log}
        """


rule missdata_bed:
    """Create bed file containing only sites passing all missing data thresholds"""
    input:
        pos=lambda w: expand(
            "results/datasets/{{dataset}}/filters/missdata/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}_over{{miss}}.pos.gz",
            chunk=chunklist,
        ),
        genbed="results/ref/{ref}/beds/genome.bed",
        gensum="results/ref/{ref}/beds/genome.bed.sum",
    output:
        bed="results/datasets/{dataset}/filters/missdata/{dataset}.{ref}_{population}{dp}_under{miss}.bed",
        sum="results/datasets/{dataset}/filters/missdata/{dataset}.{ref}_{population}{dp}_under{miss}.bed.sum",
        tmp=temp(
            "results/datasets/{dataset}/filters/missdata/{dataset}.{ref}_{population}{dp}_under{miss}.bed.tmp"
        ),
    log:
        "logs/{dataset}/filters/missdata_bed/{dataset}.{ref}_{population}{dp}_under{miss}.log",
    benchmark:
        "benchmarks/{dataset}/filters/missdata_bed/{dataset}.{ref}_{population}{dp}_under{miss}.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        r"""
        # generate bed
        (> {output.tmp}
        for i in {input.pos}; do
            zcat $i | tail -n +2 >> {output.tmp}
        done
        
        awk '{{print $1"\t"$2-1"\t"$2}}' {output.tmp} > {output.bed}
        bedtools merge -i {output.bed} > {output.tmp}
        bedtools subtract -a {input.genbed} -b {output.tmp} > {output.bed}
        
        # summarize bed
        perc=$(echo {wildcards.miss} | awk '{{print $1*100}}')
        len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' {output.bed})
        echo $len $(awk -F "\t" '{{print $2}}' {input.gensum}) $perc | \
            awk '{{print "Min "$3"% data ({wildcards.population})\t"$2-$1"\t" \
            ($2-$1)/$2*100}}' > {output.sum}) 2> {log}
        """


rule combine_beds:
    """
    Subtract all the BED files produced above from the whole genome BED to get a list
    of filtered sites to use for analyses
    """
    input:
        unpack(get_bed_filts),
    output:
        bed="results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_allsites-filts.bed",
        lis="results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_allsites-filts.list",
        sit="results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_allsites-filts.sites",
        sum="results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_allsites-filts.sum",
    log:
        "logs/{dataset}/filters/combine/{dataset}.{ref}{dp}_combine_beds.log",
    benchmark:
        "benchmarks/{dataset}/filters/combine/{dataset}.{ref}{dp}_combine_beds.log"
    conda:
        "../envs/bedtools.yaml"
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        runtime=240,
    shell:
        r"""
        (printf '%s\n' {input.filt} > {output.lis}
        cat {input.gen} > {output.bed}

        printf "Name\tLength(bp)\tPercent\n" > {output.sum}
        cat {input.sum} >> {output.sum}

        for i in {input.filt}; do
            bedtools subtract -a {output.bed} -b $i > {output.bed}.tmp
            mv {output.bed}.tmp {output.bed}
        done

        for i in {input.sums}; do
            cat $i >> {output.sum}
        done

        filtlen=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' \
            {output.bed})
        echo $filtlen $(awk -F "\t" '{{print $2}}' {input.sum}) | \
            awk '{{print "Combined\t"$1"\t"$1/$2*100}}' >> {output.sum}

        awk '{{print $1"\t"$2+1"\t"$3}}' {output.bed} > {output.sit}.tmp
        sort -V {output.sit}.tmp > {output.sit}
        rm {output.sit}.tmp) 2> {log}
        """


rule user_sites:
    """
    When users provide subsets of the genome to provide analyses on, create a bed and
    sites file for each subset, to limit analyses with.
    """
    input:
        gen="results/ref/{ref}/beds/genome.bed",
        gensum="results/ref/{ref}/beds/genome.bed.sum",
        newfilt=get_newfilt,
        allfilt="results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_allsites-filts.bed",
        sum="results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_allsites-filts.sum",
    output:
        bed="results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_{sites}-filts.bed",
        tmp=temp(
            "results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_{sites}-filts.bed.tmp"
        ),
        sit="results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_{sites}-filts.sites",
        sum="results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_{sites}-filts.sum",
    wildcard_constraints:
        sites="|".join(list(config["filter_beds"].keys())),
    log:
        "logs/{dataset}/filters/user_sites/{dataset}.{ref}{dp}_{sites}-filt.log",
    benchmark:
        "benchmarks/{dataset}/filters/user_sites/{dataset}.{ref}{dp}_{sites}-filt.log"
    conda:
        "../envs/bedtools.yaml"
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        runtime=240,
    shell:
        r"""
        (sed \$d {input.sum} > {output.sum}

        siteslen=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' {input.newfilt})
        echo $siteslen $(awk -F "\t" '{{print $2}}' {input.gensum}) | \
            awk '{{print "{wildcards.sites}-filts\t"$1"\t"$1/$2*100}}' \
            >> {output.sum}
        bedtools subtract -a {input.gen} -b {input.newfilt} > {output.tmp}
        bedtools subtract -a {input.allfilt} -b {output.tmp} > {output.bed}

        filtlen=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' {output.bed})
        echo $filtlen $(awk -F "\t" '{{print $2}}' {input.gensum}) | \
            awk '{{print "Combined\t"$1"\t"$1/$2*100}}' >> {output.sum}
        
        awk '{{print $1"\t"$2+1"\t"$3}}' {output.bed} | sort -V > {output.sit}) 2> {log}
        """


rule nofilts_sites_file:
    input:
        gen="results/ref/{ref}/beds/genome.bed",
    output:
        sites="results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_no-filts.sites"
    conda:
        "../envs/shell.yaml"
    shell:
        """
        awk '{{print $1"\t"$2+1"\t"$3}}' {input.gen} > {output.sites}
        """


rule filter_summary_table:
    """Produce table from filter summary to incorporate into report"""
    input:
        "results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_{sites}-filts.sum",
    output:
        report(
            "results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_{sites}-filts.html",
            category="Quality Control",
            subcategory="Filtering Summary",
            labels={"Filter": "{sites}", "Type": "Table"},
        ),
    log:
        "logs/{dataset}/filters/combine/{dataset}.{ref}{dp}_{sites}-filts_tsv2html.log",
    benchmark:
        "benchmarks/{dataset}/filters/combine/{dataset}.{ref}{dp}_{sites}-filts_tsv2html.log"
    conda:
        "../envs/r-rectable.yaml"
    script:
        "../scripts/tsv2html.Rmd"
