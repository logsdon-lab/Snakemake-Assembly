rule compress_homopolymers:
    input:
        illumina_dir=lambda wc: DATA_DIRS[wc.sm][f"illumina_{wc.hap}"],
    output:
        hpc_reads=join(OUTPUT_DIR, "meryl", "{sm}", "{hap}_hpc.fastq.gz"),
    params:
        fq_glob=lambda wc: multi_flags(
            *dtype_glob(str(wc.sm), f"illumina_{wc.hap}"), pre_opt="-name"
        ),
    log:
        join(LOG_DIR, "meryl", "{sm}", "{hap}_compress_homopolymers.log"),
    benchmark:
        join(BENCHMARK_DIR, "meryl", "{sm}", "{hap}_compress_homopolymers.tsv")
    conda:
        "../envs/verkko.yaml"
    shell:
        """
        {{ seqtk hpc <(find {input.illumina_dir}/ {params.fq_glob} -size +0 -exec zcat -f {{}} +) | bgzip ;}}> {output} 2> {log}
        """


rule count_kmers_meryl:
    input:
        hpc_reads=rules.compress_homopolymers.output,
    output:
        mer_db=directory(join(OUTPUT_DIR, "meryl", "{sm}", "{hap}_compress.meryl")),
    resources:
        mem=lambda wc: get_dtype_config(wc.sm, f"illumina_{wc.hap}").get("mem", "30GB"),
    threads: lambda wc: config["samples"][str(wc.sm)]["threads"] // 2
    log:
        join(LOG_DIR, "meryl", "{sm}", "{hap}_count_kmers.log"),
    benchmark:
        join(BENCHMARK_DIR, "meryl", "{sm}", "{hap}_count_kmers.tsv")
    params:
        kmer_size=lambda wc: get_dtype_config(wc.sm, f"illumina_{wc.hap}").get(
            "kmer_size", 30
        ),
        mem=lambda wc, resources: str(resources.mem).replace("GB", ""),
    conda:
        "../envs/verkko.yaml"
    shell:
        """
        meryl count compress \
        k={params.kmer_size} \
        threads={threads} \
        memory={params.mem} \
        {input.hpc_reads} \
        output {output.mer_db} 2> {log}
        """


rule generate_hapmers:
    input:
        mat_kmers=lambda wc: expand(rules.count_kmers_meryl.output, sm=wc.sm, hap="mat"),
        pat_kmers=lambda wc: expand(rules.count_kmers_meryl.output, sm=wc.sm, hap="pat"),
    output:
        mat_db=directory(join(OUTPUT_DIR, "meryl", "{sm}", "mat_compress.only.meryl")),
        pat_db=directory(join(OUTPUT_DIR, "meryl", "{sm}", "pat_compress.only.meryl")),
    log:
        join(LOG_DIR, "meryl", "{sm}", "generate_hapmers.log"),
    benchmark:
        join(BENCHMARK_DIR, "meryl", "{sm}", "generate_hapmers.tsv")
    params:
        output_dir=lambda wc, output: dirname(output.mat_db),
    conda:
        "../envs/verkko.yaml"
    shell:
        """
        log=$(realpath {log})
        cd {params.output_dir}
        $MERQURY/trio/hapmers.sh \
            {input.mat_kmers} \
            {input.pat_kmers} &> ${{log}}
        """
