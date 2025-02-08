

rule compress_homopolymers:
    input:
        illumina_dir=lambda wc: DATA_DIRS[wc.sm][f"illumina_{wc.hap}"],
    output:
        hpc_reads=join("results", "meryl", "{sm}", "{hap}_hpc.fastq.gz"),
    params:
        fq_glob=lambda wc: multi_flags(
            *dtype_glob(str(wc.sm), f"illumina_{wc.hap}"), opt="-name"
        ),
    log:
        join("logs", "meryl", "{sm}", "{hap}_compress_homopolymers.log"),
    benchmark:
        join("benchmarks", "meryl", "{sm}", "{hap}_compress_homopolymers.tsv")
    conda:
        "../envs/verkko.yaml"
    shell:
        """
        {{ seqtk hpc <(find {input.illumina_dir}/ {params.fq_glob} -size +0 -exec zcat -f {{}} +) | bgzip ;}}> {output} 2> {log}
        """


rule count_kmers:
    input:
        hpc_reads=rules.compress_homopolymers.output,
    output:
        mer_db=directory(join("results", "meryl", "{sm}", "{hap}_compress.meryl")),
    resources:
        mem=lambda wc: get_dtype_config(wc.sm, f"illumina_{wc.hap}").get("mem", "30GB"),
    threads: lambda wc: config["samples"][str(wc.sm)]["threads"] // 2
    log:
        join("logs", "meryl", "{sm}", "{hap}_count_kmers.log"),
    benchmark:
        join("benchmarks", "meryl", "{sm}", "{hap}_count_kmers.tsv")
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
        mat_kmers=lambda wc: expand(rules.count_kmers.output, sm=wc.sm, hap="mat"),
        pat_kmers=lambda wc: expand(rules.count_kmers.output, sm=wc.sm, hap="pat"),
    output:
        mat_db=directory(join("results", "meryl", "{sm}", "mat_compress.only.meryl")),
        pat_db=directory(join("results", "meryl", "{sm}", "pat_compress.only.meryl")),
    log:
        join("logs", "meryl", "{sm}", "generate_hapmers.log"),
    benchmark:
        join("benchmarks", "meryl", "{sm}", "generate_hapmers.tsv")
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


def phasing_data(wc):
    if "illumina_mat" in DATA_DIRS[wc.sm] and "illumina_pat" in DATA_DIRS[wc.sm]:
        return {"hap_kmers": expand(rules.generate_hapmers.output, sm=wc.sm)}
    else:
        return {}


def phasing_data_args(wc, inputs):
    if "hap_kmers" in inputs.keys():
        return f"--hap-kmers {inputs.hap_kmers} trio"
    else:
        return ""


rule run_verkko:
    input:
        # TODO: hic
        unpack(phasing_data),
        ont_dir=lambda wc: DATA_DIRS[wc.sm]["ont"],
        hifi_dir=lambda wc: DATA_DIRS[wc.sm]["hifi"],
    output:
        join("results", "verkko", "{sm}", "assembly.fasta"),
    conda:
        "../envs/verkko.yaml"
    threads: lambda wc: config["samples"][str(wc.sm)]["threads"]
    resources:
        mem=lambda wc: config["samples"][str(wc.sm)]["mem"] + "GB",
    log:
        "logs/run_verkko_{sm}.log",
    benchmark:
        "benchmarks/run_verkko_{sm}.tsv"
    params:
        output_dir=lambda wc, output: dirname(output[0]),
        ont_glob=lambda wc: multi_flags(*dtype_glob(str(wc.sm), "ont"), opt="-name"),
        hifi_glob=lambda wc: multi_flags(*dtype_glob(str(wc.sm), "hifi"), opt="-name"),
        phasing_data_args=lambda wc, input: phasing_data_args(wc, input),
        snakeopts=lambda wc: "--snakeopts "
        + config["samples"][str(wc.sm)].get("snakeopts", ""),
    shell:
        """
        verkko -d {params.output_dir} \
        --hifi $(find {input.hifi_dir}/ {params.hifi_glob}) \
        --nano $(find {input.ont_dir}/ {params.ont_glob}) \
        {params.phasing_data_args} \
        {params.snakeopts} &> {log}
        """


use rule generate_summary_stats as generate_summary_stats_verkko with:
    input:
        fa=rules.run_verkko.output,
    output:
        summary=join("results", "verkko", "{sm}", "assembly_stats.tsv"),
    log:
        "logs/generate_summary_stats_{sm}.log",


rule verkko_all:
    input:
        rules.run_verkko.output,
        rules.generate_summary_stats_verkko.output,
