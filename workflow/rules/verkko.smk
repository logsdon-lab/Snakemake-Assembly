

rule compress_homopolymers:
    input:
        illumina_dir=lambda wc: DATA_DIRS[wc.sm][f"illumina_{wc.hap}"],
    output:
        hpc_reads=join("results", "meryl", "{sm}", "{hap}_hpc.fastq.gz"),
    params:
        fq_glob=lambda wc: multi_flags(
            *dtype_glob(str(wc.sm), f"illumina_{wc.hap}"), pre_opt="-name"
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


rule count_kmers_meryl:
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
        mat_kmers=lambda wc: expand(rules.count_kmers_meryl.output, sm=wc.sm, hap="mat"),
        pat_kmers=lambda wc: expand(rules.count_kmers_meryl.output, sm=wc.sm, hap="pat"),
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


# TODO: Make assembly a wildcard.
def phasing_data_verkko(wc) -> dict[str, Any]:
    if "illumina_mat" in DATA_DIRS[wc.sm] and "illumina_pat" in DATA_DIRS[wc.sm]:
        return {"hap_kmers": expand(rules.generate_hapmers.output, sm=wc.sm)}
    elif "hic_mat" in DATA_DIRS[wc.sm] and "hic_pat" in DATA_DIRS[wc.sm]:
        return {
            "mat_hic_dir": DATA_DIRS[wc.sm]["hic_mat"],
            "pat_hic_dir": DATA_DIRS[wc.sm]["hic_pat"],
        }
    else:
        return {}


def phasing_data_verkko_args(wc, inputs) -> str:
    if "hap_kmers" in inputs.keys():
        return f"--hap-kmers {inputs.hap_kmers} trio"
    elif "mat_hic_dir" in inputs.keys() and "pat_hic_dir" in inputs.keys():
        # Construct find command.
        mat_find_cmd = find_sep_command(
            wc.sm, "hic_mat", inputs.mat_hic_dir[0], sep=" "
        )
        pat_find_cmd = find_sep_command(
            wc.sm, "hic_pat", inputs.pat_hic_dir[0], sep=" "
        )
        return f"--hic1 $({mat_find_cmd}) --hic2 $({pat_find_cmd})"
    else:
        return ""


def input_verkko_args(wc, inputs) -> str:
    args = ["--hifi", f"$(cat {inputs.hifi_fofn})"]
    if inputs["ont_fofn"]:
        args.append("--nano")
        args.append(f"$(cat {inputs.ont_fofn})")
    return " ".join(args)


rule run_verkko:
    input:
        unpack(phasing_data_verkko),
        ont_fofn=lambda wc: (
            expand(
                rules.generate_dtype_fofn.output, asm="verkko", dtype="ont", sm=wc.sm
            )
            if config["samples"][str(wc.sm)]["data"].get("ont")
            else []
        ),
        hifi_fofn=expand(
            rules.generate_dtype_fofn.output, asm="verkko", dtype="hifi", sm="{sm}"
        ),
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
        phasing_data_args=lambda wc, input: phasing_data_verkko_args(wc, input),
        input_args=lambda wc, input: input_verkko_args(wc, input),
        snakeopts=lambda wc: (
            f'--snakeopts {config["samples"][str(wc.sm)]["snakeopts"]}'
            if config["samples"][str(wc.sm)].get("snakeopts")
            else ""
        ),
    shell:
        """
        verkko -d {params.output_dir} \
        {params.input_args} \
        {params.phasing_data_args} \
        {params.snakeopts} &> {log}
        """


use rule generate_summary_stats as generate_summary_stats_verkko with:
    input:
        fa=rules.run_verkko.output,
    output:
        summary=join("results", "verkko", "{sm}", "assembly_stats.tsv"),
    log:
        "logs/generate_summary_stats_{sm}_verkko.log",


rule verkko_all:
    input:
        rules.run_verkko.output,
        rules.generate_summary_stats_verkko.output,
