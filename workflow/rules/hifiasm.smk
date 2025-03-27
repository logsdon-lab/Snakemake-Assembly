rule count_kmers_yak:
    input:
        illumina_dir=lambda wc: DATA_DIRS[wc.sm][f"illumina_{wc.hap}"],
    output:
        kmers=join("results", "yak", "{sm}", "{hap}.yak"),
    params:
        fq_glob=lambda wc: multi_flags(
            *dtype_glob(str(wc.sm), f"illumina_{wc.hap}"), pre_opt="-name"
        ),
        kmer_size=lambda wc: config["samples"][str(wc.sm)]["data"][
            f"illumina_{wc.hap}"
        ].get("kmer_size", 31),
        bloom_filter_size=37,
    log:
        join("logs", "yak", "{sm}", "{hap}_count_kmers_yak.log"),
    threads: lambda wc: config["samples"][str(wc.sm)]["threads"] // 2
    resources:
        mem=lambda wc: config["samples"][str(wc.sm)]["mem"] // 2,
    shell:
        """
        yak count \
        -k {params.kmer_size} \
        -b {params.bloom_filter_size} \
        -t {threads} \
        -o {output.kmers} \
        <(find {input.illumina_dir}/ {params.fq_glob} -size +0 -exec zcat -f {{}} +) 2> {log}
        """


def phasing_data_hifiasm(wc) -> dict[str, Any]:
    if "illumina_mat" in DATA_DIRS[wc.sm] and "illumina_pat" in DATA_DIRS[wc.sm]:
        return {
            "mat_kmers": expand(rules.count_kmers_yak.output, sm=wc.sm, hap="mat"),
            "pat_kmers": expand(rules.count_kmers_yak.output, sm=wc.sm, hap="pat"),
        }
    elif "hic_mat" in DATA_DIRS[wc.sm] and "hic_pat" in DATA_DIRS[wc.sm]:
        return {
            "mat_hic_dir": DATA_DIRS[wc.sm]["hic_mat"],
            "pat_hic_dir": DATA_DIRS[wc.sm]["hic_pat"],
        }
    else:
        return {}


def phasing_data_hifiasm_args(wc, inputs) -> str:
    if "mat_kmers" in inputs.keys() and "pat_kmers" in inputs.keys():
        return f"-1 {inputs.mat_kmers} -2 {inputs.pat_kmers}"
    elif "mat_hic_dir" in inputs.keys() and "pat_hic_dir" in inputs.keys():
        # Construct find command.
        mat_find_cmd = find_sep_command(
            wc.sm, "hic_mat", inputs.mat_hic_dir[0], sep=","
        )
        pat_find_cmd = find_sep_command(
            wc.sm, "hic_pat", inputs.pat_hic_dir[0], sep=","
        )
        return f"--h1 $({mat_find_cmd}) --h2 $({pat_find_cmd})"
    else:
        return ""


checkpoint run_hifiasm:
    input:
        unpack(phasing_data_hifiasm),
        ont_fofn=expand(rules.generate_dtype_fofn.output, asm="hifiasm", dtype="ont", sm="{sm}"),
        hifi_fofn=expand(rules.generate_dtype_fofn.output, asm="hifiasm", dtype="hifi", sm="{sm}"),
    output:
        # GFA name changes based on phasing data so cannot get path.
        # Don't use output directory as will delete directory if fail.
        touch(join("results", "hifiasm", "{sm}.done")),
    conda:
        "../envs/hifiasm.yaml"
    threads: lambda wc: config["samples"][str(wc.sm)]["threads"]
    resources:
        mem=lambda wc: config["samples"][str(wc.sm)]["mem"],
        lsf_extra=lambda wc: "-M " + config["samples"][str(wc.sm)]["mem"],
    log:
        "logs/run_hifiasm_{sm}.log",
    benchmark:
        "benchmarks/run_hifiasm_{sm}.tsv"
    params:
        output_dir=lambda wc, output: splitext(output[0])[0],
        phasing_data_args=lambda wc, input: phasing_data_hifiasm_args(wc, input),
    shell:
        """
        logpath=$(realpath {log})
        ont_fofn_path=$(realpath {input.ont_fofn})
        hifi_fofn_path=$(realpath {input.hifi_fofn})
        mkdir -p {params.output_dir} && cd {params.output_dir}
        hifiasm \
        -o "{wildcards.sm}" \
        {params.phasing_data_args} \
        --ul $(paste -sd "," "${{ont_fofn_path}}") \
        -t {threads} \
        $(paste -sd " " "${{hifi_fofn_path}}") 2> "${{logpath}}"
        """


checkpoint convert_gfa_to_fa:
    input:
        join("results", "hifiasm", "{sm}", "{sm}.{mdata}.p_ctg.gfa"),
    output:
        fa=join("results", "hifiasm", "{sm}", "{sm}.{mdata}.p_ctg.fa"),
        fai=join("results", "hifiasm", "{sm}", "{sm}.{mdata}.p_ctg.fa.fai"),
    log:
        "logs/convert_gfa_to_fa_{sm}_{mdata}.log",
    conda:
        "../envs/hifiasm.yaml"
    shell:
        """
        awk '/^S/{{ print ">"$2; print $3 }}' {input} > {output.fa} 2> {log}
        samtools faidx {output.fa} 2>> {log}
        """


def primary_contigs(wc):
    chkpt = checkpoints.run_hifiasm.get(**wc).output[0]
    output_dir = splitext(chkpt)[0]
    # Only look at primary contigs
    # https://hifiasm.readthedocs.io/en/latest/interpreting-output.html#interpreting-output
    wcs = glob_wildcards(join(output_dir, f"{wc.sm}.{{mdata}}.p_ctg.gfa"))
    mdata = wcs.mdata

    wildcard_constraints:
        mdata="|".join(mdata),

    return expand(rules.convert_gfa_to_fa.output.fa, sm=wc.sm, mdata=mdata)


def hap_contigs(wc):
    chkpt = checkpoints.run_hifiasm.get(**wc).output[0]
    output_dir = splitext(chkpt)[0]
    # Only look at primary contigs
    # https://hifiasm.readthedocs.io/en/latest/interpreting-output.html#interpreting-output
    wcs = glob_wildcards(join(output_dir, f"{wc.sm}.{{mdata}}.p_ctg.gfa"))

    return [
        checkpoints.convert_gfa_to_fa.get(**wc, mdata=mdata).output.fa
        for mdata in wcs.mdata
        if "hap1" in mdata or "hap2" in mdata
    ]


rule merge_haps:
    input:
        fa=hap_contigs
    output:
        fa=join("results", "hifiasm", "{sm}", "assembly.fasta"),
    log:
        "logs/merge_haps_{sm}.log",
    conda:
        "../envs/hifiasm.yaml"
    shell:
        """
        cat {input.fa} > {output.fa}
        samtools faidx {output.fa} 2> {log}
        """


use rule generate_summary_stats as generate_summary_stats_hifiasm with:
    input:
        fa=lambda wc: [*primary_contigs(wc), *expand(rules.merge_haps.output, sm=wc.sm)],
    output:
        summary=join("results", "hifiasm", "{sm}", "assembly_stats.tsv"),
    log:
        "logs/generate_summary_stats_{sm}_hifiasm.log",


rule hifiasm_all:
    input:
        rules.run_hifiasm.output,
        rules.merge_haps.output,
        rules.generate_summary_stats_hifiasm.output,
