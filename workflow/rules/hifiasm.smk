rule count_kmers_yak:
    input:
        illumina_dir=lambda wc: DATA_DIRS[wc.sm][f"illumina_{wc.hap}"],
    output:
        kmers=join("results", "yak", "{sm}", "{hap}.yak"),
    params:
        fq_glob=lambda wc: multi_flags(
            *dtype_glob(str(wc.sm), f"illumina_{wc.hap}"), opt="-name"
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


def phasing_data_hifiasm(wc):
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
        raise ValueError("Not implemented or missing phasing data.")


def phasing_data_hifiasm_args(wc, inputs):
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
        raise ValueError("Not implemented or missing phasing data.")


rule run_hifiasm:
    input:
        unpack(phasing_data_hifiasm),
        ont_dir=lambda wc: DATA_DIRS[wc.sm]["ont"],
        hifi_dir=lambda wc: DATA_DIRS[wc.sm]["hifi"],
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
        output_dir=lambda wc, output: dirname(output[0]),
        ont_glob=lambda wc: multi_flags(*dtype_glob(str(wc.sm), "ont"), opt="-name"),
        hifi_glob=lambda wc: multi_flags(*dtype_glob(str(wc.sm), "hifi"), opt="-name"),
        phasing_data_args=lambda wc, input: phasing_data_hifiasm_args(wc, input),
    shell:
        """
        logpath=$(realpath {log})
        ont_dir=$(realpath {input.ont_dir})
        hifi_dir=$(realpath {input.hifi_dir})
        mkdir -p {params.output_dir} && cd {params.output_dir}
        hifiasm \
        -o "{wildcards.sm}" \
        {params.phasing_data_args} \
        --ul $(find ${{ont_dir}} {params.ont_glob} | paste -sd ",") \
        -t {threads} \
        $(find ${{hifi_dir}} {params.hifi_glob} | paste -sd " ") 2> "${{logpath}}"
        """


# rule convert_gfa_to_fa:
#     input:
#     output:
#     shell:
#         """
#         awk '/^S/{{print ">"$2;print $3}}'
#         """


rule hifiasm_all:
    input:
        rules.run_hifiasm.output,
        # rules.generate_summary_stats_verkko.output,
