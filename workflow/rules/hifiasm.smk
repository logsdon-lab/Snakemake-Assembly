# rule count_kmers:
#     input:

#     output:
#         "{hap}.yak"
#     params:
#         kmer_size="",
#     shell:
#         """
#         # For paired end provide, two identical streams.
#         yak count -b37 -t16 -o pat.yak <(cat pat_1.fq.gz pat_2.fq.gz) <(cat pat_1.fq.gz pat_2.fq.gz)
#         """

# def strand_data(wc):
#     if "illumina_mat" in DATA_DIRS[wc.sm] and "illumina_pat" in DATA_DIRS[wc.sm]:
#         return {"hap_kmers": expand(rules.generate_hapmers.output, sm=wc.sm)}
#     else:
#         return {}


# def strand_data_args(wc, inputs):
#     if "hap_kmers" in inputs.keys():
#         return f"--hap-kmers {inputs.hap_kmers} trio"
#     else:
#         return ""


rule run_hifiasm:
    input:
        # TODO: trios
        hic_dir=lambda wc: DATA_DIRS[wc.sm]["hic"],
        ont_dir=lambda wc: DATA_DIRS[wc.sm]["ont"],
        hifi_dir=lambda wc: DATA_DIRS[wc.sm]["hifi"],
    output:
        directory(join("results", "hifiasm", "{sm}")),
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
        # TODO: Split into hic1 and hic2 in config.
        hic_h1_rgx=".*R1.*fastq.gz",
        hic_h2_rgx=".*R2.*fastq.gz",
        hic_glob=lambda wc: multi_flags(*dtype_glob(str(wc.sm), "hic"), opt="-name"),
        ont_glob=lambda wc: multi_flags(*dtype_glob(str(wc.sm), "ont"), opt="-name"),
        hifi_glob=lambda wc: multi_flags(*dtype_glob(str(wc.sm), "hifi"), opt="-name"),
    shell:
        """
        logpath=$(realpath {log})
        hic_dir=$(realpath {input.hic_dir})
        ont_dir=$(realpath {input.ont_dir})
        hifi_dir=$(realpath {input.hifi_dir})
        mkdir -p {output} && cd {output}
        hifiasm \
        -o "{wildcards.sm}" \
        --h1 $(find ${{hic_dir}} {params.hic_glob} | grep -P {params.hic_h1_rgx} | paste -sd ",") \
        --h2 $(find ${{hic_dir}} {params.hic_glob} | grep -P {params.hic_h2_rgx} | paste -sd ",") \
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
