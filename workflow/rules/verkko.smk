HAPS = ["mat", "pat"]


rule count_kmers:
    input:
        illumina_dir=lambda wc: DATA_DIRS[wc.sm][f"illumina_{wc.hap}"],
    output:
        mer_db=join("results", "meryl", "{sm}", "{hap}_compress.meryl"),
    resources:
        mem=30,
    threads: lambda wc: config["samples"][str(wc.sm)]["threads"] // 2
    params:
        fq_glob=lambda wc, input: join(str(input.illumina_dir), ".*fastq.gz"),
        kmer_size=lambda wc: config["samples"][str(wc.sm)]["data"][
            f"illumina_{wc.hap}"
        ].get("kmer_size", 30),
    conda:
        "../envs/verkko.yaml"
    shell:
        """
        meryl count compress k={params.kmer_size} threads={threads} memory={resources.mem} {params.fq_glob} output {output.mer_db}
        """


rule generate_hapmers:
    input:
        mat_kmers=lambda wc: expand(rules.count_kmers.output, sm=wc.sm, hap="mat"),
        pat_kmers=lambda wc: expand(rules.count_kmers.output, sm=wc.sm, hap="pat"),
    output:
        mat_db=join("results", "meryl", "{sm}", "mat_compress.hapmer.meryl"),
        pat_db=join("results", "meryl", "{sm}", "pat_compress.hapmer.meryl"),
    conda:
        "../envs/verkko.yaml"
    shell:
        """
        $MERQURY/trio/hapmers.sh \
            {input.mat_kmers} \
            {input.pat_kmers} \
        """


def strand_data(wc):
    if "illumina_mat" in DATA_DIRS[wc.sm] and "illumina_pat" in DATA_DIRS[wc.sm]:
        return {"hap_kmers": expand(rules.generate_hapmers.output, sm=wc.sm)}
    else:
        return {}


def strand_data_args(wc, inputs):
    if "hap_kmers" in inputs.keys():
        return f"--hap-kmers {inputs.hap_kmers} trio"
    else:
        return ""


rule run_verkko:
    input:
        # TODO: trios
        unpack(strand_data),
        ont_dir=lambda wc: DATA_DIRS[wc.sm]["ont"],
        hifi_dir=lambda wc: DATA_DIRS[wc.sm]["hifi"],
    output:
        directory(join("results", "verkko", "{sm}")),
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
        ont_glob=lambda wc: multi_flags(*dtype_glob(str(wc.sm), "ont"), opt="-name"),
        hifi_glob=lambda wc: multi_flags(*dtype_glob(str(wc.sm), "hifi"), opt="-name"),
        strand_args=lambda wc, input: strand_data_args(wc, input),
        snakeopts=lambda wc: "--snakeopts '"
        + " ".join(
            [
                "-j",
                        str(config["samples"][str(wc.sm)]["threads"]),
                        "-k",
                        "--restart-times",
                        "1",
                    ]
                )
        + "'",
    shell:
        """
        verkko -d {output} \
        --hifi $(find ${{hifi_dir}} {params.hifi_glob}) \
        --nano $(find ${{ont_dir}} {params.ont_glob}) \
        {params.strand_args} \
        {params.snakeopts}
        """


# rule convert_gfa_to_fa:
#     input:
#     output:
#     shell:
#         """
#         awk '/^S/{{print ">"$2;print $3}}'
#         """
