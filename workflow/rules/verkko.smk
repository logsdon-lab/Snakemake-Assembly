HAPS = ["mat", "pat"]


rule count_kmers:
    input:
        illumina_dir=lambda wc: DATA_DIRS[wc.sm][f"illumina_{wc.hap}"],
    output:
        mer_db=directory(join("results", "meryl", "{sm}", "{hap}_compress.meryl")),
    resources:
        mem="30GB",
    threads: lambda wc: config["samples"][str(wc.sm)]["threads"] // 2
    log:
        join("logs", "meryl", "{sm}", "{hap}_count_kmers.log")
    benchmark:
        join("benchmarks", "meryl", "{sm}", "{hap}_count_kmers.tsv")
    params:
        fq_glob=lambda wc: multi_flags(*dtype_glob(str(wc.sm), f"illumina_{wc.hap}"), opt="-name"),
        kmer_size=lambda wc: config["samples"][str(wc.sm)]["data"][
            f"illumina_{wc.hap}"
        ].get("kmer_size", 30),
        mem=lambda wc, resources: str(resources.mem).replace("GB", "")
    conda:
        "../envs/verkko.yaml"
    shell:
        """
        meryl count compress \
        k={params.kmer_size} \
        threads={threads} \
        memory={params.mem} \
        $(find {input.illumina_dir}/ {params.fq_glob}) \
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
        join("logs", "meryl", "{sm}", "generate_hapmers.log")
    benchmark:
        join("benchmarks", "meryl", "{sm}", "generate_hapmers.tsv")
    params:
        output_dir=lambda wc, output: dirname(output.mat_db)
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
        # TODO: trios
        unpack(phasing_data),
        ont_dir=lambda wc: DATA_DIRS[wc.sm]["ont"],
        hifi_dir=lambda wc: DATA_DIRS[wc.sm]["hifi"],
    output:
        directory(join("results", "verkko", "{sm}")),
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
        ont_glob=lambda wc: multi_flags(*dtype_glob(str(wc.sm), "ont"), opt="-name"),
        hifi_glob=lambda wc: multi_flags(*dtype_glob(str(wc.sm), "hifi"), opt="-name"),
        phasing_data_args=lambda wc, input: phasing_data_args(wc, input),
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
        --hifi $(find {input.hifi_dir}/ {params.hifi_glob}) \
        --nano $(find {input.ont_dir}/ {params.ont_glob}) \
        {params.phasing_data_args} \
        {params.snakeopts}
        """


# rule convert_gfa_to_fa:
#     input:
#     output:
#     shell:
#         """
#         awk '/^S/{{print ">"$2;print $3}}'
#         """
