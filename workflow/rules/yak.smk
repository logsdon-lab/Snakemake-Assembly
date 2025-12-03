rule count_kmers_yak:
    input:
        illumina_dir=lambda wc: DATA_DIRS[wc.sm][f"illumina_{wc.hap}"],
    output:
        kmers=join(OUTPUT_DIR, "yak", "{sm}", "{hap}.yak"),
    params:
        fq_glob=lambda wc: multi_flags(
            *dtype_glob(str(wc.sm), f"illumina_{wc.hap}"), pre_opt="-name"
        ),
        kmer_size=lambda wc: config["samples"][str(wc.sm)]["data"][
            f"illumina_{wc.hap}"
        ].get("kmer_size", 31),
        bloom_filter_size=37,
    conda:
        "../envs/hifiasm.yaml"
    log:
        join(LOG_DIR, "yak", "{sm}", "{hap}_count_kmers_yak.log"),
    threads: lambda wc: config["samples"][str(wc.sm)]["threads"] // 2
    resources:
        mem=lambda wc: config["samples"][str(wc.sm)]["mem"],
    shell:
        """
        yak count \
        -k {params.kmer_size} \
        -b {params.bloom_filter_size} \
        -t {threads} \
        -o {output.kmers} \
        <(find {input.illumina_dir}/ {params.fq_glob} -size +0 -exec zcat -f {{}} +) 2> {log}
        """
