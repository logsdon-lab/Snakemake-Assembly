
rule aws_sync:
    output:
        # NOTE: This will delete files on failed runs. Probably better way to trigger reruns.
        directory(join("data", "{dtype}", "{sm}")),
    params:
        uri=lambda wc: get_data_config(str(wc.sm)).get(str(wc.dtype), {}).get("uri"),
        include=lambda wc: multi_flags(
            *get_dtype_config(str(wc.sm), str(wc.dtype)).get("include", []),
            pre_opt="--include",
        ),
        exclude=lambda wc: multi_flags(
            *get_dtype_config(str(wc.sm), str(wc.dtype)).get("exclude", []),
            pre_opt="--exclude",
        ),
    threads: 1
    log:
        join(LOG_DIR, "download_{dtype}_{sm}.log"),
    benchmark:
        join(BENCHMARK_DIR, "download_{dtype}_{sm}.tsv")
    conda:
        "../envs/data.yaml"
    shell:
        """
        aws s3 --no-sign-request sync {params.uri} {output} \
        {params.include} {params.exclude} &> {log}
        """


def get_data_dirs() -> defaultdict[str, defaultdict[str, list[str]]]:
    outputs = defaultdict(lambda: defaultdict(list))
    for sm, sm_settings in config["samples"].items():
        for dtype, settings in sm_settings["data"].items():
            try:
                DataType(dtype)
            except ValueError:
                raise ValueError(f"Invalid datatype: {dtype}")

            if settings.get("uri"):
                output_dir = expand(rules.aws_sync.output, sm=sm, dtype=dtype)
            elif settings.get("path"):
                output_dir = settings["path"]

            if isinstance(output_dir, str):
                outputs[sm][dtype].append(output_dir)
            elif isinstance(output_dir, list):
                outputs[sm][dtype].extend(output_dir)
            else:
                raise ValueError(f"Invalid {output_dir} for {sm}")

    return outputs


# Global state.
DATA_DIRS = get_data_dirs()


rule generate_original_dtype_fofn:
    input:
        dtype_dir=lambda wc: DATA_DIRS[wc.sm][wc.dtype],
    output:
        fofn=join(OUTPUT_DIR, "{asm}", "{sm}_{dtype}_original.fofn"),
    params:
        glob=lambda wc: multi_flags(*dtype_glob(str(wc.sm), wc.dtype), pre_opt="-name"),
    shell:
        """
        find {input.dtype_dir}/ {params.glob} > {output.fofn}
        """


# Issue with gzip so uncompressed.
# https://github.com/marbl/verkko/issues/320
rule generate_dtype_fofn:
    input:
        rules.generate_original_dtype_fofn.output,
    output:
        join(OUTPUT_DIR, "{asm}", "{sm}_{dtype}_final.fofn"),
    threads: 16
    resources:
        processes=4,
        threads_per_file=4,
        tmp_dir=lambda wc: abspath(join(OUTPUT_DIR, "tmp")),
    conda:
        "../envs/data.yaml"
    log:
        join(LOG_DIR, "generate_dtype_fofn_{asm}_{dtype}_{sm}.log"),
    shell:
        """
        mkdir -p {resources.tmp_dir}
        cat {input} | xargs -P {resources.processes} -I {{}} bash -c '
            if [[ {{}} == *.bam ]]; then
                md5_checksum=$(md5sum {{}} | cut -c -32)
                tmp_outfq="{resources.tmp_dir}/${{md5_checksum}}.fastq";
                if [[ -f "${{tmp_outfq}}" ]]; then
                    echo "${{tmp_outfq}}";
                    exit 0;
                fi
                samtools bam2fq --threads {resources.threads_per_file} {{}} > "${{tmp_outfq}}" 2> {log}
                # Index file to ensure non-empty
                samtools faidx "${{tmp_outfq}}" -o /dev/null
                echo "${{tmp_outfq}}"
            else
                echo "{{}}"
            fi
        ' > {output}
        """
