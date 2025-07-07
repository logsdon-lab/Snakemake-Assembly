
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
        "logs/download_{dtype}_{sm}.log",
    benchmark:
        "benchmarks/download_{dtype}_{sm}.tsv"
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


checkpoint generate_original_dtype_fofn:
    input:
        dtype_dir=lambda wc: DATA_DIRS[wc.sm][wc.dtype],
    output:
        fofn=join("results", "{asm}", "{sm}_{dtype}_original.fofn"),
    params:
        glob=lambda wc: multi_flags(*dtype_glob(str(wc.sm), wc.dtype), pre_opt="-name"),
    shell:
        """
        find {input.dtype_dir}/ {params.glob} > {output.fofn}
        """


def check_for_bams(wc):
    fofn = checkpoints.generate_original_dtype_fofn.get(**wc).output
    all_files = []
    with open(fofn[0], "rt") as fh:
        for file in fh.readlines():
            file = file.strip()
            indir, fname = split(file)
            if fname.endswith(".bam"):
                bname, ext = splitext(fname)
                new_fname = join(indir, bname)
                all_files.append(
                    expand(rules.convert_bam_to_fastq.output, fname=new_fname)
                )
            elif fname.endswith(".tmp.fastq"):
                raise ValueError(f"{fname} contains a reserved suffix (.tmp.fastq). Please change this to something else.")
            else:
                all_files.append(file)
    return all_files


# Issue with gzip so uncompressed.
# https://github.com/marbl/verkko/issues/320
rule convert_bam_to_fastq:
    input:
        bam="{fname}.bam",
    output:
        fastq="{fname}.tmp.fastq",
    threads: 8
    log:
        "{fname}.to_fastq.log",
    conda:
        "../envs/data.yaml"
    shell:
        """
        samtools bam2fq --threads {threads} {input.bam} > {output.fastq} 2> {log}
        """


rule generate_dtype_fofn:
    input:
        unpack(check_for_bams),
    output:
        join("results", "{asm}", "{sm}_{dtype}_final.fofn"),
    shell:
        """
        ls {input} > {output}
        """
