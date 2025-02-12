class DataType(Enum):
    ONT = "ont"
    HIFI = "hifi"
    HIC_MAT = "hic_mat"
    HIC_PAT = "hic_pat"
    ILLUMINA_MAT = "illumina_mat"
    ILLUMINA_PAT = "illumina_pat"


rule aws_sync:
    output:
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
        "../envs/aws.yaml"
    shell:
        """
        aws s3 --no-sign-request sync {params.uri} {output} \
        {params.include} {params.exclude} &> {log}
        """


def get_data_dirs() -> list:
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
