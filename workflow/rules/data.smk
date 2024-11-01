class DataType(StrEnum):
    ONT = auto()
    HIFI = auto()
    HIC = auto()


rule aws_sync:
    output:
        directory(join("data", "{dtype}", "{sm}")),
    params:
        uri=lambda wc: get_data_config(str(wc.sm)).get(str(wc.dtype), {})["uri"],
        include=lambda wc: multi_flags(
            *get_dtype_config(str(wc.sm), str(wc.dtype)).get("include", []),
            opt="--include",
        ),
        exclude=lambda wc: multi_flags(
            *get_dtype_config(str(wc.sm), str(wc.dtype)).get("exclude", []),
            opt="--exclude",
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


rule symlink_data:
    input:
        data_dir=lambda wc: get_data_config(str(wc.sm), str(wc.dtype))["path"],
    output:
        join("data", "{dtype}", "{sm}"),
    log:
        "logs/symlink_data_{dtype}_{sm}.log",
    shell:
        """
        ln -s {input} {output} 2> {log}
        """


ruleorder: symlink_data > aws_sync


def get_data_dirs(wc) -> list:
    outputs = []
    sm = str(wc.sm)
    for dtype, settings in get_data_config(sm).items():
        dtype = DataType(dtype)

        if settings.get("uri"):
            output_dir = expand(rules.aws_sync.output, sm=sm, dtype=dtype)
        elif settings.get("path"):
            output_dir = expand(rules.symlink_data.output, sm=sm, dtype=dtype)

        if isinstance(output_dir, str):
            outputs.append(output_dir)
        elif isinstance(output_dir, list):
            outputs.extend(output_dir)
        else:
            raise ValueError(f"Invalid {output_dir} for {sm}")

    return outputs


rule get_input_dirs:
    input:
        get_data_dirs,
    output:
        touch(temp("/tmp/input_dirs_{sm}.done")),


rule get_input_dirs_only:
    input:
        expand(rules.get_input_dirs.output, sm=SAMPLES),
