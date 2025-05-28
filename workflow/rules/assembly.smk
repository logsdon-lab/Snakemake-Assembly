include: "yak.smk"
include: "meryl.smk"


def phasing_data(wc) -> dict[str, Any]:
    if "illumina_mat" in DATA_DIRS[wc.sm] and "illumina_pat" in DATA_DIRS[wc.sm]:
        if wc.asm == "verkko":
            return {"hap_kmers": expand(rules.generate_hapmers.output, sm=wc.sm)}
        else:
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


def phasing_data_args(wc, inputs) -> str:
    if "hap_kmers" in inputs.keys():
        return f"--hap-kmers {inputs.hap_kmers} trio"
    elif "mat_kmers" in inputs.keys() and "pat_kmers" in inputs.keys():
        mat_kmers = abspath(inputs.mat_kmers[0])
        pat_kmers = abspath(inputs.pat_kmers[0])
        return f"-1 {mat_kmers} -2 {pat_kmers}"
    elif "mat_hic_dir" in inputs.keys() and "pat_hic_dir" in inputs.keys():
        mat_hic_dir = abspath(inputs.mat_hic_dir[0])
        pat_hic_dir = abspath(inputs.pat_hic_dir[0])
        # Construct find command.
        mat_find_cmd = find_sep_command(wc.sm, "hic_mat", mat_hic_dir, sep=",")
        pat_find_cmd = find_sep_command(wc.sm, "hic_pat", pat_hic_dir, sep=",")
        if wc.asm == "verkko":
            return f"--hic1 $({mat_find_cmd}) --hic2 $({pat_find_cmd})"
        else:
            return f"--h1 $({mat_find_cmd}) --h2 $({pat_find_cmd})"
    else:
        return ""


def input_args(wc, inputs) -> str:
    if wc.asm == "hifiasm":
        args = []
        if inputs["ont_fofn"] and not inputs["hifi_fofn"]:
            ont_fofn = abspath(inputs.ont_fofn[0])
            # ONT only assembly.
            # https://github.com/chhylp123/hifiasm?tab=readme-ov-file#assembling-ont-reads
            args.append("--ont")
            args.append(f'$(paste -sd " " {ont_fofn})')
        elif inputs["ont_fofn"]:
            ont_fofn = abspath(inputs.ont_fofn[0])
            args.append(f'--ul $(paste -sd "," {ont_fofn})')
        if inputs["hifi_fofn"]:
            hifi_fofn = abspath(inputs.hifi_fofn[0])
            args.append(f'$(paste -sd " " {hifi_fofn})')
        args_str = " ".join(args)
        if not args_str:
            raise ValueError(f"Require at least ont or hifi data for {wc.sm}.")
    else:
        args = ["--hifi", f"$(cat {inputs.hifi_fofn})"]
        if inputs["ont_fofn"]:
            args.append("--nano")
            args.append(f"$(cat {inputs.ont_fofn})")
        args_str = " ".join(args)
    return args_str


checkpoint run_assembler:
    input:
        unpack(phasing_data),
        ont_fofn=lambda wc: (
            expand(
                rules.generate_dtype_fofn.output,
                zip,
                asm=wc.asm,
                dtype="ont",
                sm=wc.sm,
            )
            if config["samples"][str(wc.sm)]["data"].get("ont")
            else []
        ),
        hifi_fofn=lambda wc: (
            expand(
                rules.generate_dtype_fofn.output,
                zip,
                asm=wc.asm,
                dtype="hifi",
                sm=wc.sm,
            )
            if config["samples"][str(wc.sm)]["data"].get("hifi")
            else []
        ),
    output:
        # hifiasm:
        # GFA name changes based on phasing data so cannot get path.
        # Don't use output directory as will delete directory if fail.
        touch(join("results", "{asm}", "{sm}.done")),
    conda:
        "../envs/{asm}.yaml"
    threads: lambda wc: config["samples"][str(wc.sm)]["threads"]
    resources:
        mem=lambda wc: config["samples"][str(wc.sm)]["mem"],
        lsf_extra=lambda wc: "-M " + config["samples"][str(wc.sm)]["mem"],
    log:
        abspath("logs/run_{asm}_{sm}.log"),
    benchmark:
        "benchmarks/run_{asm}_{sm}.tsv"
    params:
        added_args=lambda wc: config["samples"][str(wc.sm)].get("added_args", ""),
        output_dir=lambda wc, output: splitext(output[0])[0],
        phasing_data_args=lambda wc, input: phasing_data_args(wc, input),
        input_args=lambda wc, input: input_args(wc, input),
    shell:
        """
        if [[ "{wildcards.asm}" == "verkko" ]]; then
            verkko -d {params.output_dir} \
            {params.input_args} \
            {params.phasing_data_args} \
            {params.added_args} &> {log}
        else
            mkdir -p {params.output_dir} && cd {params.output_dir}
            hifiasm \
            -t {threads} \
            -o "{wildcards.sm}" \
            {params.phasing_data_args} \
            {params.input_args} \
            {params.added_args} 2> {log}
        fi
        """


checkpoint convert_gfa_to_fa:
    input:
        join("results", "{asm}", "{sm}", "{sm}.{mdata}.p_ctg.gfa"),
    output:
        fa=join("results", "{asm}", "{sm}", "{sm}.{mdata}.p_ctg.fa"),
        fai=join("results", "{asm}", "{sm}", "{sm}.{mdata}.p_ctg.fa.fai"),
    log:
        "logs/convert_gfa_to_fa_{asm}_{sm}_{mdata}.log",
    wildcard_constraints:
        asm="hifiasm",
    conda:
        "../envs/{asm}.yaml"
    shell:
        """
        awk '/^S/{{ print ">"$2; print $3 }}' {input} > {output.fa} 2> {log}
        samtools faidx {output.fa} 2>> {log}
        """


def primary_contigs(wc):
    chkpt = checkpoints.run_assembler.get(**wc).output[0]
    output_dir = splitext(chkpt)[0]
    # Only look at primary contigs
    # https://hifiasm.readthedocs.io/en/latest/interpreting-output.html#interpreting-output
    wcs = glob_wildcards(join(output_dir, f"{wc.sm}.{{mdata}}.p_ctg.gfa"))
    mdata = wcs.mdata

    return expand(rules.convert_gfa_to_fa.output.fa, **wc, mdata=mdata)


def hap_contigs(wc):
    chkpt = checkpoints.run_assembler.get(**wc).output[0]
    output_dir = splitext(chkpt)[0]
    # Only look at primary contigs
    # https://hifiasm.readthedocs.io/en/latest/interpreting-output.html#interpreting-output
    wcs = glob_wildcards(join(output_dir, f"{wc.sm}.{{mdata}}.p_ctg.gfa"))

    return [
        checkpoints.convert_gfa_to_fa.get(**wc, mdata=mdata).output.fa
        for mdata in wcs.mdata
        if "hap1" in mdata or "hap2" in mdata
    ]


rule hifiasm_output:
    input:
        fa=hap_contigs,
    output:
        fa=join("results", "{asm}", "{sm}", "assembly.fasta"),
    log:
        "logs/merge_haps_{asm}_{sm}.log",
    wildcard_constraints:
        asm="hifiasm",
    conda:
        "../envs/{asm}.yaml"
    shell:
        """
        cat {input.fa} > {output.fa}
        samtools faidx {output.fa} 2> {log}
        """


rule verkko_output:
    input:
        join("results", "{asm}", "{sm}.done"),
    output:
        join("results", "{asm}", "{sm}", "assembly.fasta"),
    wildcard_constraints:
        asm="verkko",


def asm_output(wc):
    outputs = []
    if wc.asm == "verkko":
        outputs.extend(expand(rules.verkko_output.output, sm=wc.sm, asm=wc.asm))
    else:
        outputs.extend(primary_contigs(wc))
        outputs.extend(expand(rules.hifiasm_output.output, sm=wc.sm, asm=wc.asm))
    return outputs


rule generate_summary_stats:
    input:
        fa=asm_output,
    output:
        summary=join("results", "{asm}", "{sm}", "assembly_stats.tsv"),
    log:
        "logs/generate_summary_stats_{sm}_{asm}.log",
    params:
        script=workflow.source_path("../scripts/seq_stats.py"),
    threads: 1
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {params.script} {input.fa} -r -t {threads} > {output} 2> {log}
        """


rule asm_all:
    input:
        rules.run_assembler.output,
        rules.generate_summary_stats.output,
