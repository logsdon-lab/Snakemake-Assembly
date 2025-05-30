CFG_NUCFLAG_SAMPLES = {sm["name"]: sm for sm in config["nucflag"].get("samples", [])}


def get_nucflag_inputs(sm: str) -> dict[str, Any]:
    inputs = {
        "name": sm,
        "asm_fa": expand(
            (
                rules.verkko_output.output
                if config["samples"][sm]["assembler"] == "verkko"
                else rules.hifiasm_output.output
            ),
            sm=sm,
            asm=config["samples"][sm]["assembler"],
        )[0],
    }
    CFG_NUCFLAG_SAMPLES_SM = CFG_NUCFLAG_SAMPLES.get(sm, {})
    if CFG_NUCFLAG_SAMPLES_SM.get("read_dir"):
        inputs = {
            **inputs,
            "read_dir": join(CFG_NUCFLAG_SAMPLES_SM["read_dir"], sm),
            "read_rgx": CFG_NUCFLAG_SAMPLES_SM["read_rgx"],
            **CFG_NUCFLAG_SAMPLES_SM,
        }
    elif CFG_NUCFLAG_SAMPLES_SM.get("read_fofn"):
        inputs = {
            **inputs,
            "read_fofn": join(CFG_NUCFLAG_SAMPLES_SM["read_fofn"], f"{sm}.fofn"),
            **CFG_NUCFLAG_SAMPLES_SM,
        }
    else:
        cfg = get_dtype_config(sm, "hifi")
        args = [
            arg
            for arg in find_sep_command(
                sm, "hifi", cfg["path"], with_paste=False
            ).split(" ")
            if arg
        ]
        # We cannot past hifi fofn because has not been created yet.
        # Must find files and pass.
        # TODO: Currently incompatible with S3 hifi input.
        # * Would need to setup checkpoint and wait for files to download. Then could use hifi fofn.
        #  expand(
        #       rules.generate_dtype_fofn.output,
        #       zip,
        #       asm=wc.asm,
        #       dtype="hifi",
        #       sm=wc.sm,
        #   )
        out = subprocess.run(args, check=True, capture_output=True, text=True)
        files = [file for file in out.stdout.split("\n") if file]
        if out.returncode != 0 or not files:
            raise FileNotFoundError(
                f"Unable to find reads for {sm} with command: {args.join(' ')}"
            )

        inputs = {**inputs, "reads": files}
    return inputs


NUCFLAG_CFG = {
    "samples": [get_nucflag_inputs(sm) for sm in SAMPLES],
    # Dump other config.
    **{k: v for k, v in config["nucflag"].items() if k != "samples"},
}


module NucFlag:
    snakefile:
        "Snakemake-NucFlag/workflow/Snakefile"
    config:
        NUCFLAG_CFG


use rule * from NucFlag


outputs.extend(rules.nucflag.input)
