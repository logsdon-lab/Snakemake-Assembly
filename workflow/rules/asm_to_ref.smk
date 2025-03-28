ALN_CFG = {
    "ref": config["asm_to_ref"]["ref"],
    "sm": {
        sm: expand(
            (
                rules.run_verkko.output
                if assembler == "verkko"
                else rules.merge_haps_hifiasm.output
            ),
            sm=sm,
        )
        for sm in SAMPLES
    },
    "aln_threads": config["asm_to_ref"].get("threads", 8),
    "aln_mem": config["asm_to_ref"].get("mem", "60GB"),
    "temp_dir": config["asm_to_ref"].get("temp_dir", "results/asm_to_ref/temp"),
    "output_dir": config["asm_to_ref"].get("output_dir", "results/asm_to_ref"),
    "logs_dir": config["asm_to_ref"].get("log_dir", "logs/asm_to_ref"),
    "benchmarks_dir": config["asm_to_ref"].get("bmk_dir", "benchmarks/asm_to_ref"),
    "mm2_opts": config["asm_to_ref"].get(
        "mm2_opts", "-x asm20 --secondary=no -s 25000 -K 8G"
    ),
}


module asm_to_ref:
    snakefile:
        github(
            "koisland/asm-to-reference-alignment",
            path="workflow/Snakefile",
            branch="ideogram",
        )
    config:
        ALN_CFG


use rule * from asm_to_ref as asm_ref_*


for mode in config["asm_to_ref"].get("mode", ["saffire"]):
    if mode == "saffire":
        outputs.extend(rules.asm_ref_all.input)

    elif mode == "ideogram":
        outputs.extend(rules.asm_ref_ideogram.input)
