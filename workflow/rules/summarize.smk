
rule generate_summary_stats:
    input:
        fa="",
    output:
        summary="",
    log:
        "logs/generate_summary_stats.log",
    params:
        script=workflow.source_path("../scripts/seq_stats.py"),
    threads: 1
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {params.script} {input.fa} -r -t {threads} > {output} 2> {log}
        """
