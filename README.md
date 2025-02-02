# Snakemake-Assembly
Workflow to generate `hifiasm` and `verkko` assemblies.

### Why?
No clear, simple alternatives.
Need basic AWS and local path support.


### Usage
```bash
snakemake -np --use-conda --configfile config.yaml --workflow-profile none
```

### Config
Data source either locally or AWS:
* `{sm}.data.{dtype}.path` will symlink to directory.
* `{sm}.data.{dtype}.uri` with sync from the specified S3 uri.

```yaml
samples:
  mPanTro3:
    threads: 40
    mem: 250GB
    assembler: hifiasm
    data:
      ont:
        path: /project/logsdon_shared/data/PrimateT2T/ont/mPanTro3
        include: ["*.fq.gz"]
        exclude: ["*fast5/*"]
      hifi:
        path: /project/logsdon_shared/data/PrimateT2T/hifi_data/mPanTro3
        include: ["*.hifi_reads.fq.gz"]
        exclude: ["*.bam*"]
      hic:
        path: /project/logsdon_shared/data/PrimateT2T/dovetail_hic/mPanTro3
        include: ["*.fastq.gz"]
```
