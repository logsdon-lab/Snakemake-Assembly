# Snakemake-Assembly
Workflow to generate [`hifiasm`](https://github.com/chhylp123/hifiasm) and [`verkko`](https://github.com/marbl/verkko) assemblies.

## Why?
No well-documented, simple `Snakemake` workflows.
Need basic AWS and local path support.


## Usage
```bash
snakemake -np --use-conda --configfile config.yaml --workflow-profile none
```


## Config
Each sample is contained with a block in `samples`.
```yaml
samples:
  sample_name:
    threads: ... # Number of threads
    mem: ... # In GB. ex. "200GB"
    assembler: ... # Assembler. Either "verkko" or "hifiasm"
    data: ...
```

### Assembler
Either:
* `verkko`
* `hifiasm`

See `workflow/envs/(verkko|hifiasm).yaml` for version information.

### Data

#### Types
The following data types are supported for `{sm}.data.{dtype}`.
* `"ont"`
* `"hifi"`
* `"hic_mat"`
* `"hic_pat"`
* `"illumina_mat"`
* `"illumina_pat"`

#### Sources
Data sources can be either local or on AWS:

`path`
* `{sm}.data.{dtype}.path` will get data from local directory.

`uri`
* `{sm}.data.{dtype}.uri` will `aws sync` from the specified S3 uri.

##### Local
```yaml
samples:
  mPanTro3:
    threads: 40
    mem: 250GB
    assembler: hifiasm
    data:
      ont:
        path: /project/logsdon_shared/data/PrimateT2T/ont/mPanTro3
        # Include files to use.
        # Exclude not currently supported.
        include: ["*.fq.gz"]
      hifi:
        path: /project/logsdon_shared/data/PrimateT2T/hifi_data/mPanTro3
        include: ["*.hifi_reads.fq.gz"]
```

##### S3
```yaml
samples:
  mPanTro3:
    threads: 32
    mem: 250GB
    assembler: hifiasm
    data:
      ont:
        uri: s3://genomeark/species/Pan_troglodytes/mPanTro3/genomic_data/ont/
        # Include files
        include: ["*.fq.gz"]
        # Exclude files to download if include not specific enough.
        exclude: ["*old-guppy-runs/*", "*.bam*", "*fast5/*"]
      hifi:
        uri: s3://genomeark/species/Pan_troglodytes/mPanTro3/genomic_data/pacbio_hifi/
        include: ["*.hifi_reads.fq.gz"]
        exclude: ["*previous-versions/*", "*.bam*", "*ccs*"]
```

### Examples
For more examples, see the `examples/` directory.
