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
  * **Required**
* `"hifi"`
  * **Required**
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
        include: ["*.fq.gz"]
        # Exclude files.
        exclude: ["*fail.fq.gz"]
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

### Analyses
Additional analyses can be added:

#### asm_to_ref
Align assembly to a reference genome.

```yaml
asm_to_ref:
  ref:
    CHM13: /project/logsdon_shared/projects/twins_chrY_assembly/data/reference/T2T-CHM13v2.fasta
  mm2_opts: "-x asm20 --secondary=no -s 25000 -K 8G"
  threads: 32
  mem: 250GB
  mode: ["saffire", "ideogram"]
```

One or more modes is possible.

##### Ideogram
Generate an ideogram of the assembly.

##### SafFire
Generate [`SafFire`](https://github.com/mrvollger/SafFire) beds.

#### NucFlag
> [!NOTE]
> TODO

### Examples
For more examples, see the `examples/` directory.
