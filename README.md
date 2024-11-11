# Snakemake-Assembly
Workflow to generate `hifiasm` and `verkko` assemblies.

### Why?
No clear, simple alternatives.
Need basic AWS and local path support.


### Usage
```bash
snakemake -np --use-conda --configfile config.yaml --workflow-profile none
```