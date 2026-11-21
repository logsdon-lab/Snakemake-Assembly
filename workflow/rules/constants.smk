from enum import Enum


SAMPLES = config["samples"].keys()
OUTPUT_DIR = config["output_dir"]
BENCHMARK_DIR = config["benchmark_dir"]
LOG_DIR = config["log_dir"]


class Haplotype(Enum):
    MAT = "mat"
    PAT = "pat"


class Assembler(Enum):
    VERKKO = "verkko"
    HIFIASM = "hifiasm"


class DataType(Enum):
    ONT = "ont"
    HIFI = "hifi"
    HIC_MAT = "hic_mat"
    HIC_PAT = "hic_pat"
    ILLUMINA_MAT = "illumina_mat"
    ILLUMINA_PAT = "illumina_pat"


# Store environment files for assemblers
ASM_ENVS = {}
for asm in Assembler:
    asm = asm.value
    ASM_ENVS[asm] = workflow.source_path(f"../envs/{asm}.yaml")
    ASM_ENVS[f"{asm}_template"] = workflow.source_path(f"../envs/{asm}_template.yaml")


wildcard_constraints:
    asm="|".join(a.value for a in Assembler),
    sm="|".join(SAMPLES),
    hap="|".join(h.value for h in Haplotype),
    dtype="|".join(d.value for d in DataType),
    mdata=r"[^\/]*",
