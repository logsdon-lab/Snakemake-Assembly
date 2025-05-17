SAMPLES = config["samples"].keys()


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


DATA_DIR = join("data", "{dtype}", "{sm}")


wildcard_constraints:
    asm="|".join(a.value for a in Assembler),
    sm="|".join(SAMPLES),
    hap="|".join(h.value for h in Haplotype),
    dtype="|".join(d.value for d in DataType),
