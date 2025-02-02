SAMPLES = config["samples"].keys()
HAPS = ["mat", "pat"]
DATA_DIR = join("data", "{dtype}", "{sm}")


wildcard_constraints:
    sm="|".join(SAMPLES),
    hap="|".join(HAPS),
