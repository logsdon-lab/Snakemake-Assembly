SAMPLES = config["samples"].keys()
DATA_DIR = join("data", "{dtype}", "{sm}")


wildcard_constraints:
    sm="|".join(SAMPLES),
