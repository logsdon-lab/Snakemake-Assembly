from os.path import join, dirname
from typing import NamedTuple
from collections import defaultdict
from enum import Enum


def get_data_config(sm: str) -> dict:
    return config["samples"][sm]["data"]


def get_dtype_config(sm: str, dtype: str) -> dict:
    return get_data_config(sm).get(dtype, {})


def multi_flags(*vals, opt: str) -> str:
    return " ".join(f"{opt} {val}" for val in vals)


def dtype_glob(sm: str, dtype: str, *, default=None) -> list[str]:
    if not default:
        default = ["*.fastq.gz"]
    return get_dtype_config(sm, dtype).get("include", default)
