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


def dtype_glob(
    sm: str, dtype: str, *, which: str = "include", default=None
) -> list[str]:
    if not default:
        default = ["*.fastq.gz"]
    return get_dtype_config(sm, dtype).get(which, default)


# TODO: exclude https://stackoverflow.com/a/4210072
def find_sep_command(sm: str, dtype: str, indir: str, *, sep: str = " ") -> str:
    indir = mat_dir = os.path.abspath(indir)
    include_flags = multi_flags(*dtype_glob(sm, dtype, which="include"), opt="-name")
    return f"find {indir}/ {include_flags} | paste -sd '{sep}'"
