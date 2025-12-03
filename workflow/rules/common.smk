import subprocess

from enum import Enum
from typing import NamedTuple, Any
from collections import defaultdict
from os.path import join, dirname, split, splitext, basename, abspath


def get_data_config(sm: str) -> dict:
    return config["samples"][sm]["data"]


def get_dtype_config(sm: str, dtype: str) -> dict:
    return get_data_config(sm).get(dtype, {})


def multi_flags(*vals, pre_opt: str, post_opt: str = "") -> str:
    return " ".join(f"{pre_opt} {val} {post_opt}" for val in vals)


def dtype_glob(sm: str, dtype: str, *, which: str = "include") -> list[str]:
    return get_dtype_config(sm, dtype).get(which, [])


def find_sep_command(
    sm: str, dtype: str, indir: str, *, sep: str = " ", with_paste: bool = True
) -> str:
    indir = mat_dir = os.path.abspath(indir)
    include_flags = multi_flags(
        *dtype_glob(sm, dtype, which="include"), pre_opt="-name"
    )
    # https://stackoverflow.com/a/4210072
    exclude_flags = multi_flags(
        *dtype_glob(sm, dtype, which="exclude"),
        pre_opt=r"-not \( -name",
        post_opt=r"-prune \)",
    )
    cmd = f"find {indir}/ {include_flags} {exclude_flags}"
    if with_paste:
        cmd += " | paste -sd '{sep}'"
    return cmd
