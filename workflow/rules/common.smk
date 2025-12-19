import os
import yaml
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


# We want to be able to set version programmatically for testing multiple versions.
def dynamic_assembler_conda_env(wc) -> str:
    if not config["samples"][wc.sm].get("version"):
        return workflow.source_path(f"workflow/envs/{wc.asm}.yaml")

    version = config["samples"][wc.sm]["version"]
    template_conda_env = workflow.source_path(f"workflow/envs/{wc.asm}_template.yaml")
    output_conda_dir = join(OUTPUT_DIR, wc.asm, wc.sm)
    output_conda_env = abspath(join(output_conda_dir, "env.yaml"))
    os.makedirs(output_conda_dir, exist_ok=True)
    with open(template_conda_env, "rb") as fh, open(output_conda_env, "wt") as ofh:
        conda_env = yaml.safe_load(fh)
        conda_env["dependencies"].remove(wc.asm)
        conda_env["dependencies"].append(f"{wc.asm}=={version}")
        yaml.safe_dump(conda_env, ofh)

    return output_conda_env
