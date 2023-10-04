#!/usr/bin/env python3
import argparse
import pathlib
from typing import Optional

import pandas as pd
import yaml
import tqdm
import shutil


def parse_arguments(
    key: Optional[str] = None, stage: Optional[str] = None, method: Optional[str] = None
):
    """
    Parses command line arguments
    """
    parser = argparse.ArgumentParser(description="Run Tool")
    parser.add_argument(
        "--method",
        required=False,
        default=method,
        type=str,
        help="method key from params.yaml to use",
    )
    parser.add_argument(
        "--stage",
        required=False,
        default=stage,
        type=str,
        help="method key from params.yaml to use",
    )
    parser.add_argument(
        "--key",
        required=False,
        default=key,
        type=str,
        help="method key from params.yaml to use",
    )
    args = parser.parse_args()
    params = yaml.safe_load(pathlib.Path("params.yaml").read_text())[args.key][
        args.method
    ]
    settings = dict()
    settings.update(params.get("default_system_settings", {}))
    settings.update(params.get("default_settings", {}))
    settings.update(params["stages"][args.stage].get("settings", {}))
    return (
        settings,
        params["stages"][args.stage].get("inputs", {}),
        params["stages"][args.stage].get("output", None),
    )


def main():
    params, inputs, output = parse_arguments(key="rnaseq", method="rnaseq_merge")
    output = pathlib.Path(output)
    output.mkdir(exist_ok=True, parents=True)

    folder_old = pathlib.Path(inputs["old"])
    folder_new = pathlib.Path(inputs["new"])
    folders_old = pd.Series(folder_old.rglob("*.fastq.gz"))
    folders_new = pd.Series(folder_new.rglob("*.fastq.gz"))
    folders_old.index = (
        folders_old.map(lambda x: x.stem)
        .str.extract("A-(\w-[\w]+)_small_(R[12])\.fastq")
        .apply(lambda x: x[0] + "_RNA_" + x[1], axis=1)
    )
    folders_new.index = (
        folders_new.map(lambda x: x.stem)
        .str.extract("A-(\w-[\w]+)_REPEAT_(R[12])\.fastq")
        .apply(lambda x: x[0] + "_RNA_" + x[1], axis=1)
    )
    folders_old = folders_old.loc[folders_old.index.difference(folders_new.index)]
    files = pd.concat((folders_old, folders_new))
    for name, file in tqdm.tqdm(files.items(), total=len(files)):
        outfile = output.joinpath("merged_" + name).with_suffix(".fastq.gz")
        shutil.copyfile(file, outfile)


if __name__ == "__main__":
    main()
