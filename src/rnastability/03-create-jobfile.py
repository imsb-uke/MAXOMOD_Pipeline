#!/usr/bin/env python3
import argparse
import pathlib
import subprocess
from typing import Optional

import pandas as pd
import yaml


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
    settings.update(params.get("default_system_settings"))
    settings.update(params.get("default_settings"))
    settings.update(params["stages"][args.stage].get("settings", {}))
    return (
        settings,
        params["stages"][args.stage].get("inputs", {}),
        params["stages"][args.stage].get("outputs", None),
    )


def main():
    params, inputs, outputs = parse_arguments(
        key="rnastability", method="rembrandts_jobfile"
    )
    outputs = pathlib.Path(outputs)
    outputs.parent.mkdir(exist_ok=True, parents=True)

    files = pd.Series(pathlib.Path(inputs).rglob("*.txt"))
    files.name = "File"
    files = files.to_frame()
    files = files[files.File.map(lambda x: x.stat().st_size) > 0]
    files["ReadType"] = files.File.map(lambda x: x.suffixes[0][1:])
    files["Label"] = (
        files.File.map(lambda x: x.stem).str.split(".", expand=True).iloc[:, 0]
    )
    files.ReadType = files.ReadType.astype("category")
    files["Batch"] = 1
    files.ReadType = files.ReadType.cat.rename_categories(
        {"intron": "intronic", "exon": "exonic"}
    )
    files = files[["Label", "File", "ReadType", "Batch"]]
    matched = files.Label.value_counts()
    matched = matched[matched == 2].index
    files =  files[files.Label.isin(matched)]
    files = files.sort_values(["Batch", "Label", "ReadType"])
    files.File = files.File.map(lambda x: str(x.relative_to(inputs)))
    files.to_csv(outputs, sep="\t", index=False)


if __name__ == "__main__":
    main()
