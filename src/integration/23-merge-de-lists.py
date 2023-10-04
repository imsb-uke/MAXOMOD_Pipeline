#!/usr/bin/env python3

import pathlib
import pandas as pd
import tqdm
import argparse
from typing import Optional

import numpy as np
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
        default=key,
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
        default=method,
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
        params["stages"][args.stage].get("output", None),
    )


def read_and_process(file):
    tmp = pd.read_csv(file, sep=",", index_col=0)
    tmp[["MODEL", "OMIC", "SEX"]] = file.stem.split("_")
    return tmp


def main():
    params, inputs, outputs = parse_arguments()
    files = list(pathlib.Path(inputs).glob("*.csv"))
    results = (read_and_process(file) for file in tqdm.tqdm(files))
    pathlib.Path(outputs).parent.mkdir(exist_ok=True, parents=True)
    results = pd.concat(results)
    results.MODEL = results.MODEL.str.upper()
    results.SEX = results.SEX.str.replace("males", "male")
    results["REGULATION"] = pd.NA
    results.loc[results.log2FoldChange < 0, "REGULATION"] = "down"
    results.loc[results.log2FoldChange > 0, "REGULATION"] = "up"
    if params["filter"]:
        results = results.query(params["filter"])
    results.to_csv(outputs)


if __name__ == "__main__":
    main()
