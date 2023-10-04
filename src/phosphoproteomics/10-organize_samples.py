#!/usr/bin/env python3

import argparse
import pathlib
from typing import Optional
import yaml

import pandas as pd


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
        params["stages"][args.stage].get("output", None),
    )


def main():
    settings, inputs, outputs = parse_arguments(
        key="phosphoproteomics", method="organize_samples"
    )

    outdir = pathlib.Path(outputs)
    outdir.mkdir(exist_ok=True, parents=True)

    data = pd.read_excel(inputs["inputs"], index_col=0)

    mapping = pd.read_csv(inputs["samplesheet"], sep="\t", index_col=0)

    data = data.T

    data.index = data.index.str.replace(
        settings["intensity_column_prefix"], ""
    ).str.replace("_", "-")

    data = pd.merge(data, mapping, left_index=True, right_index=True).set_index(
        "SampleID"
    )

    cohort = pd.read_csv(inputs["cohort"], index_col=0)

    exclude = settings["exclude_samples"]

    if exclude:
        cohort = cohort[~cohort.index.isin(exclude)]

    data = data.reindex(cohort.index)

    cohort.to_csv(outdir.joinpath("cohort_filtered.csv"))
    data.to_csv(outdir.joinpath("intensity_mat_filtered_imputed_log2transf_norm.csv"))


if __name__ == "__main__":
    main()
