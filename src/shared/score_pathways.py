import argparse
import pathlib
from typing import Optional

import numpy as np
import pandas as pd
import decoupler as dc
import yaml
import json
import gzip


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


def main():
    params, inputs, outputs = parse_arguments()
    outputs = pathlib.Path(outputs)
    outputs.parent.mkdir(exist_ok=True, parents=True)
    data = pd.read_csv(inputs["expression"], index_col=0)
    with gzip.open(inputs["genelists"], "rb") as file:
        pathways = json.load(file)
    pathways = pathways[params["species"]][params["idtype"]]
    pathways = pd.Series(pathways)
    pathways = pathways.explode().to_frame().reset_index()
    pathways.columns = ["name", "gene"]
    pathways["weight"] = 1.0
    pathways = pathways.drop_duplicates()

    if params["omic"] == "rnaseq":
        data.index = data.index.str.split(".", n=1, expand=True).get_level_values(0)
    elif params["omic"] == "proteomics":
        data = data.T

    ncounts = pathways.name.value_counts().to_frame()
    ncounts.columns = ["before"]
    pathways = pathways[pathways.gene.isin(data.index)]
    ncounts["after"] = pathways.name.value_counts()
    ncounts["fraction"] = ncounts["after"] / ncounts["before"]
    ncounts.to_csv(outputs.with_suffix(".ncount.csv"))

    res = dc.decouple(
        mat=data.T,
        net=pathways,
        source="name",
        target="gene",
        weight="weight",
        methods=params["methods"],
        args=params["args"],
    )
    res = pd.concat(res).to_csv(outputs.with_suffix(".csv"))


if __name__ == "__main__":
    main()
