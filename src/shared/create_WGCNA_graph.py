import argparse
import pathlib
from typing import Optional

import numpy as np
import pandas as pd
import yaml
import networkx as nx


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
    params, inputs, outputs = parse_arguments(key="WGCNA", method="graph")
    outputs = pathlib.Path(outputs)
    outputs.parent.mkdir(exist_ok=True, parents=True)
    network = pd.read_feather(inputs)
    network["source"] = network.columns
    network = network.melt(id_vars="source", var_name="target", value_name="weight")
    network = network.query(f"weight > {params['threshold']}")
    network = network[network.target != network.source]
    network = nx.from_pandas_edgelist(network, edge_attr="weight")
    nx.write_gml(network, outputs)


if __name__ == "__main__":
    main()
