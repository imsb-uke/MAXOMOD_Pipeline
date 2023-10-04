import argparse
import pathlib
from typing import Optional

import numpy as np
import pandas as pd
import yaml
import networkx as nx
import itertools


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
    params, inputs, outputs = parse_arguments(key="WGCNA", method="graph_filter")
    outputs = pathlib.Path(outputs)
    graph = nx.read_gml(inputs)
    if params["min_size"] is not None:
        components = [
            c for c in nx.connected_components(graph) if len(c) < params["min_size"]
        ]
        graph.remove_nodes_from(itertools.chain.from_iterable(components))
    if params["max_size"] is not None:
        components = [
            c for c in nx.connected_components(graph) if len(c) > params["max_size"]
        ]
        graph.remove_nodes_from(itertools.chain.from_iterable(components))
    layout = getattr(nx, f"{params['layout']}_layout")
    kwargs = params["layout_kwargs"]
    if not kwargs:
        kwargs = {}
    pos = layout(graph, weight="weight", **kwargs)
    np.savez(outputs.joinpath("layout.npz"), **pos)
    nx.write_gml(graph, outputs.joinpath("filtered.gml"))


if __name__ == "__main__":
    main()
