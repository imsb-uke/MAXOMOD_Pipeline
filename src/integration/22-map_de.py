#!/usr/bin/env python3

import argparse
import pathlib
from typing import Optional

import numpy as np
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


def map_proteomics(distances, inputs, params):
    mapping = pd.read_csv(
        inputs["mapping"],
        sep="\t",
        header=None,
        usecols=[1, 18],
        names=["Protein", "_id"],
    ).dropna()
    mapping.set_index("Protein", inplace=True)
    mapping._id = mapping._id.str.split(";", expand=True).iloc[:, 0].str.strip()
    distances = pd.merge(
        distances, mapping, left_index=True, right_index=True, how="left"
    )
    distances = map_ensembl_by_species(distances, inputs, params)
    return distances.drop_duplicates()


def map_rnaseq(distances, inputs, params):
    distances["_id"] = (
        distances.index.to_series().str.split(".", expand=True).iloc[:, 0]
    )
    distances = map_ensembl_by_species(distances, inputs, params)
    return distances.drop_duplicates()


def map_ensembl_by_species(distances, inputs, params):
    species_mapping = (
        pd.read_csv(
            params["species_mapping"],
            sep="\t",
            usecols=[
                f"Gene name ({params['target']})",
                f"Gene stable ID ({params['species']})",
                f"Gene stable ID ({params['target']})",
                f"NCBI gene (formerly Entrezgene) ID ({params['target']})",
            ],
        )
        .dropna(how="all")
        .dropna(subset=[f"Gene stable ID ({params['species']})"])
    )
    species_mapping["_id"] = species_mapping[f"Gene stable ID ({params['species']})"]
    species_mapping = species_mapping.rename(
        columns={
            f"Gene name ({params['target']})": "SYMBOL",
            f"Gene stable ID ({params['target']})": "ENSEMBL",
            f"NCBI gene (formerly Entrezgene) ID ({params['target']})": "ENTREZ",
        }
    )
    species_mapping["ENTREZ"] = pd.Categorical(species_mapping["ENTREZ"])
    species_mapping["ENTREZ"].cat.categories = species_mapping[
        "ENTREZ"
    ].cat.categories.astype(int)

    distances = pd.merge(
        distances, species_mapping, left_on="_id", right_on="_id", how="left"
    )
    distances = distances.drop(
        columns=f"Gene stable ID ({params['species']})", errors="ignore"
    )
    return distances.drop_duplicates()


def main():
    params, inputs, outputs = parse_arguments()
    outputs = pathlib.Path(outputs)
    outputs.parent.mkdir(exist_ok=True, parents=True)
    data = pd.read_csv(inputs["de"], sep=params.get("sep", ","), index_col=0)
    mapping = globals()[f"map_{params['datatype']}"]
    data = mapping(data, inputs, params).drop(columns="_id", errors="ignore")
    data.to_csv(outputs)


if __name__ == "__main__":
    main()
