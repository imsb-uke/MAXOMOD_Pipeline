import pandas as pd
import numpy as np
import pathlib
import yaml
import argparse
import decoupler as dc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--key", type=str, required=True)
    parser.add_argument("--method", type=str, required=True)
    args = parser.parse_args()

    with open("params.yaml") as file:
        settings = yaml.safe_load(file)["transcriptionsfactors"][args.method][args.key]

    outpath = pathlib.Path(settings["outputpath"])
    outpath.parent.mkdir(exist_ok=True, parents=True)

    net = pd.read_feather(settings["database"])
    net = net.query(settings["filter"])

    data = pd.read_csv(settings["inputs"], index_col=0)
    if "proteinmapping" in settings:
        data = data.T
        mapping = (
            pd.read_csv(
                settings["proteinmapping"],
                sep="\t",
                header=None,
                usecols=[1, 18],
                index_col=0,
            )
            .squeeze()
            .dropna()
            .str.split(";")
            .explode()
            .str.strip()
            .dropna()
        )
        mapping = mapping[~mapping.index.duplicated()]
        idx = data.index.intersection(mapping.index)
        data = data.loc[idx]
        data.index = mapping[idx].values

    mapping = pd.read_csv(settings["mapping"], index_col=0).squeeze()
    idx = data.index.intersection(mapping.index)
    data = data.loc[idx]
    data.index = mapping[idx].values
    data = data.groupby(level=0).mean()

    results = dc.decouple(
        mat=data.T,
        net=net,
        source="source",
        target="target",
        weight="weight",
        methods=settings["methods"],
        **settings["kwargs"],
    )

    data = pd.concat(results, axis=1, ignore_index=False)
    data = data.T.reset_index()
    data.to_feather(outpath)


if __name__ == "__main__":
    main()
