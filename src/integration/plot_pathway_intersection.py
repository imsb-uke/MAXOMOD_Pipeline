#!/usr/bin/env python3

from typing import Optional
import pandas as pd
import argparse
import yaml
import pathlib
import matplotlib.pyplot as plt
import upsetplot


def get_upset_data_frame(data, reference, query):
    sets = data.groupby(reference)[query].unique().map(set)
    all_values = data[query].unique()
    res = pd.DataFrame({value: sets.map(lambda x: value in x) for value in all_values})
    return res


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

def read_csv(*args, **kwargs):
    try:
        data = pd.read_csv(*args, **kwargs)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()
    return data


def main():
    params, inputs, output = parse_arguments(key="integration")

    files = pd.Series(inputs)
    if isinstance(params["names"], list):
        files.index = files.index.str.split("_", expand=True)
    data = files.transform(read_csv, index_col=0)
    outdir = pathlib.Path(output)
    outdir.mkdir(exist_ok=True, parents=True)
    for key, value in data.items():
        value[params["names"]] = key
    data = pd.concat(data.values)
    data.columns = data.columns.str.replace(".", "_", regex=False)
    filtered = data.query(params["filter"])
    upsetdf = get_upset_data_frame(filtered, params["names"], params["column"])
    upsetdf = upsetdf.T
    upsetdf = upsetdf.sort_index(axis=1)
    upsetplot.plot(
        upsetdf.groupby(upsetdf.columns.to_list()).size(),
        show_percentages=params["percent"],
    )
    plt.savefig(outdir.joinpath("pathway_intersection_all.png"), bbox_inches="tight")
    groups = pd.Series(pd.NA, index=upsetdf.index)
    groups[upsetdf.sum(axis=1) == upsetdf.columns.size] = "00_all"
    if isinstance(params["names"], list):
        melted = upsetdf.melt(ignore_index=False).reset_index()
        for i, key in enumerate(params["names"]):
            for name in upsetdf.columns.levels[i]:
                tmp = upsetdf.xs(name, level=key, axis=1)
                tmp = tmp[tmp.sum(axis=1) > 0]
                if tmp.columns.size < 2 or tmp.index.size < 2:
                    continue
                upsetplot.plot(
                    tmp.groupby(tmp.columns.to_list()).size(),
                    show_percentages=params["percent"],
                )
                plt.savefig(
                    outdir.joinpath(f"pathway_intersection_{key}_{name}.png"),
                    bbox_inches="tight",
                )
            tmp = melted.groupby(["index", key]).value.sum() > 0
            tmp = tmp[tmp].reset_index().drop_duplicates(keep=False, subset="index")
            groups.loc[tmp["index"]] = tmp[key].values
        tmp = melted.groupby(["index"] + params["names"]).value.all()
        tmp = tmp[tmp].reset_index().drop_duplicates(keep=False, subset="index")
        groups.loc[tmp["index"]] = (
            tmp[params["names"]].apply(lambda x: ",".join(x.tolist()), axis=1).values
        )
    else:
        tmp = melted.groupby(["index", params["names"]]).value.all()
        tmp = tmp[tmp].reset_index().drop_duplicates(keep=False, subset="index")
        groups.loc[tmp["index"]] = tmp[params["names"]].values

    groups.name = "Condition"
    groups = groups.dropna().sort_values()
    groups.to_csv(outdir.joinpath(f"common_and_unique_pathways.csv"))


if __name__ == "__main__":
    main()
