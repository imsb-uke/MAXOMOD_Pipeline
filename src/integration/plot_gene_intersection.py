#!/usr/bin/env python3

from typing import Optional
import pandas as pd
import argparse
import yaml
import pathlib
import matplotlib.pyplot as plt
import upsetplot
import tqdm


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


def read_id_mapping(filename):
    idmapping = pd.read_csv(filename, sep="\t", index_col=1, header=None)
    idmapping = idmapping[18].dropna().str.split(";", expand=True)[0].dropna()
    idmapping.name = "ENSEMBL"
    return idmapping


def mapids(data, idmapping):
    mapping = idmapping.loc[:, data.name].dropna()
    mapping.name = "ENSEMBL"
    data = pd.merge(
        data["data"], mapping, left_index=True, right_index=True, how="left"
    )
    data.set_index("ENSEMBL", inplace=True)
    return data


def map_species(data, speciesmapping):
    cols = data["data"].columns.to_list()
    tmp = data["data"].loc[~data["data"].index.isna()]
    data = pd.merge(
        tmp,
        speciesmapping,
        left_index=True,
        right_on=f"Gene stable ID ({data['type']})",
        how="left",
    )
    data = (
        data.loc[
            :, (data.columns.isin(cols) | (data.columns == "Gene stable ID (human)"))
        ]
        .drop_duplicates(subset=["Gene stable ID (human)"])
        .dropna()
    )
    return data


def read_csv(*args, **kwargs):
    try:
        data = pd.read_csv(*args, **kwargs)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()
    return data


def remove_ensembl_suffix(data):
    data.index = data.index.str.replace("\.\d+$", "", regex=True)


def filtering(data, filtervalue):
    data.columns = data.columns.str.replace(".", "_", regex=False)
    filtered = data.query(filtervalue)
    return filtered


def main():
    params, inputs, output = parse_arguments(key="integration")

    data = pd.DataFrame.from_dict(params["metadata"], orient="columns")
    data["data"] = pd.Series(inputs).transform(read_csv, index_col=0)
    outdir = pathlib.Path(output)
    outdir.mkdir(exist_ok=True, parents=True)
    tqdm.tqdm.pandas(desc="Filter")
    data["data"] = data["data"].progress_apply(filtering, filtervalue=params["filter"])

    if params["idtype"] == "protein":
        tqdm.tqdm.pandas(desc="Map Proteins")
        idmapping = pd.Series(params["mapping"]).transform(read_id_mapping).T
        data["data"] = data.progress_apply(mapids, idmapping=idmapping, axis=1)
    else:
        tqdm.tqdm.pandas(desc="Remove Ensembl version")
        data["data"].progress_apply(remove_ensembl_suffix)

    tqdm.tqdm.pandas(desc="Map Species")

    species_mapping = pd.read_csv(params["speciesmapping"], sep="\t", index_col=0)
    data["data"] = data.progress_apply(
        map_species, speciesmapping=species_mapping, axis=1
    )
    data.set_index(params["names"], inplace=True)
    for key, value in data["data"].items():
        value[params["names"]] = key

    data = pd.concat(data["data"].values)

    upsetdf = get_upset_data_frame(data, params["names"], "Gene stable ID (human)")
    upsetdf = upsetdf.T
    upsetdf = upsetdf.sort_index(axis=1)
    upsetplot.plot(
        upsetdf.groupby(upsetdf.columns.to_list()).size(),
        show_percentages=params["percent"],
    )
    plt.savefig(outdir.joinpath("gene_intersection_all.png"), bbox_inches="tight")
    groups = pd.Series(pd.NA, index=upsetdf.index)
    groups[upsetdf.sum(axis=1) == upsetdf.columns.size] = "00_all"
    if isinstance(params["names"], list):
        melted = upsetdf.melt(ignore_index=False).reset_index()
        for i, key in enumerate(params["names"]):
            for name in upsetdf.columns.levels[i]:
                tmp = upsetdf.xs(name, level=key, axis=1)
                tmp = tmp[tmp.sum(axis=1) > 0]
                try:
                    upsetplot.plot(
                        tmp.groupby(tmp.columns.to_list()).size(),
                        show_percentages=params["percent"],
                    )
                except AttributeError:
                    continue
                else:
                    plt.savefig(
                        outdir.joinpath(f"gene_intersection_{key}_{name}.png"),
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
    groups.to_csv(outdir.joinpath(f"common_and_unique_genes.csv"))


if __name__ == "__main__":
    main()
