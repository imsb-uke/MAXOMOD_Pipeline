import argparse
import pathlib
from typing import Optional

import numpy as np
import pandas as pd
import yaml
import networkx as nx
import json
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns


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
    parser.add_argument(
        "--omic",
        required=True,
        type=str,
        help="omic key from params.yaml to use",
    )
    args = parser.parse_args()
    params = yaml.safe_load(pathlib.Path("params.yaml").read_text())[args.key][
        args.method
    ]
    settings = dict()
    settings.update(params.get("default_system_settings"))
    settings.update(params.get("default_settings"))
    params = params[args.omic]
    settings.update(params["stages"][args.stage].get("settings", {}))
    return (
        settings,
        params["stages"][args.stage].get("inputs", {}),
        params["stages"][args.stage].get("outputs", None),
    )


def map_proteomics(genes, uniprotmapping, ensemblmapping):
    if genes.empty:
        genes["gene_name"] = None
        return genes
    mapping = pd.read_csv(uniprotmapping, sep="\t", usecols=[1, 18], header=None)
    mapping = mapping.set_index(1)[18].dropna().str.split(";").explode().str.strip()
    mapping.name = "ENSEMBL"
    genes = pd.merge(genes, mapping, left_on="UniProtName", right_index=True)
    mapping = pd.read_csv(ensemblmapping, index_col="ensembl_gene_id").rename(
        columns={"hgnc_symbol": "gene_name", "mgi_symbol": "gene_name"}
    )
    if genes.empty:
        genes["ENSEMBL"] = None
    genes = pd.merge(genes, mapping, left_on="ENSEMBL", right_index=True)
    if genes.empty:
        genes["gene_name"] = None
    return genes


def map_regulation(genes):
    if genes.empty:
        genes["model"] = None
    regulation = (
        genes.pivot_table(
            index="gene_name",
            columns="model",
            values="regulation",
            fill_value="",
            aggfunc="sum",
        )
        .reindex(["male", "female"], axis=1)
        .fillna("")
    )
    regulation["regulation"] = pd.NA
    regulation.loc[regulation.male == "upregulated", "regulation"] = "male-upregulated"
    regulation.loc[
        regulation.male == "downregulated", "regulation"
    ] = "male-downregulated"
    regulation.loc[
        regulation.female == "upregulated", "regulation"
    ] = "female-upregulated"
    regulation.loc[
        regulation.female == "downregulated", "regulation"
    ] = "female-downregulated"
    regulation.loc[
        (regulation.male == "upregulated") & (regulation.female == "upregulated"),
        "regulation",
    ] = "common-upregulated"
    regulation.loc[
        (regulation.male == "downregulated") & (regulation.female == "downregulated"),
        "regulation",
    ] = "common-downregulated"
    regulation.loc[
        (regulation.male == "upregulated") & (regulation.female == "downregulated"),
        "regulation",
    ] = "female-downregulated,male-upregulated"
    regulation.loc[
        (regulation.female == "upregulated") & (regulation.male == "downregulated"),
        "regulation",
    ] = "female-upregulated,male-downregulated"
    merged = pd.merge(
        genes.drop(columns=["regulation", "model"]),
        regulation.drop(columns=["male", "female"]),
        left_on="gene_name",
        right_index=True,
    )
    return merged.drop_duplicates(subset=["gene_name", "regulation"]).rename(
        columns={"regulation": "model"}
    )


def main():
    params, inputs, outputs = parse_arguments(key="WGCNA", method="graph_plot")
    outputs = pathlib.Path(outputs)

    graph = nx.read_gml(inputs["graph"])
    pos = dict(np.load(inputs["layout"]))
    genes = json.loads(pathlib.Path(inputs["modules"]).read_bytes())
    genes = pd.Series(genes).explode().reset_index()
    genes = genes.rename(columns={"index": "Module", 0: "name"})
    genes["ColorName"] = genes.Module.str.slice(
        2,
    )
    genes["color"] = genes["ColorName"].map(mcolors.CSS4_COLORS)
    colors = genes.set_index("name").color.reindex(graph.nodes)
    colormapping = (
        genes[["Module", "color"]].drop_duplicates().set_index("Module").color.dropna()
    )
    handles = [Patch(facecolor=value) for name, value in colormapping.items()]
    plt.figure(figsize=(12, 12))
    nx.draw(
        graph,
        pos=pos,
        node_color=colors.fillna("none"),
        node_size=10,
        edge_color="gray",
    )
    plt.legend(
        handles,
        colormapping.to_dict(),
        bbox_to_anchor=(1.1, 1.1),
        bbox_transform=plt.gcf().transFigure,
        loc="upper right",
    )
    plt.savefig(outputs.joinpath("graph_module.pdf"), dpi=300, bbox_inches="tight")
    plt.close()

    signature = pd.read_csv(
        inputs["celltype_marker"], usecols=["celltype", params["species"]]
    )
    signature = signature.rename(columns={params["species"]: "name"})
    signature = signature[~signature.name.duplicated(keep=False)]
    mapping = {
        k: sns.color_palette("tab20")[i]
        for i, k in enumerate(signature.celltype.unique())
    }
    signature["color"] = signature.celltype.map(mapping)
    colors = signature.set_index("name").color.reindex(graph.nodes)
    colormapping = (
        signature[["celltype", "color"]]
        .drop_duplicates()
        .set_index("celltype")
        .color.dropna()
    )
    handles = [Patch(facecolor=value) for name, value in colormapping.items()]
    plt.figure(figsize=(12, 12))
    nx.draw(
        graph,
        pos=pos,
        node_color=colors.fillna("none"),
        node_size=10,
        edge_color="gray",
    )
    plt.legend(
        handles,
        colormapping.to_dict(),
        bbox_to_anchor=(1.1, 1.1),
        bbox_transform=plt.gcf().transFigure,
        loc="upper right",
    )
    plt.savefig(outputs.joinpath("graph_celltype.pdf"), dpi=300, bbox_inches="tight")
    plt.close()

    genes_female = pd.read_csv(inputs["deg_female"])
    genes_male = pd.read_csv(inputs["deg_male"])
    genes_female["model"] = "female"
    genes_male["model"] = "male"
    genes = pd.concat((genes_male, genes_female))
    if "gene_name" not in genes.columns:
        genes = map_proteomics(
            genes, inputs["uniprotmapping"], inputs["ensemblmapping"]
        )
    genes = map_regulation(genes)

    mapping = {k: sns.color_palette()[i] for i, k in enumerate(genes.model.unique())}
    genes["color"] = genes.model.map(mapping)
    genes = genes.rename(columns={"gene_name": "name"})
    genes = genes.drop_duplicates()
    colors = genes.set_index("name").color.reindex(graph.nodes)
    colormapping = (
        genes[["model", "color"]].drop_duplicates().set_index("model").color.dropna()
    )
    handles = [Patch(facecolor=value) for name, value in colormapping.items()]
    plt.figure(figsize=(12, 12))
    nx.draw(
        graph,
        pos=pos,
        node_color=colors.fillna("none"),
        node_size=10,
        edge_color="gray",
    )
    plt.legend(
        handles,
        colormapping.to_dict(),
        bbox_to_anchor=(1.1, 1.1),
        bbox_transform=plt.gcf().transFigure,
        loc="upper right",
    )
    plt.savefig(outputs.joinpath("graph_degs.pdf"), dpi=300, bbox_inches="tight")
    plt.close()

    genes_dict = {}
    for sex in ("male", "female"):
        genes = pd.read_csv(inputs[f"deg_full_{sex}"])
        if "gene_name" not in genes.columns:
            genes = map_proteomics(
                genes, inputs["uniprotmapping"], inputs["ensemblmapping"]
            )
        genes = genes.set_index("gene_name").log2FoldChange
        genes = genes[~genes.index.duplicated()]
        genes_dict[sex] = genes.reindex(graph.nodes).fillna(0)

    vmin = min(genes_dict["male"].min(), genes_dict["female"].min())
    vmax = max(genes_dict["male"].max(), genes_dict["female"].max())
    vcenter = max(abs(vmin), abs(vmax))

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(16, 12))
    fig.set_facecolor("w")
    axes = axes.flatten()
    for i, sex in enumerate(("male", "female")):
        nx.draw(
            graph,
            pos,
            ax=axes[i],
            vmin=-vcenter,
            vmax=vcenter,
            node_color=genes_dict[sex],
            cmap="seismic",
            node_size=5,
            edge_color="gray",
        )
        axes[i].set_title(sex.capitalize())
    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    sm = plt.cm.ScalarMappable(
        cmap="seismic", norm=plt.Normalize(vmin=-vcenter, vmax=vcenter)
    )
    sm._A = []
    plt.colorbar(sm, cax=cbar_ax)
    plt.savefig(
        outputs.joinpath(f"graph_deg_foldchange.pdf"), dpi=300, bbox_inches="tight"
    )
    plt.close()


if __name__ == "__main__":
    main()
