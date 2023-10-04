import pandas as pd
import pathlib
import yaml
import json
import gzip


def read(filename):
    try:
        data = pd.read_csv(filename, usecols=["GENE NAME"]).squeeze()
    except ValueError:
        data = pd.read_csv(filename, usecols=[0], header=None).squeeze()
    return data


def main():
    files = pd.Series(pathlib.Path("database/genesets/raw_lists/").glob("*.csv"))
    files.index = (
        files.map(lambda x: x.stem)
        .str.extract("List of GO Genes -\s*([^(\.]+)")[0]
        .str.strip()
    )
    data = files.transform(read).stack().reset_index(level=-1, drop=True)
    data.name = "SYMBOL"
    data = data.to_frame()
    mapping = pd.read_csv(
        "database/uniprot/uniprot_HUMAN_9606_idmapping_selected.tab.gz",
        sep="\t",
        header=None,
        usecols=[1, 18],
        names=["UNIPROT", "ENSEMBL"],
    )
    mapping = (
        mapping.dropna()
        .set_index("UNIPROT")
        .ENSEMBL.str.split(";")
        .explode()
        .str.strip()
        .reset_index()
    )
    mapping = pd.merge(
        mapping,
        pd.read_csv("database/biomaRt/ensembl_to_symbol_human.csv"),
        left_on="ENSEMBL",
        right_on="ensembl_gene_id",
        how="outer",
    ).drop(columns=["ENSEMBL"])
    data = data.reset_index().rename({0: "pathway"}, axis=1)
    data = pd.merge(data, mapping, left_on="SYMBOL", right_on="hgnc_symbol").drop(
        columns=["SYMBOL", "hgnc_symbol"]
    )
    mapping = pd.read_csv(
        "human_rat_mouse_protein_gene_link.csv.gz",
        usecols=["Gene stable ID (human)", "Gene stable ID (mouse)"],
        sep="\t",
    ).dropna()
    mapping = mapping.drop_duplicates()
    mapping2 = pd.read_csv(
        "database/uniprot/uniprot_MOUSE_10090_idmapping_selected.tab.gz",
        sep="\t",
        header=None,
        usecols=[1, 18],
        names=["UNIPROT", "ENSEMBL"],
    )
    mapping2 = mapping2.dropna(subset=["ENSEMBL"])
    mapping = pd.merge(
        mapping,
        mapping2,
        left_on="Gene stable ID (mouse)",
        right_on="ENSEMBL",
        how="left",
    ).drop(columns=["ENSEMBL"])
    data = pd.merge(
        data,
        mapping,
        left_on=["ensembl_gene_id"],
        right_on=["Gene stable ID (human)"],
        suffixes=("_human", "_mouse"),
        how="left",
    )
    data.duplicated().value_counts()
    data = data.drop_duplicates()
    data = data.drop(columns=["Gene stable ID (human)"]).rename(
        columns={
            "ensembl_gene_id": "ENSEMBL_human",
            "Gene stable ID (mouse)": "ENSEMBL_mouse",
        }
    )
    data = (
        data.groupby("pathway")
        .apply(lambda x: x.apply(lambda y: set(y.dropna())))
        .drop(columns="pathway")
        .stack()
        .reset_index()
    )
    data["genes"] = data[0].map(list)
    data[["IDtype", "Species"]] = data.level_1.str.split("_", expand=True)
    data = data.drop(columns=["level_1", 0])
    results = {}
    for species, g1 in data.groupby("Species"):
        results[species] = {}
        for ids, g2 in g1.groupby("IDtype"):
            results[species][ids] = g2.set_index("pathway").genes.to_dict()

    with gzip.open("database/genesets/all_genesets.json.gz", "wt") as file:
        json.dump(results, file)


if __name__ == "__main__":
    main()
