#!/usr/bin/env python


import pathlib

import pandas as pd
import tqdm
from statsmodels.stats.multitest import multipletests


def rename_columns(data, namings):
    cols = data.columns.to_series()
    cols = cols.replace(namings)
    data.columns = cols.values


def merge_other(data, columns):
    tmp = data.loc[
        :, ~((data.columns.isin(columns)) | (data.columns.str.contains("Unnamed:")))
    ]
    data = data.loc[:, data.columns.isin(columns)].copy()
    data["other"] = tmp.apply(
        lambda x: ";".join((f"{k}={v}" for k, v in x.items())), axis=1
    )
    return data


def apply_padjust(data):
    if "p.adjust" in data.columns:
        return

    pvalues = data.pvalue.dropna()
    pvalues = pd.Series(
        multipletests(pvalues.values, method="fdr_bh")[1], index=pvalues.index
    )
    data["p.adjust"] = pvalues.reindex(data.index)

def parse_das_dtu(folder, source):
    folder = pathlib.Path(folder)
    files = pd.Series(folder.rglob("ORA/ORA_*.csv"))
    files = files.to_frame()
    files.columns = ["file"]
    files["model"] = pd.NA
    models = files.file.map(lambda x: x.parts[2])
    files.loc[models == "FUS-mouse", "model"] = "FUS"
    files.loc[models == "C9orf72-mouse", "model"] = "C9orf72"
    files.loc[models == "TDP43-mouse", "model"] = "TDP43"
    files.loc[models == "SOD1-mouse", "model"] = "SOD1"
    files.loc[models == "human", "model"] = "HUMAN"
    files["omic"] = "rnaseq"
    files["source"] = source

    files["enrichment"] = "OverRepresentation"
    files["sex"] = (
        files.file.map(lambda x: x.stem)
        .str.split("_", expand=True)
        .iloc[:, 1]
        .str.lower()
        .replace({"fem": "female"})
    )
    files["type"] = "GO"
    cols = files.columns[files.columns != "file"]
    files = files.set_index(cols.to_list())
    data = files.file.map(pd.read_csv)
    data = data[~data.map(lambda x: x.empty)].copy()

    namings = {
        "ONTOLOGY": "subtype",
        "ID": "ID",
        "Pathway": "Description",
        "pvalue": "pvalue",
        "geneID": "genes",
        "Count": "Gene_count",
        "score": "score",
        "p.adjust": "p.adjust",
    }
    data.transform(rename_columns, namings=namings)
    data = data.transform(merge_other, columns=namings.values())
    for i, value in data.items():
        value[data.index.names] = i
    data = pd.concat(data.values)
    data.genes = data.genes.str.strip().str.replace("\s+", ",", regex=True)
    pathway_id = data.ID.str.extract("(GO:.+$)|(\?(mmu.+$))")
    pathway_id[0].update(pathway_id[2])
    data.ID = pathway_id[0]
    return data


def parse_das(folder):
    return parse_das_dtu(folder, "DAS")


def parse_dtu(folder):
    return parse_das_dtu(folder, "DTU")


def parse_wgcna(folder):
    folder = pathlib.Path(folder)
    files = pd.Series(folder.rglob("*.csv"))
    files = files.to_frame()
    files.columns = ["file"]
    filename = files.file.astype(str).str.lower()
    files = files[~files.file.astype(str).str.contains("_All_")]
    files = files[files.file.astype(str).str.contains("_ME")]
    files["omic"] = "proteomics"
    files.loc[filename.str.contains("rnaseq"), "omic"] = "rnaseq"
    files["model"] = pd.NA

    files.loc[filename.str.contains("fus"), "model"] = "FUS"
    files.loc[filename.str.contains("c9orf72"), "model"] = "C9orf72"
    files.loc[filename.str.contains("c9orf"), "model"] = "C9orf72"
    files.loc[filename.str.contains("tdp43"), "model"] = "TDP43"
    files.loc[filename.str.contains("sod1"), "model"] = "SOD1"
    files.loc[filename.str.contains("human"), "model"] = "HUMAN"

    files["type"] = pd.NA
    files.loc[filename.str.contains("go_enr"), "type"] = "GO"
    files.loc[filename.str.contains("kegg_enr"), "type"] = "KEGG"
    files["subtype"] = pd.NA
    files.loc[filename.str.contains("_cc_"), "subtype"] = "CC"
    files.loc[filename.str.contains("_bp_"), "subtype"] = "BP"
    files.loc[filename.str.contains("_mf_"), "subtype"] = "MF"
    files["enrichment"] = "OverRepresentation"
    files["source"] = "WGCNA"

    files["ModuleID"] = files.file.astype(str).str.extract("(ME\w+)")[0]

    files["sex"] = pd.NA
    files.loc[filename.str.contains("_males"), "sex"] = "male"
    files.loc[filename.str.contains("_females"), "sex"] = "female"

    cols = files.columns[files.columns != "file"]
    files = files.set_index(cols.to_list())
    data = files.file.map(pd.read_csv)
    data = data[~data.map(lambda x: x.empty)]
    namings = {
        "gs_id": "ID",
        "gs_description": "Description",
        "gs_pvalue": "pvalue",
        "gs_genes": "genes",
        "z_score": "score",
        "geneID": "genes",
        "gs_de_count": "Gene_count",
        "Count": "Gene_count",
        "p.adjust": "p.adjust",
    }
    data.transform(rename_columns, namings=namings)
    data = data.transform(merge_other, columns=namings.values())
    for i, value in data.items():
        value[data.index.names] = i
    return pd.concat(data.values)


def parse_rnaseq(*args, **kwargs):
    files = pd.Series(
        pathlib.Path("datasets/consortium/").rglob(
            "rnaseq/gene-enrichment/gene_set_enrichment/only_*males/*.csv"
        )
    )
    files = files.to_frame()
    files.columns = ["file"]
    files["omic"] = "rnaseq"
    files["sex"] = (
        files.file.map(lambda x: x.parent.name)
        .str.replace("only_", "")
        .str.replace("males", "male")
    )
    files["model"] = (
        files.file.map(lambda x: x.parent.parent.parent.parent.parent.parent.name)
        .str.replace("-mouse", "")
        .str.replace("-datasets", "")
        .replace({"human": "HUMAN"})
    )
    files["enrichment"] = "GSEA"
    files["type"] = "GO"
    files["source"] = "DEG"
    files["subtype"] = (
        files.file.map(lambda x: x.stem).str.rsplit("-", n=1, expand=True).iloc[:, -1]
    )
    cols = files.columns[files.columns != "file"]
    files = files.set_index(cols.to_list())
    data = files.file.map(pd.read_csv)
    data = data[~data.map(lambda x: x.empty)]
    namings = {
        "gs_id": "ID",
        "gs_description": "Description",
        "gs_pvalue": "pvalue",
        "core_enrichment": "genes",
        "NES": "score",
        "geneID": "genes",
        "gs_de_count": "Gene_count",
        "setSize": "Gene_count",
        "p.adjust": "p.adjust",
    }
    data.transform(rename_columns, namings=namings)
    data = data.transform(merge_other, columns=namings.values())
    for i, value in data.items():
        value[data.index.names] = i
    return pd.concat(data.values)


def parse_proteomics(*args, **kwargs):
    files = pd.Series(
        pathlib.Path("datasets/consortium/").rglob(
            "proteomics/gsea_enrichment/gsea_results/only_*males-*.csv"
        )
    )
    files = files.to_frame()
    files.columns = ["file"]
    files["omic"] = "proteomics"
    files["sex"] = (
        files.file.map(lambda x: x.stem)
        .str.replace("only_", "")
        .str.split("-", n=1, expand=True)
        .iloc[:, 0]
        .str.replace("males", "male")
    )
    files["model"] = (
        files.file.map(lambda x: x.parent.parent.parent.parent.parent.name)
        .str.replace("-mouse", "")
        .str.replace("-datasets", "")
        .replace({"human": "HUMAN"})
    )
    files["enrichment"] = "GSEA"

    files["source"] = "DEG"
    files["subtype"] = (
        files.file.map(lambda x: x.stem)
        .str.rsplit("-", n=1, expand=True)
        .iloc[:, -1]
        .str.replace("_gsea_results", "")
    )
    files["type"] = "GO"
    files.loc[files["subtype"] == "KEGG", "type"] = "KEGG"
    files["subtype"] = files["subtype"].replace({"KEGG": pd.NA})
    cols = files.columns[files.columns != "file"]
    files = files.set_index(cols.to_list())
    data = files.file.map(pd.read_csv)
    data = data[~data.map(lambda x: x.empty)]
    namings = {
        "gs_id": "ID",
        "gs_description": "Description",
        "gs_pvalue": "pvalue",
        "core_enrichment": "genes",
        "NES": "score",
        "geneID": "genes",
        "gs_de_count": "Gene_count",
        "setSize": "Gene_count",
        "p.adjust": "p.adjust",
    }
    data.transform(rename_columns, namings=namings)
    data = data.transform(merge_other, columns=namings.values())
    for i, value in data.items():
        value[data.index.names] = i
    return pd.concat(data.values)


def parse_dcloc(*args, **kwargs):
    files = pd.Series(
        pathlib.Path("datasets/consortium/").rglob("networks/dcloc/enrichment/**/*.csv")
    )
    files = files.to_frame()
    files.columns = ["file"]
    files["omic"] = files.file.map(
        lambda x: x.parent.parent.parent.parent.parent.parent.name
    )
    files["model"] = files.file.map(
        lambda x: x.parent.parent.parent.parent.parent.parent.parent.parent.name
    )
    files["sex"] = files.file.map(lambda x: x.parent.parent.name)
    files["model"] = (
        files["model"]
        .str.replace("-mouse", "")
        .str.replace("-datasets", "")
        .replace({"human": "HUMAN"})
    )
    files["type"] = files.file.map(lambda x: x.parent.name).str.upper()
    files["enrichment"] = pd.NA
    files.loc[files["type"] == "KEGG", "enrichment"] = "OverRepresentation"
    files.loc[files["type"] == "GSEA", "enrichment"] = "GSEA"
    files.loc[files["type"] == "GO", "enrichment"] = "OverRepresentation"
    files.loc[files["type"] == "GSEA", "type"] = "KEGG"
    files["subtype"] = (
        files.file.map(lambda x: x.stem)
        .str.replace("pathways", "")
        .str.replace("_", "")
        .replace({"": pd.NA})
    )
    files["source"] = "DCLOC"
    cols = files.columns[files.columns != "file"]
    files = files.set_index(cols.to_list())
    files = files[files.file.map(lambda x: x.stat().st_size) > 0]
    data = files.file.map(pd.read_csv)
    data = data[~data.map(lambda x: x.empty)]
    namings = {
        "gs_id": "ID",
        "gs_description": "Description",
        "gs_pvalue": "pvalue",
        "core_enrichment": "genes",
        "NES": "score",
        "geneID": "genes",
        "gs_de_count": "Gene_count",
        "Count": "Gene_count",
        "setSize": "Gene_count",
        "p.adjust": "p.adjust",
    }
    data.transform(rename_columns, namings=namings)
    data = data.transform(merge_other, columns=namings.values())
    for i, value in data.items():
        value[data.index.names] = i
    return pd.concat(data.values)


def parse_wgcna_proteomics(*args,**kwargs):
    return parse_wgcna(*args,**kwargs)

def parse_wgcna_rnaseq(*args,**kwargs):
    return parse_wgcna(*args,**kwargs)

def parse_genetonic(folder):
    files = pd.Series(
        pathlib.Path("datasets/consortium/").rglob("rnaseq/genetonic/*_de_topGO_*.csv")
    )
    files = files.to_frame()
    files.columns = ["file"]
    files.file = files.file.astype(str)
    files["model"] = pd.NA
    files["source"] = "GeneTonic"
    files.loc[files.file.str.contains("FUS"), "model"] = "FUS"
    files.loc[files.file.str.contains("C9orf72"), "model"] = "C9orf72"
    files.loc[files.file.str.contains("TDP43"), "model"] = "TDP43"
    files.loc[files.file.str.contains("SOD1"), "model"] = "SOD1"
    files.loc[files.file.str.contains("human"), "model"] = "HUMAN"
    files["omic"] = "rnaseq"
    files["type"] = "GO"
    files["subtype"] = pd.NA
    files.loc[files.file.str.contains("topGO_BP"), "subtype"] = "BP"
    files.loc[files.file.str.contains("topGO_MF"), "subtype"] = "MF"
    files.loc[files.file.str.contains("topGO_CC"), "subtype"] = "CC"
    files["enrichment"] = "OverRepresentation"
    files["sex"] = pd.NA
    files.loc[files.file.str.contains("male"), "sex"] = "male"
    files.loc[files.file.str.contains("female"), "sex"] = "female"
    cols = files.columns[files.columns != "file"]
    files = files.set_index(cols.to_list())
    data = files.file.map(pd.read_csv)
    data = data[~data.map(lambda x: x.empty)]

    namings = {
        "gs_id": "ID",
        "gs_description": "Description",
        "gs_pvalue": "p.adjust",
        "gs_genes": "genes",
        "z_score": "score",
        "gs_de_count": "Gene_count",
    }
    data.transform(rename_columns, namings=namings)
    data.map(apply_padjust)
    data = data.transform(merge_other, columns=namings.values())
    for i, value in data.items():
        value[data.index.names] = i
    return pd.concat(data.values)



def main():

    functions = {
        "rnaseq": None,
        "proteomics": None,
        "wgcna_proteomics": "results/WGCNA",
        "wgcna_rnaseq": "results/WGCNA_rnaseq",
        "dcloc": None,
        "dtu": "results/DTU/",
        "das": "results/DAS",
        "genetonic": None,
    }
    data = [
        globals()[f"parse_{method}"](arg)
        for method, arg in tqdm.tqdm(functions.items())
    ]
    pd.concat(data).to_csv("results/integration/pathways/all_pathways.csv", index=False)


if __name__ == "__main__":
    main()
