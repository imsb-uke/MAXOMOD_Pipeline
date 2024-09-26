import argparse
import csv
import pathlib

import pandas as pd
import tqdm


def get_introns(data):
    gene = data.query("Type == 'gene'").iloc[0]
    start = gene.Start
    results = []
    for _, exon in (
        data.query("Type != 'gene'").sort_values("Start", ascending=True).iterrows()
    ):
        results.append({"Start": start, "End": exon.Start})
        start = exon.End + 1
    results.append({"Start": start, "End": gene.End})
    results = pd.DataFrame(results)
    results["Type"] = "intron"
    for col in gene.index:
        if col not in results.columns:
            results[col] = gene[col]
    return results


_COL_ORDER = [
    "Chrom",
    "Source",
    "Type",
    "Start",
    "End",
    ".1",
    "Strand",
    ".2",
    "other",
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=str)
    parser.add_argument("output", type=str)
    tqdm.tqdm.pandas()
    args = parser.parse_args()
    outputdir = pathlib.Path(args.output)
    outputdir.mkdir(exist_ok=True, parents=True)
    data = pd.read_csv(
        args.input,
        sep="\t",
        comment="#",
        header=None,
    )
    data.columns = _COL_ORDER
    data.other = data.other.str.split(";", n=1, expand=True).iloc[:, 0] + ";"
    data = data[data.Type.isin(("gene", "exon"))]
    genes = data.groupby("other").Type.value_counts().unstack().query("gene == 1")
    data = data[data.other.isin(genes.index)]
    introns = data.groupby("other").progress_apply(get_introns).reset_index(drop=True)
    introns = introns[introns.End > introns.Start]
    data[_COL_ORDER].to_csv(
        outputdir.joinpath("exons.gtf"),
        sep="\t",
        header=False,
        index=False,
        quoting=csv.QUOTE_NONE,
        doublequote=False,
    )
    introns[_COL_ORDER].to_csv(
        outputdir.joinpath("introns.gtf"),
        sep="\t",
        header=False,
        index=False,
        quoting=csv.QUOTE_NONE,
        doublequote=False,
    )


if __name__ == "__main__":
    main()
