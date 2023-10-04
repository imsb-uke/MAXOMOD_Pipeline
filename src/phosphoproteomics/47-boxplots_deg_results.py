#!/usr/bin/env python3

import argparse
import pathlib
from typing import Optional
import yaml

import pandas as pd
import seaborn as sns
from matplotlib.patches import Rectangle


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
        params["stages"][args.stage].get("output", None),
    )


def subset_data(vals, data, meta=None):
    data = data[vals.index].melt(
        ignore_index=False, var_name="PhosphoSite", value_name="Abundance"
    )
    if meta is not None:
        data = pd.merge(data, meta, left_index=True, right_index=True)
    return data


def main():
    settings, inputs, outputs = parse_arguments(
        key="phosphoproteomics", method="deg_boxplots"
    )
    outdir = pathlib.Path(outputs)
    outdir.mkdir(exist_ok=True, parents=True)

    raw = pd.read_csv(inputs["raw"], index_col=0)

    annot = pd.read_csv(inputs["metadata"], index_col=0)

    data = pd.Series(pathlib.Path(inputs["degpath"]).glob(settings["filepattern"]))

    data.index = data.map(lambda x: x.stem)

    data = data.to_frame()

    data.columns = ["filename"]

    data[settings["splitby"]] = data.index.str.extract(settings["pattern"])[0].values

    data["data"] = data.filename.transform(pd.read_csv, index_col=0)

    data["raw"] = data["data"].transform(subset_data, data=raw, meta=annot)
    col_wrap = settings["catplot_kwargs"].pop("col_wrap", 4)
    for name, row in data.iterrows():
        if row.raw.empty:
            continue
        g = sns.catplot(
            x=settings["splitby"],
            y="Abundance",
            col="PhosphoSite",
            data=row.raw,
            **settings["catplot_kwargs"],
            col_wrap=min(col_wrap, row.raw.PhosphoSite.nunique()),
        )
        g.map_dataframe(
            sns.stripplot,
            x=settings["splitby"],
            y="Abundance",
            **settings["stripplot_kwargs"],
        )
        ax = g.axes.flatten()[-1]
        ticks = ax.get_xticklabels()
        diff = (ticks[1].get_position()[0] - ticks[0].get_position()[0]) / 2
        mapping = {t.get_text(): t.get_position()[0] for t in ticks}
        for site, ax in g.axes_dict.items():
            ylim = ax.get_ylim()
            xlim = ax.get_xlim()
            ax.add_patch(
                Rectangle(
                    (mapping[row[settings["splitby"]]] - diff, ylim[0]),
                    2 * diff,
                    ylim[1] - ylim[0],
                    fc="none",
                    color="red",
                    lw=2,
                )
            )
            ax.set_title(
                f"{site} \n(l2fc:  {row.data.loc[site, 'log2FoldChange']:.2f}, padj: {row.data.loc[site, 'padj']:.1e})"
            )
        g.figure.savefig(
            outdir.joinpath(name).with_suffix(".png"), dpi=200, bbox_inches="tight"
        )


if __name__ == "__main__":
    main()
