import pandas as pd
import pathlib
import yaml
import argparse
import tempfile
import tqdm.auto as tqdm
import subprocess
from typing import Optional


def run_gofigure(data, commandlineargs, outputdir):
    if len(data) == 0:
        return None
    name = data.name
    if not isinstance(name, str):
        name = "_".join(map(str, name))
    data = data.drop_duplicates(subset=["ID"], keep=False)
    with tempfile.NamedTemporaryFile() as tmpfile:
        data.to_csv(tmpfile.name, sep="\t", header=False, index=False)
        defaultsargs = [
            "--input",
            tmpfile.name,
            "--input_type",
            "standard",
            "--output",
            str(outputdir.joinpath(name)),
        ]
        commandline = (
            ["python3", "src/shared/GoFigure/gofigure.py"]
            + defaultsargs
            + commandlineargs
        )
        return subprocess.run(commandline, shell=False, check=True)


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
    settings.update(params.get("default_system_settings", {}))
    settings.update(params.get("default_settings", {}))
    settings.update(params["stages"][args.stage].get("settings", {}))
    return (
        settings,
        params["stages"][args.stage].get("inputs", {}),
        params["stages"][args.stage].get("outputs", None),
    )


def main():
    params, inputs, outputs = parse_arguments()
    outpath = pathlib.Path(outputs)
    outpath.mkdir(exist_ok=True, parents=True)
    
    files = pd.Series(pathlib.Path(inputs).rglob("*.csv"))
    files = files.to_frame()
    files.columns = ["filename"]
    files["sex"] = files.filename.astype(str).str.extract("only.(\w+)s")
    files = files.dropna().copy()
    files["GOtype"] = files.filename.astype(str).str.extract("(MF|BP|CC)")[0]
    files = files.dropna().copy()
    files = files.set_index(["sex", "GOtype"]).filename

    data = pd.concat(files.transform(pd.read_csv, usecols=["ID", "p.adjust"]).to_dict())
    data = data.reset_index(level=-1, drop=True)
    tqdm.tqdm.pandas()
    print(
        data.groupby(level=params["level"]).apply(
            run_gofigure,
            commandlineargs=params["commandline"] if params["commandline"] else [],
            outputdir=outpath,
        )
    )


if __name__ == "__main__":
    main()
