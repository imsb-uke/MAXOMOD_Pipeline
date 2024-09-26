#!/usr/bin/env python3
import argparse
import pathlib
from typing import Optional

import yaml
import subprocess


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
        params["stages"][args.stage].get("outputs", None),
    )


def main():
    params, inputs, outputs = parse_arguments()
    outputs = pathlib.Path(outputs)
    outputs.parent.mkdir(exist_ok=True, parents=True)

    command = [
        "htseq-count",
        "-m",
        params["method"],
        "-f",
        "bam",
        "-t",
        params["type"],
        "-s",
        params["stranded"],
        inputs["bam"],
        inputs["annotation"],
    ]
    with outputs.open("w") as file:
        subprocess.check_call(command, stdout=file, shell=False)


if __name__ == "__main__":
    main()
