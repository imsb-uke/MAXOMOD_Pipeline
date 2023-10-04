#!/usr/bin/env python3
import argparse
import collections.abc
import os
import re
import shutil

import yaml


def update_dict_recursively(d, u):
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            d[k] = update_dict_recursively(d.get(k, {}), v)
        else:
            d[k] = v
    return d

def parse_arguments():
    """
    Parses command line arguments
    """
    parser = argparse.ArgumentParser(description="Run srna organize_samples stage")
    # add expected arguments
    parser.add_argument(
        "--dataset-name",
        dest="dataset_name",
        default=None,
        help="Under params.yaml, one of the keys under srna.organize_samples.stages (e.g. sod1)",
    )
    # parse args
    args = parser.parse_args()
    # Uncomment to debug the output of argparse:
    # raise ValueError(args)
    return args


def get_settings(omic="srna", stage="organize_samples"):
    args = parse_arguments()
    if args.dataset_name is None:
        raise ValueError("Please set --dataset-name")
    with open("params.yaml") as fh:
        params = yaml.safe_load(fh)
    our_settings = {}

    def_settings = params[omic][stage]["default_settings"]
    our_settings = update_dict_recursively(our_settings, def_settings)

    stage_settings = params[omic][stage]["stages"][args.dataset_name].get("settings", {})
    our_settings = update_dict_recursively(our_settings, stage_settings)
    
    def_syssettings = params[omic][stage]["default_system_settings"]
    our_settings = update_dict_recursively(our_settings, def_syssettings)

    syssettings = params[omic][stage]["stages"][args.dataset_name].get("system_settings", {})
    our_settings = update_dict_recursively(our_settings, syssettings)

    return our_settings

def main(settings):
    input_files = os.listdir(settings["input_directory"])
    input_name_pattern = re.compile(settings["input_name_pattern"])
    os.makedirs(settings["output_directory"], exist_ok=True)
    for file in input_files:
        match = input_name_pattern.match(file)
        if match is None:
            continue
        new_filename = re.sub(
            input_name_pattern, settings["output_name_replacement"], file
        )
        source_file = os.path.join(settings["input_directory"], file)
        dest_file = os.path.join(settings["output_directory"], new_filename)
        shutil.copy(source_file, dest_file)


if __name__ == "__main__":
    settings = get_settings()
    main(settings)
