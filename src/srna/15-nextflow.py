#!/usr/bin/env python3
import argparse
import os
import re
import shutil
import subprocess
import yaml
import collections.abc



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
    parser = argparse.ArgumentParser(description="Run srna nextflow stage")
    # add expected arguments
    parser.add_argument(
        "--dataset-name",
        dest="dataset_name",
        default=None,
        help="Under params.yaml, one of the keys under srna.nextflow.stages (e.g. sod1)",
    )
    # parse args
    args = parser.parse_args()
    # Uncomment to debug the output of argparse:
    # raise ValueError(args)
    return args


def get_settings(omic="srna", stage="nextflow"):
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


def write_nf_config(output_directory, cpus=20, memory="600.GB"):
    nf_conf = "\n".join([
        "process.cpus = " + str(cpus),
        "process.memory = " + str(memory),
        "plugins {\n  id 'nf-amazon'\n}\n",
        "docker.runOptions = " + "'" + " ".join([
            "-u " + str(os.getuid()) + ":" + str(os.getgid()),
        ]) + "'",
    ]) + "\n"
    with open(os.path.join(output_directory, "nextflow.config"), "w") as fh:
        fh.write(nf_conf)


def main(settings):
    basedir = os.getcwd()
    os.makedirs(settings["output_directory"], exist_ok=True)
    write_nf_config(
        settings["output_directory"], cpus=settings["cpus"], memory=settings["memory"]
    )
    os.chdir(settings["output_directory"])
    absolute_input_glob = os.path.join(
        basedir, settings["input_directory"], settings["input_glob"]
    )
    command = [
        settings["nextflow"],
        "run",
        settings["nf_pipeline"],
        "-r",
        settings["nf_pipeline_version"],
        "-profile",
        settings["nf_profile"],
    ]
    command.extend(["--max_memory", settings["memory"]])
    command.extend(["--max_cpus", str(settings["cpus"])])
    command.extend(["--genome", settings["genome"]])
    if settings["mirna_gtf"] is not None:
        mirna_gtf_path = os.path.join(basedir, settings["mirna_gtf"])
        command.extend(["--mirna_gtf", mirna_gtf_path])
    command.extend(["--reads", absolute_input_glob])
    if settings.get("resume"):
        command.append("-resume")
    print("Running command:...")
    print(command)
    # workaround https://github.com/nextflow-io/nextflow/issues/1962
    envir = os.environ
    envir["LC_ALL"] = "C"
    # end of workaround
    subprocess.run(command, check=True, env=envir)


if __name__ == "__main__":
    settings = get_settings()
    main(settings)
