#!/usr/bin/env python3
import argparse
import copy
import os
import re
import subprocess
import pathlib

import pandas as pd
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
    parser = argparse.ArgumentParser(description="Run nextflow rnaseq pipeline")
    # add expected arguments
    parser.add_argument(
        "--dataset-name",
        dest="dataset_name",
        default=None,
        help="Under params.yaml, one of the keys under rnaseq.nextflow.stages (e.g. sod1)",
    )
    # parse args
    args = parser.parse_args()
    # Uncomment to debug the output of argparse:
    # raise ValueError(args)
    return args


def get_settings():
    args = parse_arguments()
    if args.dataset_name is None:
        raise ValueError("Please set --dataset-name")
    with open("params.yaml") as fh:
        params = yaml.safe_load(fh)
    our_settings = {}

    def_settings = params["rnaseq"]["nextflow"]["default_settings"]
    our_settings = update_dict_recursively(our_settings, def_settings)

    stage_settings = params["rnaseq"]["nextflow"]["stages"][args.dataset_name].get("settings", {})
    our_settings = update_dict_recursively(our_settings, stage_settings)

    def_syssettings = params["rnaseq"]["nextflow"]["default_system_settings"]
    our_settings = update_dict_recursively(our_settings, def_syssettings)

    syssettings = params["rnaseq"]["nextflow"]["stages"][args.dataset_name].get("system_settings", {})
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


def prepare_samplesheet(settings, basedir):
    sample_annot = pd.read_csv(settings["cohort"])
    sample_filename_base = [
        re.sub(
            pattern=settings["sample_to_filename_pattern"],
            repl=settings["sample_to_filename_replacement"],
            string=x,
        )
        for x in sample_annot["SampleID"]
    ]
    fastq_1 = [
        os.path.join(
            basedir,
            settings["input_directory"],
            x + settings["sample_to_filename_suffix_R1"],
        )
        for x in sample_filename_base
    ]
    if settings["sample_to_filename_suffix_R2"] is None:
        fastq_2 = pd.Series(len(sample_annot) * [""])
    else:
        fastq_2 = [
            os.path.join(
                basedir,
                settings["input_directory"],
                x + settings["sample_to_filename_suffix_R2"],
            )
            if pathlib.Path(basedir,
                settings["input_directory"],
                x + settings["sample_to_filename_suffix_R2"]).exists() else ""
            for x in sample_filename_base
        ]
    samplesheet = pd.DataFrame(
        {
            "group": sample_annot["Condition"],
            "replicate": sample_annot.groupby("Condition")["Condition"].transform(
                lambda x: list(range(1, len(x) + 1))
            ),
            "fastq_1": fastq_1,
            "fastq_2": fastq_2,
            "strandedness": pd.Series(len(sample_annot) * [settings["strandedness"]]),
        }
    )
    out_fn = os.path.join(basedir, settings["output_directory"], "samplesheet.csv")
    samplesheet.to_csv(out_fn, index=False)
    out2_fn = os.path.join(
        basedir, settings["output_directory"], "nextflow_ids_to_samples.csv"
    )
    nextflow_ids_to_samples = pd.DataFrame(
        {
            "SampleID": sample_annot["SampleID"],
            "nf_id": [
                str(x) + "_R" + str(y)
                for (x, y) in list(zip(samplesheet.group, samplesheet.replicate))
            ],
        }
    )
    nextflow_ids_to_samples.to_csv(out2_fn, index=False)
    return out_fn


def main(settings):
    basedir = os.getcwd()
    os.makedirs(settings["output_directory"], exist_ok=True)
    write_nf_config(
        settings["output_directory"], cpus=settings["cpus"], memory=settings["memory"]
    )
    command = [settings["nextflow"]] + settings["nextflow_extra_args"]
    command.extend([
        "run",
        settings["nf_pipeline"],
        "-r",
        settings["nf_pipeline_version"],
        "-profile",
        settings["nf_profile"],
    ])
    samplesheet_fn = prepare_samplesheet(settings, basedir)
    command.extend(["--input", os.path.join(basedir, samplesheet_fn)])
    if settings.get("genome") is not None:
        command.extend(["--genome", settings["genome"]])
    else:
        command.extend(["--igenomes_ignore", "true"])
    for setting in ("fasta", "gtf", "star_index", "bed12", "transcript_fasta"):
        if settings.get(setting) is not None:
            command.extend(["--" + setting, os.path.join(basedir, settings[setting])])
    if settings["salmon_transcripts"]:
        command.extend(["--pseudo_aligner", "salmon"])
    if settings.get("save_reference"):
        command.append("--save_reference")
    if settings.get("gencode"):
        command.append("--gencode")
    if settings.get("skip_alignment"):
        command.append("--skip_alignment")
    if settings.get("deseq2_vst"):
        command.append("--deseq2_vst")
    if settings.get("resume"):
        command.append("-resume")
    print("Running command:...")
    print(" ".join(command))
    os.chdir(settings["output_directory"])
    # workaround https://github.com/nextflow-io/nextflow/issues/1962
    envir = os.environ
    envir["LC_ALL"] = "C"
    # end of workaround
    subprocess.run(command, check=True, env=envir)


if __name__ == "__main__":
    settings = get_settings()
    main(settings)
