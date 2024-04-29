import pandas as pd
import pathlib
import gzip
import shutil
import tqdm.auto as tqdm
import multiprocessing as mp


def _remove_prefix(string):
    if string.startswith('"'):
        string = string[1:]
    if string.endswith('"'):
        string = string[:-1]
    return string


def parse_series_matrix(file):
    res = {}
    with gzip.open(file, "rt") as gfile:
        for line in gfile:
            line = line[1:].split("\t")
            if line[0].startswith("Sample"):
                res[line[0][7:]] = line[1:]
    res = pd.DataFrame(res).apply(lambda x: x.str.strip())
    return res.applymap(_remove_prefix)


def gzip_copy(files):
    infile, outfile = files
    with open(infile, "rb") as f_in:
        with gzip.open(outfile, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


def main():
    mapping = pd.read_csv("datasets/consortium/00_download/SraRunTable.txt")
    data = pd.concat(
        (
            parse_series_matrix(
                "datasets/consortium/00_download/GSE234246-GPL17021_series_matrix.txt.gz"
            ),
            parse_series_matrix(
                "datasets/consortium/00_download/GSE234246-GPL24247_series_matrix.txt.gz"
            ),
        )
    )
    data = pd.merge(mapping, data, left_on="Library Name", right_on="geo_accession")
    data[["Model", "Condition", "Sex", "SampleID"]] = data.title.str.split(
        ",", expand=True
    )
    data = data[
        ["Assay Type", "Run", "Model", "Condition", "Sex", "SampleID", "geo_accession"]
    ]
    pbar = tqdm.tqdm(total=len(data))
    for model, group in data.groupby("Model"):
        samplesheet = group.drop(columns=["Run", "Model"]).set_index("SampleID")
        samplesheet = samplesheet.groupby(level=0).apply(
            lambda x: x.apply(lambda y: ";".join(y.drop_duplicates().values))
        )
        rawpath = pathlib.Path(
            "datasets/consortium/", f"{model}-mouse", "01_received_data"
        )
        rawpath.mkdir(exist_ok=True, parents=True)
        samplesheetpath = rawpath.joinpath("cohort", "sample_annotations.csv")
        samplesheetpath.parent.mkdir(exist_ok=True, parents=True)
        samplesheet.to_csv(samplesheetpath)
        for assaytype, foldername in {"miRNA-Seq": "srna", "RNA-Seq": "rnaseq"}.items():
            tmp = group[group["Assay Type"] == assaytype].set_index("Run")["SampleID"]
            folderpath = rawpath.joinpath(foldername)
            folderpath.mkdir(exist_ok=True, parents=True)
            args = []
            for runid, sampleid in tmp.items():
                for file in pathlib.Path("datasets/consortium/00_download/fastq").glob(
                    f"{runid}*"
                ):
                    new_name = folderpath.joinpath(
                        file.name.replace(runid, sampleid) + ".gz"
                    )
                    args.append((file, new_name))
            with mp.Pool(processes=8) as pool:
                list(pool.imap_unordered(gzip_copy, args))
            pbar.update(len(tmp))
    pbar.close()


if __name__ == "__main__":
    main()
