import pandas as pd
import yaml
import pathlib


def read_files(file):
    return len(pd.read_csv(file, index_col=0))


def main():

    with open("params.yaml", "r") as file:
        params = yaml.safe_load(file)["phosphoproteomics"]["deg_count"]
    files = pd.Series(params)
    counts = files.map(read_files).to_dict()
    outpath = pathlib.Path("results/integration/phosphoproteomics/deg_counts.yaml")
    outpath.parent.mkdir(exist_ok=True, parents=True)
    with outpath.open("w") as file:
        yaml.dump(counts, file)


if __name__ == "__main__":
    main()
