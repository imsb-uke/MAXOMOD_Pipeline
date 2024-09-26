import decoupler as dc
import pathlib
import sys



def main():
    species = sys.argv[1].strip()
    net = dc.get_dorothea(organism=species, levels=["A", "B", "C", "D"])
    net.source = net.source.astype("category")
    net.target = net.target.astype("category")
    net.confidence = net.confidence.astype("category")
    outpath = pathlib.Path(f"database/dorothea/{species}.feather")
    outpath.parent.mkdir(exist_ok=True)
    net.reset_index(drop=True).to_feather(outpath)

if __name__ == "__main__":
    main()
