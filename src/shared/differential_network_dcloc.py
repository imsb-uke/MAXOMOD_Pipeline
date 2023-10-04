#!/usr/bin/env python3

import pandas as pd
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import rpy2.robjects as ro
import argparse
import yaml
import pathlib


DCloc = """
#########################################
##########  Helper functions  ###########
#########################################

### Significance test for Fisher-transformed correlation coefficients
# n: number of samples
# x: first correlation to compare
# y: second correlation
#
c.test <- function(n, x, y) {
    stan.dev <- sqrt(2/(n - 3))
    dist <- x-y
    z <- abs(dist/stan.dev)
    p <- 2*(1 - pnorm(z))
    return(c(p,z))
}

### Fisher transformation
fisherr2z <- function(r) {
  z <- log((1 + r)/(1 - r)) / 2
  return(z)
}

### Inverse Fisher transformation
fisherz2r <- function(z) {
  r <- (exp(2 * z) - 1)/(1 + exp(2 * z))
  return(r)
}

#########################################
##########  DCloc()  ####################
#########################################
# DESCRIPTION
# This functions calculates the differential correlation between disease
# condition A and B by comparison of the local topology of correlation
# networks.
#
# ARGUMENTS
# Mat.A: gene expression matrix for disease condition A.
# Mat.B: gene expression matrix for disease condition B.
# Mat.A and Mat.B should have the same dimension.
# Rows = genes, columns = samples.
# n.supp: number of supporting points between r.min and r.max.
# r.min: minimal correlation to be considered.
# r.max: maximal correlation to be considered.
# min.neigh: minimum number of neigbors of a gene to be considered as
# differentially correlated.
#
# VALUE
# The function returns a matrix with the differential correlation of each gene.
# The colums of the matrix contain the absolute topological dissimilarity,
# the topological dissimilarity and the mean number of network neighbors for each
# of the disease conditions.
#
DCloc  <- function(cor.A, cor.B, nsupp=100, minneigh=3) {
    n.supp=nsupp
    min.neigh = minneigh
    r.min=0
    r.max=fisherz2r(2.5)

    ngenes <- nrow(cor.A)
    for(i in 1:ngenes){
        cor.A[i,i]<-0
        cor.B[i,i]<-0
    }
    cor.A <- fisherr2z(cor.A)
    cor.B <- fisherr2z(cor.B)
    z.min <- fisherr2z(r.min)
    z.max <- fisherr2z(r.max)
    z <- seq(z.min, z.max, length.out=n.supp)

    #Computing the topological dissimilarity
    result <- list()
    for (x in c("nA", "nB", "one", "both", "d")) {
        result[[x]] <- matrix(nrow=ngenes, ncol=n.supp, 0)
        rownames(result[[x]]) <- colnames(cor.A)
        colnames(result[[x]]) <- z
    }

    for (j in 1:n.supp) {
        z.thres <- z[j]
        network.A <- (cor.A >= z.thres)
        network.B <- (cor.B >= z.thres)
        for (i in 1:ngenes) {
            index.A <- which(network.A[i, ]==1)
            index.B <- which(network.B[i, ]==1)
            result$nA[i, j] <- length(index.A)
            result$nB[i, j] <- length(index.B)
            result$both[i, j] <- length(intersect(index.A, index.B))
            result$one[i, j] <- length(union(index.A, index.B))

            if( result$one[i, j] >= min.neigh ) {
                 result$d[i, j] <-(( 1 - result$both[i, j]/result$one[i, j] )* sign(result$nA[i, j] - result$nB[i, j]))
             }

        }
    }

    ret <- matrix(ncol=4,nrow=ngenes,0)
    colnames(ret) <- c("abs.top.dissim", "top.dissim", "mean.neigh.A", "mean.neigh.B")
    rownames(ret) <- colnames(cor.A)
    ret[ ,1] <- apply(result$d, 1, function(x){mean(abs(x))})
    ret[ ,2] <- apply(result$d, 1, mean)
    ret[ ,3] <- apply(result$nA, 1, mean)
    ret[ ,4] <- apply(result$nB, 1, mean)
    return(ret)
}
"""

DCloc = ro.r(DCloc)


def parse_arguments():
    """
    Parses command line arguments
    """
    parser = argparse.ArgumentParser(description="Run Cosifer")
    # add expected arguments
    parser.add_argument(
        "--dataset-name",
        help="dataset key from params.yaml to use",
    )
    parser.add_argument(
        "--omic",
        help="omics from params.yaml to use",
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()
    params = yaml.safe_load(pathlib.Path("params.yaml").read_text())["cosifernetwork"][
        args.omic
    ]["dcloc"]

    file1 = pathlib.Path(params["stages"][args.dataset_name]["infile_control"])
    file2 = pathlib.Path(params["stages"][args.dataset_name]["infile_treated"])

    dataA = pd.read_csv(file2)

    dataB = pd.read_csv(file1)

    dataA = dataA.pivot_table(
        index="e1", columns="e2", values="intensity", fill_value=0.0
    )

    dataB = dataB.pivot_table(
        index="e1", columns="e2", values="intensity", fill_value=0.0
    )

    idx = dataA.index.union(dataA.columns).union(dataB.index).union(dataB.columns)

    dataA = dataA.reindex(idx, axis=0).reindex(idx, axis=1).fillna(0)

    dataB = dataB.reindex(idx, axis=0).reindex(idx, axis=1).fillna(0)

    dataA = dataA + dataA.T

    dataB = dataB + dataB.T

    with localconverter(ro.default_converter + pandas2ri.converter):
        matA = ro.conversion.py2rpy(dataA)
    with localconverter(ro.default_converter + pandas2ri.converter):
        matB = ro.conversion.py2rpy(dataB)

    res = DCloc(matA, matB, nsupp=params["nsupp"], minneigh=params["minneigh"])

    with localconverter(ro.default_converter + pandas2ri.converter):
        pyres = ro.conversion.rpy2py(res)

    results = pd.DataFrame(
        pyres,
        index=dataA.index,
        columns=("abs.top.dissim", "top.dissim", "mean.neigh.A", "mean.neigh.B"),
    )

    output = pathlib.Path(params["stages"][args.dataset_name]["output_file"])
    output.parent.mkdir(exist_ok=True, parents=True)
    results.to_csv(output)


if __name__ == "__main__":
    main()
