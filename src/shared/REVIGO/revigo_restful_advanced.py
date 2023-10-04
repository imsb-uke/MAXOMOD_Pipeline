#!/usr/bin/env python3

# Adapted from http://revigo.irb.hr/CodeExamples/revigo_restful.py.txt (excessed 2022-06-09)

# Python script for programtic access to Revigo. Run it with (last output file name is optional):
# python3 revigo.py example.csv result.csv

import requests
import time
import sys
import pandas
import yaml
import argparse
import pathlib

def parse_arguments():
    """
    Parses command line arguments
    """
    parser = argparse.ArgumentParser(description="Run Tool")
    parser.add_argument(
        "--dataset_name",
        required=True,
        type=str,
        help="Dataset name for individual mouse model or human data.",
    )
    parser.add_argument(
        "--omics_type",
        required=True,
        type=str,
        help="Dataset name for individual mouse model or human data.",
    )
    args = parser.parse_args()
    print(args)
    return args.dataset_name, args.omics_type

def main():
    dataset_name, omics_type = parse_arguments()
    print(dataset_name)

    #getting settings
    with open("params.yaml", "r") as file:
        params = yaml.safe_load(file)[omics_type]["revigo"]["stages"][dataset_name]["settings"]
    
    for analysis_name in params["input"].keys():
        print(analysis_name)
        output = pathlib.Path(params["output_directory"]).joinpath(analysis_name)
        output.mkdir(exist_ok=True, parents=True)
        
        for ontology_name in params["input"][analysis_name].keys():
            print(ontology_name)
            
            # Read enrichments file
            #userData = open(sys.argv[1],'r').read()
            #userData = open('./revigo_example.txt','r').read()
            userData = pandas.read_csv(params["input"][analysis_name][ontology_name])

            if userData.shape[0] == 0:
                continue 

            userData= userData[["ID","p.adjust"]]
            userData=userData.to_csv(header=False,index=False,sep='\t')
            print(userData)

            # Submit job to Revigo
            #The list of currently supported species:
            #NCBI taxon ID	Species name
            #0	Whole UniProt database (default)
            #3702	Arabidopsis thaliana
            #224308	Bacillus subtilis subsp. subtilis str. 168
            #9913	Bos taurus
            #6239	Caenorhabditis elegans
            #3055	Chlamydomonas reinhardtii
            #3827	Cicer arietinum
            #7955	Danio rerio
            #44689	Dictyostelium discoideum
            #7227	Drosophila melanogaster
            #83333	Escherichia coli K-12
            #9031	Gallus gallus
            #9606	Homo sapiens
            #10090	Mus musculus
            #83332	Mycobacterium tuberculosis H37Rv
            #39947	Oryza sativa Japonica Group
            #36329	Plasmodium falciparum 3D7
            #208964	Pseudomonas aeruginosa PAO1
            #10116	Rattus norvegicus
            #559292	Saccharomyces cerevisiae S288C
            #284812	Schizosaccharomyces pombe 972h-
            #4081	Solanum lycopersicum
            #1111708	Synechocystis sp. PCC 6803 substr. Kazusa
            #31033	Takifugu rubripes
            #8364	Xenopus tropicalis
            #4577	Zea mays
        
            if params["organism"]=="hsapiens":
                speciesTaxon='9606'
            elif params["organism"]=="mmusculus":
                speciesTaxon='10090'
            else:
                speciesTaxon='0'

            payload = {'cutoff':'0.7', 'valueType':'pvalue', 'speciesTaxon':speciesTaxon, 'measure':'SIMREL', 'goList':userData}
            r = requests.post("http://revigo.irb.hr/StartJob.aspx", data=payload)

            jobid = r.json()['jobid']

            # Check job status
            running = 1
            while (running!=0):
                r = requests.post("http://revigo.irb.hr/QueryJobStatus.aspx", data={'jobid':jobid})
                running = r.json()['running']
                time.sleep(1)

            # Fetch results
            #Gather the job results at http://revigo.irb.hr/ExportJob.aspx with the following parameters:

            #    "jobid" - Job ID that you collected in the first step.
            #    "namespace" - [1, 2, 3] The namespace for which you are collecting results. 1 represents BIOLOGICAL_PROCESS, 2 represents CELLULAR_COMPONENT, 3 represents MOLECULAR_FUNCTION.
            #    "type" - [csvtable, table, rtable, xgmml, csvtree, rtree, simmat] The type of the output that you require.
            #    - CSVTable (obsolete) gets the resulting table with Scatterplot data in a comma separated values format.
            #    - Table gets the resulting table with scatter plot data (2D and 3D) in a tab separated values format.
            #    - RTable gets the R script for Scatterplot.
            #    - Xgmml gets the xgmml for Cytoscape.
            #    - CSVTree gets the Tree Map data in a comma separated values format.
            #    - RTree gets the R script for Tree Map.
            #    - SimMat gets the Term Similarity Matrix in a tab separated values format.

            #ExportJob method returns text object in the format you specified, or a single line beginning with "Error:" if the error occured, parameters are invalid or no data is available for a given namespace.
            #The server will keep the results for 15 minutes since your last Job command.
            
            if ontology_name=="BP":
                namespace='1'
            elif ontology_name=="CC":
                namespace='2'
            elif ontology_name=="MF":
                namespace='3'
            else:
                raise Exception("Not correct namespace! Must be of ['BP','CC','MF']")
                
            for type_export in ['RTable','Xgmml']:
                r = requests.post("http://revigo.irb.hr/ExportJob.aspx", data={'jobid':jobid, 'namespace': namespace, 'type': type_export})
                if type_export=='RTable':
                    outfile_name = f"{dataset_name}_{analysis_name}_{ontology_name}_revigo_scatterplot.r"
                elif type_export=='Xgmml':
                    outfile_name = f"{dataset_name}_{analysis_name}_{ontology_name}_revigo_cytoscape.xgmml"
                else:
                    raise Exception("type_export not supported! One of [csvtable, table, rtable, xgmml, csvtree, rtree, simmat]")
                # Write results to a file - if file name is not provided the default is result.csv
                with open(output.joinpath(outfile_name) ,'w') as f:
                    f.write(r.text)

if __name__ == "__main__":
    main()
