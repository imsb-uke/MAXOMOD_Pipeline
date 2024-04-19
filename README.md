# Multiomic ALS signatures highlight sex differences and molecular subclusters and identify the MAPK pathway as therapeutic target

Reproducible analysis pipeline for the [MAXOMOD](https://www.gesundheitsforschung-bmbf.de/de/maxomod-multi-omische-analyse-axono-synaptischer-degeneration-bei-motoneuronerkrankungen-9409.php) project.

Please check the paper on [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.08.14.553180v1).

## Starting to contribute/work with the pipeline

### Clone the git repository:

    git clone https://github.com/imsb-uke/MAXOMOD_Pipeline.git ./maxomod

Enter the cloned directory:

    cd maxomod

### Download the expected data from public repositories:

Human sequencing data: [EGAS00001007318](https://ega-archive.org/datasets/) [due to patient data, access is restricted]

Mouse sequencing data: [GSE234246](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234246)

Proteomics data: [PXD043300](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD043300)

Phosphoproteomics data: [PXD043297](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD043297)

### Organize the data:

All data should be organized in a datasets using the following structure
`datasets/consortium/<model>/01_received_data/<omic>`
and a metadata directory
`datasets/consortium/<model>/01_received_data/cohort`.

The pipeline expects `fastq.gz` files for the sequencing data und `txt` files for the proteomics data and `csv` files for the phosphoproteomics data.

Please, use [DVC](https://dvc.org/) to see which exact files names are required:

    dvc status srna_organize_samples proteomics_organize_samples phosphoproteomics_organize_samples rnaseq_nextflow

### Reproduce the results

Run [DVC](https://dvc.org/) repro to reproduce the results:

    dvc repro

