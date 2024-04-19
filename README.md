# README: MAXOMOD Pipeline

## Multiomic ALS signatures highlight sex differences and molecular subclusters and identify the MAPK pathway as therapeutic target
This repository contains the analytical pipeline for the MAXOMOD project, which focuses on the multi-omic analysis of axono-synaptic degeneration in motor neuron diseases, particularly ALS. The project explores sex differences and molecular subclusters in ALS and investigates the MAPK pathway as a potential therapeutic target.

For a detailed understanding of the scientific background and the findings, refer to our paper published on bioRxiv [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.08.14.553180v1).

## Table of Contents

- [Getting Started](#getting-started)
  - [Prerequisites](#prerequisites)
- [Data Preparation](#data-preparation)
  - [Organize Data](#organize-data)
- [Reproducing Results](#reproducing-results)
- [Contributing](#contributing)
- [Contact](#contact)


## Getting Started
### Prerequisites
- Git
- DVC

## Data Preparation

Human sequencing data: [EGAS00001007318](https://ega-archive.org/datasets/) [due to patient data, access is restricted]

Mouse sequencing data: [GSE234246](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234246)

Proteomics data: [PXD043300](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD043300)

Phosphoproteomics data: [PXD043297](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD043297)

### Organize Data

All data should be organized in datasets using the following structures

```
datasets/
    consortium/
        <model>/
            01_received_data/
                <omic>/
            cohort/
```

The pipeline expects `fastq.gz` files for the sequencing data und `txt` files for the proteomics data and `csv` files for the phosphoproteomics data.

Please, use [DVC](https://dvc.org/) to see which exact file names are required:

```bash
dvc status srna_organize_samples proteomics_organize_samples phosphoproteomics_organize_samples rnaseq_nextflow
```


## Reproducing Results

To reproduce the analysis results, execute the following command:

```bash
dvc repro
```

This command will run the predefined pipelines to process and analyze the data according to the methodology described in the associated publication.

## Contributing

We welcome contributions to enhance the reproducibility and scope of the analysis. Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct and the process for submitting pull requests.

Clone the git repository:

    ```
    git clone https://github.com/imsb-uke/MAXOMOD_Pipeline.git ./maxomod
    ```

Enter the cloned directory:

   ``` cd maxomod ```

## Contact

For questions or collaboration offers, please contact the project's principal investigators via email provided on the MAXOMOD project page: [MAXOMOD Contact Information](https://www.gesundheitsforschung-bmbf.de/de/maxomod-multi-omische-analyse-axono-synaptischer-degeneration-bei-motoneuronerkrankungen-9409.php).


