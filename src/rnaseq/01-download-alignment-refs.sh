#!/bin/sh

# https://github.com/mikelove/tximeta/blob/master/inst/extdata/hashtable.csv

mkdir -p "database/align-references/hsapiens"
mkdir -p "database/align-references/mmusculus"

export ftp_proxy=${http_proxy}

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.p13.genome.fa.gz \
 -O database/align-references/hsapiens/GRCh38.p13.genome.fa.gz

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.primary_assembly.genome.fa.gz \
 -O database/align-references/hsapiens/GRCh38.primary_assembly.genome.fa.gz

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.transcripts.fa.gz -O database/align-references/hsapiens/gencode.v37.transcripts.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz -O database/align-references/hsapiens/gencode.v37.annotation.gtf.gz

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M26/GRCm39.genome.fa.gz \
 -O database/align-references/mmusculus/GRCm39.genome.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M26/GRCm39.primary_assembly.genome.fa.gz \
 -O database/align-references/mmusculus/GRCm39.primary_assembly.genome.fa.gz

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M26/gencode.vM26.transcripts.fa.gz -O database/align-references/mmusculus/gencode.vM26.transcripts.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M26/gencode.vM26.annotation.gtf.gz -O database/align-references/mmusculus/gencode.vM26.annotation.gtf.gz

