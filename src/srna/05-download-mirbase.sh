#!/bin/sh

mkdir -p "database/mirbase/hsapiens"
mkdir -p "database/mirbase/mmusculus"

export ftp_proxy=${http_proxy}

wget "ftp://mirbase.org/pub/mirbase/22.1/genomes/mmu.gff3" -O "database/mirbase/mmusculus/mirbase_mmu.gff3"

wget "ftp://mirbase.org/pub/mirbase/22.1/genomes/hsa.gff3" -O "database/mirbase/hsapiens/mirbase_hsa.gff3"
