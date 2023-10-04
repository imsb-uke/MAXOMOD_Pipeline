#!/bin/sh

mkdir -p "database/stringdb"

echo "downloading uniprot to string id mappings:"
export ftp_proxy=${http_proxy}
wget \
  "https://string-db.org/mapping_files/uniprot/mouse.uniprot_2_string.2018.tsv.gz" \
  -O "database/stringdb/mouse.uniprot_2_string.2018.tsv.gz"

wget \
  "https://string-db.org/mapping_files/uniprot/human.uniprot_2_string.2018.tsv.gz" \
  -O "database/stringdb/human.uniprot_2_string.2018.tsv.gz"

wget \
  "https://stringdb-static.org/download/protein.links.v11.0/10090.protein.links.v11.0.txt.gz" \
  -O "database/stringdb/10090.protein.links.txt.gz"

wget \
  "https://stringdb-static.org/download/protein.links.v11.0/9606.protein.links.v11.0.txt.gz" \
  -O "database/stringdb/9606.protein.links.txt.gz"
