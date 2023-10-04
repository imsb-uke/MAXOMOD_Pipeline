#!/bin/sh

mkdir -p "database/uniprot"

export ftp_proxy=${http_proxy}

wget \
  "https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz" \
  -O "database/uniprot/uniprot_HUMAN_9606_idmapping_selected.tab.gz"

wget \
  "https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/idmapping/by_organism/MOUSE_10090_idmapping_selected.tab.gz" \
  -O "database/uniprot/uniprot_MOUSE_10090_idmapping_selected.tab.gz"
