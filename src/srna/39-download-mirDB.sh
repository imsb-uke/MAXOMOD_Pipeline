#!/bin/sh

mkdir -p "database/miRDB"

export ftp_proxy=${http_proxy}
wget \
  "http://mirdb.org/download/miRDB_v6.0_prediction_result.txt.gz" \
  -O "database/miRDB/miRDB_v6.0_prediction_result.txt.gz"
