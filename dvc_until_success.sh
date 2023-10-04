#!/bin/bash

# Sometimes the pipeline crashes because some R pkg has
# some memory corruption issues.
# Because the problem is random, retrying usually works
# This script runs `dvc repro` non-stop until it finishes
# successfully

echo `date -R`

while true; do
  dvc repro

  if [ "$?" = "0" ]; then
    exit 0
  fi
  echo "Retrying after 5 seconds..."
  sleep 5
  echo `date -R`	

done
