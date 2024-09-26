#!/bin/bash

zcat $1 | \
    awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$14,$10,$7}}' | \
    tr -d '";' | \
    bgzip > $2