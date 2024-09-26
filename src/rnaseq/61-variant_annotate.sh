#!/bin/bash

bcftools annotate -a $1 -c CHROM,FROM,TO,GENE -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') $2 > $3