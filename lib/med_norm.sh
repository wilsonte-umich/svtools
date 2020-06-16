#!/bin/bash

BGZ_FILE=$1
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

if [ -f "$BGZ_FILE" ]; then
    zcat $BGZ_FILE |
    Rscript $DIR/med_norm_bed.R |
    bgzip -c > $BGZ_FILE.tmp
    mv $BGZ_FILE.tmp $BGZ_FILE
    tabix -p bed $BGZ_FILE
else
    echo "med_norm.sh: bad file"
    exit 100
fi
