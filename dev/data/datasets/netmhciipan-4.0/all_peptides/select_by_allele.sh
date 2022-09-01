#!/bin/bash

FILE=$1
ALLELE=$2
 
awk -v ALLELE=${ALLELE} '$3==ALLELE {print $0}' ${FILE}
