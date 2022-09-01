#!/bin/bash

FILE=$1
LENGTH=$2
 
awk -v LENGTH=${LENGTH} 'length($1)==LENGTH {print $0}' ${FILE}
