#!/bin/bash

./select_by_allele.sh train_BA1_size-15.txt DRB1_0101 > train_BA1_size-15-DRB1_0101.txt
./select_by_allele.sh train_BA2_size-15.txt DRB1_0101 > train_BA2_size-15-DRB1_0101.txt
./select_by_allele.sh train_BA3_size-15.txt DRB1_0101 > train_BA3_size-15-DRB1_0101.txt
./select_by_allele.sh train_BA4_size-15.txt DRB1_0101 > train_BA4_size-15-DRB1_0101.txt
./select_by_allele.sh train_BA5_size-15.txt DRB1_0101 > train_BA5_size-15-DRB1_0101.txt

./select_by_allele.sh test_BA1_size-15.txt DRB1_0101 > test_BA1_size-15-DRB1_0101.txt
./select_by_allele.sh test_BA2_size-15.txt DRB1_0101 > test_BA2_size-15-DRB1_0101.txt
./select_by_allele.sh test_BA3_size-15.txt DRB1_0101 > test_BA3_size-15-DRB1_0101.txt
./select_by_allele.sh test_BA4_size-15.txt DRB1_0101 > test_BA4_size-15-DRB1_0101.txt
./select_by_allele.sh test_BA5_size-15.txt DRB1_0101 > test_BA5_size-15-DRB1_0101.txt
