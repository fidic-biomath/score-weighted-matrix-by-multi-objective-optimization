#!/bin/bash

./select_by_size.sh train_BA1.txt 15 > train_BA1_size-15.txt
./select_by_size.sh train_BA2.txt 15 > train_BA2_size-15.txt
./select_by_size.sh train_BA3.txt 15 > train_BA3_size-15.txt
./select_by_size.sh train_BA4.txt 15 > train_BA4_size-15.txt
./select_by_size.sh train_BA5.txt 15 > train_BA5_size-15.txt

./select_by_size.sh test_BA1.txt 15 > test_BA1_size-15.txt
./select_by_size.sh test_BA2.txt 15 > test_BA2_size-15.txt
./select_by_size.sh test_BA3.txt 15 > test_BA3_size-15.txt
./select_by_size.sh test_BA4.txt 15 > test_BA4_size-15.txt
./select_by_size.sh test_BA5.txt 15 > test_BA5_size-15.txt
