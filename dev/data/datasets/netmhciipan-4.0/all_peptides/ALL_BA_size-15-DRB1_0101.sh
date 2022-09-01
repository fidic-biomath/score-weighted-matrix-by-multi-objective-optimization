cat test_BA1_size-15-DRB1_0101.txt test_BA2_size-15-DRB1_0101.txt test_BA3_size-15-DRB1_0101.txt test_BA4_size-15-DRB1_0101.txt test_BA5_size-15-DRB1_0101.txt train_BA1_size-15-DRB1_0101.txt train_BA2_size-15-DRB1_0101.txt train_BA3_size-15-DRB1_0101.txt train_BA4_size-15-DRB1_0101.txt train_BA5_size-15-DRB1_0101.txt > ALL_BA_size-15-DRB1_0101.txt

cat ALL_BA_size-15-DRB1_0101.txt | sort | uniq -c | awk '{print $2}' > peptide_BA_size-15-DRB1_0101.txt
