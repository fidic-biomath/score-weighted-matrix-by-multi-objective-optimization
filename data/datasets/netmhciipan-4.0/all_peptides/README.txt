http://www.cbs.dtu.dk/suppl/immunology/NAR_NetMHCpan_NetMHCIIpan/
consulted: June 30th 2021

NetMHCIIpan-4.0
Here, you will find the data set used for training of NetMHCIIpan-4.0.

NetMHCIIpan_train.tar.gz

Download the file and untar the content using

cat NetMHCIIpan_train.tar.gz | uncompress | tar xvf -
This will creat the directory called NetMHCIIpan_train. In this directory you will find 22 files. 10 files (train_BA?.txt, test_BA?.txt) with partitions with binding affinity data, and 10 files (train_EL?.txt, test_EL?.txt) with eluted ligand data. The format for each file is

AAASVPAADKFKTFE 0.203668 HLA-DPA10103-DPB10201 XXXAAATFEXXX 
APEVKYTVFETALKK 0.838333 HLA-DPA10103-DPB10201 XXXAPELKKXXX 
ATFEAMYLGTCKTLT 0.325328 HLA-DPA10103-DPB10201 XXXATFTLTXXX 
AVWVDGKARTAWVDS  0.14783 HLA-DPA10103-DPB10201 XXXAVWVDSXXX 
ELYYAIYKASPTLAF 0.617078 HLA-DPA10103-DPB10201 XXXELYLAFXXX 
ENVIDVKLVDANGKL 0.173508 HLA-DPA10103-DPB10201 XXXENVGKLXXX 
where the different columns are peptide, target value, MHC_molecule/cell-line, and context. In cases where the 3rd columns is a cell-line ID, the MHC molecules expressed in the cell-line are listed in the allelelist.txt file.
The allelelist.txt file contains the information about alleles expressed in each MA cell-line data set, and pseudosequence.2016.all.X.dat the MHC pseudo sequenes for each MHC molecule.
