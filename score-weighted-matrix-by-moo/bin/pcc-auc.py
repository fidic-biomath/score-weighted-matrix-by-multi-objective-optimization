import array
import random
import json
import numpy as np
import pandas as pd
#from math import sqrt

import sys

import corr1 # 
import corr2 #
import corr3 #

global array_r    #
global data_peptide_core_chars #
global index_n    #
global label      #
global hasha      #
global exp_scores #
global dataset_size #  Number of peptides for training  PCC and AUC

# Funcion de fitness: multiobjetivo 1) pcc 2) auc
#def performance(individual):
    # CORR1
    # ????????
    # inputs:
    #   np.array(individual):
    #   label:
    #   data_peptide_core_chars:
    #   hasha:
    #   array_r:
    #   exp_scores:
    #   index_n:
    #   dataset_size:
    #
    # outputs:
    #   pcc:
    #   auc:
#    pcc, auc = corr1.corr_nsga2(np.array(individual), label, data_peptide_core_chars, hasha, array_r, exp_scores, index_n, dataset_size)
    #pcc = np.person
    #auc = np.roc.auc
#    return pcc, auc

def main():
    print("*")
    print("*Starting evaluation")
    print("")

    data_peptides = sys.argv[1]
    score_matrix = sys.argv[2]
    #pocket_weights = sys.argv[3]
    individual=[1,1,1,1,1,1,1,1,1]

    # PEPTIDE FILE SIZE
    dataset_size = len(open(data_peptides).readlines())
    print("*")
    print("*Number of data peptides: ", dataset_size)
    print("*")
    print("***********Binding Predictor**************")
    print("*Performance indicators: PCC and ROC(AUC)*")
    print("*")

###########################################################################
    array_r = np.zeros((9,20), dtype='f')

    # CORR2 
    # ??????
    # inputs:
    #   dataset_size:
    #   data_peptides:
    #   score_matrix:
    #
    # outputs:
    #  array_r <class 'numpy.ndarray'>:  score_matrix
    #   label <class 'bytes'>:   columns labels of array_r
    #   index_n <class 'int'>:  number of 9-mer possible peptides
    array_r, label, index_n = corr2.store_nsga2(dataset_size, data_peptides, score_matrix)
    print( array_r )
    print( label )
    print( index_n )

    #myFile = pd.read_csv(score_matrix, header=None, delim_whitespace=True)
    #myFile = pd.read_csv(score_matrix, delim_whitespace=True)
    #print(myFile)


    data_peptide_core_chars = np.array((index_n,9), dtype='U')
    hasha = np.array((index_n), dtype='U')
    exp_scores = np.zeros((dataset_size), dtype='f')

    # CORR3
    # ??????
    # inputs:
    #   dataset_size:
    #   index_n:
    #   data_peptides: 
    #   score_matrix:
    # outputs:
    #   data_peptide_core_chars <class 'numpy.ndarray'>:  9-mer possible peptides by chars
    #   hasha <class 'numpy.ndarray'>:  ????
    #   exp_scores <class 'numpy.ndarray'>:  experimental peptide binding scores 
    data_peptide_core_chars, hasha, exp_scores = corr3.store_nsga2(dataset_size, index_n, data_peptides, score_matrix)
    np.set_printoptions(threshold=sys.maxsize)
    print( type(data_peptide_core_chars) )
    print( data_peptide_core_chars.shape )
    print( data_peptide_core_chars )
    print( type(hasha) )
    print( hasha.shape )
    print( hasha )
    print( type(exp_scores) )
    print( exp_scores.shape )
    print( exp_scores )


# CORR1                                                                                                                             
# ????????                                                                                                                          
# inputs:                                                                                                                           
#   np.array(individual):                                                                                                           
#   label:                                                                                                                          
#   data_peptide_core_chars:                                                                                                        
#   hasha:                                                                                                                          
#   array_r:                                                                                                                        
#   exp_scores:                                                                                                                     
#   index_n:                                                                                                                        
#   dataset_size:                                                                                                                   
#                                                                                                                                   
# outputs:                                                                                                                          
#   pcc:                                                                                                                            
#   auc: 

    pcc, auc = corr1.corr_nsga2(np.array(individual), label, data_peptide_core_chars, hasha, array_r, exp_scores, index_n, dataset_size)
    return pcc, auc

if __name__ == "__main__":
    #individual=[1,1,1,1,1,1,1,1,1]
    pcc, auc = main()
    print("PCC value: %d",pcc)
    print("AUC value: %d",auc)
