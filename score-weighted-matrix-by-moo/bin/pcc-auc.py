import corr1
import corr2
import corr3
#import pdb
import sys
import numpy as np
import pandas as pd
import matplotlib
from sklearn.metrics import roc_curve, accuracy_score, roc_auc_score
import matplotlib
matplotlib.use('Qt5Agg')
import  matplotlib.pyplot as plt

global auc
global pcc
global pred_scores

def main():
    print("*")
    print("*Starting evaluation")

    data_peptides_filename = sys.argv[1]
    score_matrix_filename = sys.argv[2]
    #pocket_weights = sys.argv[3]
    individual=[1,1,1,1,1,1,1,1,1]

    # PEPTIDE FILE SIZE
    dataset_size = len(open(data_peptides_filename).readlines())
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
    #   data_peptides_filename:
    #   score_matrix_filename:
    #
    # outputs:
    #   array_r <class 'numpy.ndarray'>:  score_matrix
    #   label <class 'bytes'>:   columns labels of array_r
    #   index_n <class 'int'>:  number of 9-mer possible peptides
    array_r, label, index_n = corr2.store_nsga2(dataset_size, data_peptides_filename, score_matrix_filename)
    #print( array_r )
    #print( label )
    #print( index_n )

    data_peptide_core_chars = np.array((index_n,9), dtype='U')
    hasha = np.array((index_n), dtype='U')
    exp_scores = np.zeros((dataset_size), dtype='f')

    # CORR3
    # ??????
    # inputs:
    #   dataset_size:
    #   index_n:
    #   data_peptides_filename: 
    #   score_matrix_filename:
    # outputs:
    #   data_peptide_core_chars <class 'numpy.ndarray'>:  9-mer possible peptides by chars
    #   hasha <class 'numpy.ndarray'>:  ????
    #   exp_scores <class 'numpy.ndarray'>:  experimental peptide binding scores 
    data_peptide_core_chars, hasha, exp_scores = corr3.store_nsga2(dataset_size, index_n, data_peptides_filename, score_matrix_filename)
    #np.set_printoptions(threshold=sys.maxsize)
    #print( type(data_peptide_core_chars) )
    #print( data_peptide_core_chars.shape )
    #print( data_peptide_core_chars )
    #print( type(hasha) )
    #print( hasha.shape )
    #print( hasha )
    #print( type(exp_scores) )
    print( exp_scores.shape )
    #print( exp_scores )
    # CORR1
    # 
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
    #
    # Debugging
    #pdb.set_trace()
    pcc, auc, pred_scores = corr1.corr_nsga2(individual, label, data_peptide_core_chars, hasha, array_r, exp_scores, index_n, dataset_size)
    #print( type(pred_scores) )
    print( pred_scores.shape )
    #print( pred_scores )
    print("PCC_value: %f" % pcc)
    print("AUC_value: %f" % auc)
    #score_matrix = pd.read_csv(score_matrix_filename, delim_whitespace=True)
    #print(score_matrix)
    #data_peptides = pd.read_csv(data_peptides_filename, delim_whitespace=True, header=None)
    #print(data_peptides)

    ###################################################
    my_array = np.array([exp_scores,pred_scores])
    #print(my_array.T)
    df = pd.DataFrame(my_array.T, columns=['exp_scores', 'pred_scores'])

    #..."where the first column gives the peptide, the second column the log50k transformed binding affinity (i.e. 1 - log50k( aff nM)), and the last column the class II allele. When classifying the peptides into binders and non-binders for calculation of the AUC values for instance, a threshold of 500 nM is used. This means that peptides with log50k transformed binding affinity values greater than 0.426 are classified as binders." Source: https://services.healthtech.dtu.dk/service.php?NetMHCIIpan-3.2  (march 24th 2022)
    df.loc[df['exp_scores'] > 0.426,'exp_binders' ]=1
    df.loc[df['exp_scores'] <= 0.426,'exp_binders' ]=0

    # apply normalization min max 
    df['pred_scores_min_max'] = (df['pred_scores'] - df['pred_scores'].min()) / (df['pred_scores'].max() - df['pred_scores'].min())
    print(df)

    ################################################################################
    # ROC dataframe
    fpr, tpr, thresholds = roc_curve(df['exp_binders'],df['pred_scores_min_max'])
    roc_df = pd.DataFrame({'recall':tpr,'specificity': 1-fpr})
    #rint(roc_df)

    #print(df.corr().iloc[[0],[1]])
    #pcc1 = df.corr().iloc[[0],[1]].to_numpy()
    pcc1 = df.corr().iloc[[0],[1]].values
    auc1 = np.sum(roc_df.recall[:-1] * np.diff(1 - roc_df.specificity))
    auc2 = roc_auc_score( df['exp_binders'],df['pred_scores_min_max'] )
    print("pandas.PCC_value: %f" % pcc1)
    print("pandas.AUC_value1: %f" % auc1)
    print("pandas.AUC_value2: %f" % auc2)

    ##################################################################################
    # Plot ROC curve
    ax = roc_df.plot(x='specificity', y='recall', figsize=(4,4), legend=False)
    ax.set_ylim(0, 1)
    ax.set_xlim(1, 0)
    ax.plot((1, 0), (0, 1))
    ax.set_xlabel('specificity')
    ax.set_ylabel('recall')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
