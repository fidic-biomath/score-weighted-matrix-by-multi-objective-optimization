import corr1
import corr2
import corr3
import sys
import numpy as np
import pandas as pd

global pred_scores

def main():
    np.set_printoptions(threshold=sys.maxsize)
    print("*")
    print("*********Setup Parameters**************")

    data_peptides_filename = sys.argv[1]
    score_matrix_filename = sys.argv[2]
    if len(sys.argv) > 3:
        individual_filename = sys.argv[3]
        individual=np.genfromtxt(individual_filename)
        threshold_f1_score_max=float(sys.argv[4])
    else:
        individual=[1,1,1,1,1,1,1,1,1]
        threshold_f1_score_max=float(0.5)

    # PEPTIDE FILE SIZE
    dataset_size = len(open(data_peptides_filename).readlines())

    ###########################################################################
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
    #############################################################   
    np.set_printoptions(threshold=sys.maxsize)
    print("*")
    print("*Score matrix name: ", score_matrix_filename)
    print("*Score matrix:")
    for i in array_r:
        for j in i:
            print(j, end=" ")
        print()
    print("*Pocket Weights: ", individual)
    print("*Data set file: ", data_peptides_filename)
    print("*Number of data peptides: ", dataset_size)
    print("*Number of 9-mer possible peptides: ", index_n)
    print("*")
    print("**********Starting Evaluation*************")
    print("***********Binding Predictor**************")
#    print("*")
#    print("*peptide_number core_number core_sequence calculated_score exp_score")
#    print("*")
    ##########################################################
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
    pcc, auc, pred_scores, pred_core_index = corr1.corr_nsga2(individual, label, data_peptide_core_chars, hasha, array_r, exp_scores, index_n, dataset_size)
    ########################################################################
    # Print of:
    # *peptide_number core_number core_sequence calculated_score exp_score
    k=1
    k_list=[]
    i_list=[]
    core_list=[]
    pred_scores_list=[]
    #exp_scores_list=[]
    for i in pred_core_index:
        core_chars=[x.decode() for x in data_peptide_core_chars[i-1]] # Se usa i-1 por que de Fortran (corr1) viene el conteo de pred_core_index a partir de "1". 
#        print(k, i, ''.join(core_chars), pred_scores[k-1], exp_scores[k-1])
        k_list.append(k)
        i_list.append(i)
        core_list.append(''.join(core_chars))
        pred_scores_list.append(pred_scores[k-1])
    #    exp_scores_list.append(exp_scores[k-1])
        k+=1
    data_df=pd.DataFrame({'peptide_num':k_list,'core_num':i_list,'core':core_list,'pred_scores':pred_scores_list})
    data_df.to_csv("results_score-predict.log", encoding='utf-8', index=False, sep=' ')
    ##########################################################################
#    data_df = pd.DataFrame({'exp_scores':exp_scores,'pred_scores':pred_scores})
    ###########################################################################
    #..."where the first column gives the peptide, the second column the log50k transformed binding affinity (i.e. 1 - log50k( aff nM))    ,
    # and the last column the class II allele. When classifying the peptides into binders and non-binders for calculation of the
    # AUC values for instance, a threshold of 500 nM is used. This means that peptides with log50k transformed binding affinity
    # values greater than 0.426 are classified as binders."
    # Source: https://services.healthtech.dtu.dk/service.php?NetMHCIIpan-3.2  (march 24th 2022)
    #data_df.loc[data_df['exp_scores'] >  0.426,'exp_binders' ]=1
    #data_df.loc[data_df['exp_scores'] <= 0.426,'exp_binders' ]=0
    # apply normalization min max in "pred_scores"
    data_df['pred_scores_max_min'] = (data_df['pred_scores'] - data_df['pred_scores'].min()) / (data_df['pred_scores'].max() - data_df['pred_scores'].min())
    print("**********Normalization*******************")
    print(f"Pred Score: max= {data_df['pred_scores'].max()} min= {data_df['pred_scores'].min()}")
    ##################################################################################
    # Threshold F1-Score-max
    print("******************************************")
    print(f"Threshold: {threshold_f1_score_max}")
    ###############################################################################
    # Filtering resulst dataframe by f1-score-max threshold
    print("******************************************")
    print( data_df[ data_df['pred_scores_max_min'] > threshold_f1_score_max ][['peptide_num','core_num','core','pred_scores','pred_scores_max_min']].sort_values(by='pred_scores_max_min', ascending=False).to_string(index=False) )

if __name__ == "__main__":
    main()