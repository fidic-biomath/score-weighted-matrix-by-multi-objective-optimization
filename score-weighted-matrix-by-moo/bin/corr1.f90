!   INPUT:
!       indv
!       LABEL: first row of score matrix (aminoacid one-letter)
!       TRAIN_CHAR:  peptide 9-mer core sequence by chars
!       HASHA:  ????
!       ARRAY_RR: Score matrix
!       EXP_SCORES:  Experimental data
!       index_n:  number of 9-mer core sequences
!       count_size: dataset size
!   OUPUTS:
!       corr: Pearson Correlation Coeficient
!       auc:  Area Under Curve of Receiver Operating Characteristic (ROC) analysis.
!       PRED_SCORES: Calculated data
 
subroutine corr_nsga2(indv, LABEL, TRAIN_CHAR, HASHA, ARRAY_RR, EXP_SCORES, index_n, count_size, corr, auc, PRED_SCORES)
    !>
    !Initialization of matrices that depends on the sliding_window function
    !
    implicit none

    real*4  indv(9)
    intent(in) indv
    character*20, intent(in) :: LABEL
    character, dimension(index_n, 9), intent(in) :: TRAIN_CHAR
    character*50, dimension(index_n), intent(in) :: HASHA
    real, dimension(9,20), intent(in) :: ARRAY_RR
    real, dimension(count_size), intent(in) :: EXP_SCORES
    integer, intent(in) :: index_n
    integer, intent(in) :: count_size
    real, intent(out) :: corr
    real, intent(out) :: auc
    real, intent(out), dimension(count_size) :: PRED_SCORES

    integer i, aux_n, resol
    integer j, k, tmp_index
    real tmp, mean_x, mean_y
    real cov, var_x, var_y
    real min, max, diff_max, diff_min
    real tp, fn, tn, fp, diff1, diff2, a, b, h
    real threshold_exp, threshold_pred
    real, dimension(count_size) :: TPR, FPR
    real, dimension(index_n) :: PR_SC, PR_SC_N
    real, dimension(index_n,9) :: NUM_SCORES

    !!********* assigning scores*********
    i = 0
    j = 0
    k = 0
  
    DO i = 1, index_n
      DO j = 1, 9
        DO k = 1, 20
          IF (TRAIN_CHAR(i,j) == LABEL(k:k)) THEN
            NUM_SCORES(i,j) = ARRAY_RR(j,k)
          ENDIF
        END DO
      END DO
    END DO

    !*********sum of weights*********
    i = 0
    DO i = 1, index_n
        PR_SC(i) = SUM(NUM_SCORES(i,:)*indv(:))
    END DO

    ! *********stacking to original size*********
    PR_SC_N = PR_SC
    
    i = 1
    j = 1
    k = 0
    aux_n = 1
    tmp_index = 1
    
    DO WHILE (i < index_n)
        DO WHILE (HASHA(i) == HASHA(i+1))
            IF (aux_n  == 1) tmp = PR_SC(i)
            IF (tmp < PR_SC(i+1)) tmp = PR_SC(i+1)
            aux_n = aux_n + 1
            i = i + 1
        END DO
        IF (aux_n > 1) THEN
            PRED_SCORES(j) = tmp
            aux_n = 1
            j = j + 1
            IF (i == index_n) GO TO 100
            tmp = PR_SC(i+1)
            i = i + 1
        ENDIF
        DO WHILE (HASHA(i) /= HASHA(i+1))
            PRED_SCORES(j) = PR_SC(i)
            j = j + 1
            i = i + 1
        END DO
        IF (i == index_n) PRED_SCORES(j) = PR_SC(i)
    END DO
    
    100 corr = 0.0  ! ??????? WTF
    
    ! *********correlating*********
    corr = 0.0
    mean_x = SUM(EXP_SCORES) / count_size
    mean_y = SUM(PRED_SCORES) / count_size
    cov = SUM((EXP_SCORES(1:count_size) - mean_x) * (PRED_SCORES(1:count_size) - mean_y))
    var_x = SUM((EXP_SCORES(1:count_size) - mean_x) * (EXP_SCORES(1:count_size) - mean_x))
    var_y = SUM((PRED_SCORES(1:count_size) - mean_y) * (PRED_SCORES(1:count_size) - mean_y))
    corr = (cov / SQRT(var_x)) / SQRT(var_y)
    
    !*******************ROC******************
    max = 0
    DO i = 1, count_size
        IF (max<PRED_SCORES(i)) max = PRED_SCORES(i)
    END DO
    
    min = max
    DO i = 1, count_size
        IF (min>PRED_SCORES(i)) min = PRED_SCORES(i)
    END DO
      
    threshold_exp = 0.426
    resol = 40
    
    diff_max = max - min
    diff_min = diff_max / resol
    threshold_pred = min 
    
    DO i = 1, resol
        threshold_pred = min + (diff_min * i)
        tp = 0 !true positives
        fn = 0 !false negatives
        tn = 0 !true negatives
        fp = 0 !false positives
        DO j = 1, count_size
            IF (EXP_SCORES(j)>threshold_exp) THEN !below of the experimental threshold is negative
                IF (PRED_SCORES(j)>threshold_pred) THEN ! negative for predicted
                    tn = tn + 1
                ELSE ! positive for predicted
                    fp = fp + 1
                ENDIF
            ELSE !above or equal of the experimental threshold is positive
                IF (PRED_SCORES(j)>threshold_pred) THEN ! negative for predicted
                    fn = fn + 1
                ELSE ! positive for predicted
                    tp = tp + 1
                ENDIF
            ENDIF
        END DO
        TPR(i) = tp / (tp + fn)
        FPR(i) = fp / (tn + fp)
    END DO
    
    auc = 0.0
    diff1 = 0.0
    diff2 = 0.0
    DO i=1, resol
        h = FPR(i)-diff1
        a = TPR(i)
        b = diff2
        auc=(( (a+b)/2 )*h)+auc
        diff1 = FPR(i)
        diff2 = TPR(i)
    END DO
end subroutine corr_nsga2
