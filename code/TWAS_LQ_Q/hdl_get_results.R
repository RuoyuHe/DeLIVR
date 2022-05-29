for (i in 1:100){
  setwd("~/deepRIV/UKB/results/TWAS_LQ_Q/hdl")
  library(data.table)
  
  lq_names = c('gene_name','chr',
               'rsq_stage1_X1',
               'adjrsq_stage1_X1',
               'rsq_stage1_X2',
               'adjrsq_stage1_X2',
               'lq_effect_X1',
               'lq_effect_X2',
               'lq_intercept',
               'lq_global_pval',
               'lq_nonlinear_pval',
               'stage1_X1_fstat',
               'stage1_X1_pval',
               'stage1_X2_fstat',
               'stage1_X2_pval',
               'num_snps_X1',
               'num_snps_X2')
  
  q_names = c('gene_name','chr',
              'rsq_stage1_X2',
              'adjrsq_stage1_X2',
              'q_effect',
              'q_intercept',
              'q_pval',
              'stage1_X2_fstat',
              'stage1_X2_pval',
              'num_snps_X2')
  
  poly_names = c(
    'gene_name','chr',
    'rsq_stage1_X1',
    'adjrsq_stage1_X1',
    'poly_lq_effect_X1',
    'poly_lq_effect_X2',
    'poly_lq_intercept',
    'poly_lq_global_pval',
    'poly_lq_nonlinear_pval',
    'poly_q_effect',
    'poly_q_intercept',
    'poly_q_pval',
    'stage1_X1_fstat',
    'stage1_X1_pval',
    'num_snps_X1')
  
  X12_lq_names = c(
    'gene_name','chr',
    # 'rsq_stage1_X12',
    # 'adjrsq_stage1_X12',
    'X12_lq_effect_X1',
    'X12_lq_effect_X2',
    'X12_lq_intercept',
    'X12_lq_global_pval',
    'X12_lq_nonlinear_pval',
    'stage1_X12_fstat',
    'stage1_X12_pval',
    'num_snps_X12')
  
  X12_q_names = c(
    'gene_name','chr',
    # 'rsq_stage1_X12',
    # 'adjrsq_stage1_X12',
    'X12_q_effect',
    'X12_q_intercept',
    'X12_q_pval',
    'stage1_X12_fstat',
    'stage1_X12_pval',
    'num_snps_X12'
  )
  
  file_name = paste0('hdl_',i,'.RData')
  load(file_name)
  n = length(results)
  lq_results = NULL
  q_results = NULL
  poly_results = NULL
  X12_lq_results = NULL
  X12_q_results = NULL
  for (j in 1:n){
    print(paste('i, j, n:',i,j,n))
    if (is.null(results[[j]])) next()
    if ((results[[j]]$num_snps_X2 >= 2)&&(results[[j]]$num_snps_X1 >= 2)&&
        !(is.null(results[[j]]$stage1_X1_fstat)) && 
        (results[[j]]$stage1_X1_fstat >= 10)&&(results[[j]]$stage1_X2_fstat >= 10)){
        lq_global_pval = pf(results[[j]]$stage2_lq$fstatistic[1],
                     results[[j]]$stage2_lq$fstatistic[2],
                     results[[j]]$stage2_lq$fstatistic[3],lower.tail=F)
        lq_nonlinear_pval = results[[j]]$stage2_lq$coefficients[3,4]
        lq_intercept = results[[j]]$stage2_lq$coefficients[1,1]
        lq_effect_X1 = results[[j]]$stage2_lq$coefficients[2,1]
        lq_effect_X2 = results[[j]]$stage2_lq$coefficients[3,1]
        lq_results = rbind(lq_results,
                           c(results[[j]]$gene,results[[j]]$chr,
                             results[[j]]$rsq_stage1_X1,
                             results[[j]]$adjrsq_stage1_X1,
                             results[[j]]$rsq_stage1_X2,
                             results[[j]]$adjrsq_stage1_X2,
                             lq_effect_X1,
                             lq_effect_X2,
                             lq_intercept,
                             lq_global_pval,
                             lq_nonlinear_pval,
                             results[[j]]$stage1_X1_fstat,
                             results[[j]]$stage1_X1_pval,
                             results[[j]]$stage1_X2_fstat,
                             results[[j]]$stage1_X2_pval,
                             results[[j]]$num_snps_X1,
                             results[[j]]$num_snps_X2))
    }
    
    if((results[[j]]$num_snps_X2 >= 2)&&
       !is.null(results[[j]]$stage1_X2_fstat)&&
       (results[[j]]$stage1_X2_fstat >=10)){
      q_pval = pf(results[[j]]$stage2_q$fstatistic[1],
                   results[[j]]$stage2_q$fstatistic[2],
                   results[[j]]$stage2_q$fstatistic[3],lower.tail=F)
      q_intercept = results[[j]]$stage2_q$coefficients[1,1]
      q_effect = results[[j]]$stage2_q$coefficients[2,1]
      q_results = rbind(q_results,
                         c(results[[j]]$gene,results[[j]]$chr,
                           results[[j]]$rsq_stage1_X2,
                           results[[j]]$adjrsq_stage1_X2,
                           q_effect,
                           q_intercept,
                           q_pval,
                           results[[j]]$stage1_X2_fstat,
                           results[[j]]$stage1_X2_pval,
                           results[[j]]$num_snps_X2))
    }
    
    if (!is.null(results[[j]]$stage2_poly_lq) && 
        !is.null(results[[j]]$stage2_poly_q) &&
        !is.null(results[[j]]$stage1_X1_fstat) &&
        results[[j]]$stage1_X1_fstat >= 10){
      poly_lq_effect_X1 = results[[j]]$stage2_poly_lq$coefficients[2,1]
      poly_lq_effect_X2 = results[[j]]$stage2_poly_lq$coefficients[3,1]
      poly_lq_intercept = results[[j]]$stage2_poly_lq$coefficients[1,1]
      poly_lq_global_pval = pf(results[[j]]$stage2_poly_lq$fstatistic[1],
                               results[[j]]$stage2_poly_lq$fstatistic[2],
                               results[[j]]$stage2_poly_lq$fstatistic[3],lower.tail=F)
      poly_lq_nonlinear_pval = results[[j]]$stage2_poly_lq$coefficients[3,4]
      
      poly_q_effect = results[[j]]$stage2_poly_lq$coefficients[2,1]
      poly_q_intercept = results[[j]]$stage2_poly_lq$coefficients[1,1]
      poly_q_pval = pf(results[[j]]$stage2_poly_q$fstatistic[1],
                       results[[j]]$stage2_poly_q$fstatistic[2],
                       results[[j]]$stage2_poly_q$fstatistic[3],lower.tail=F)
      
      poly_results = rbind(poly_results,
                        c(results[[j]]$gene,results[[j]]$chr,
                          results[[j]]$rsq_stage1_X1,
                          results[[j]]$adjrsq_stage1_X1,
                          poly_lq_effect_X1,
                          poly_lq_effect_X2,
                          poly_lq_intercept,
                          poly_lq_global_pval,
                          poly_lq_nonlinear_pval,
                          poly_q_effect,
                          poly_q_intercept,
                          poly_q_pval,
                          results[[j]]$stage1_X1_fstat,
                          results[[j]]$stage1_X1_pval,
                          results[[j]]$num_snps_X1))
    }
  
    if (!is.null(results[[j]]$stage2_X12_lq) && 
        !is.null(results[[j]]$stage1_X12_fstat) &&
        results[[j]]$stage1_X12_fstat >= 10 &&
        !is.null(results[[j]]$stage1_X1_fstat) &&
        results[[j]]$stage1_X1_fstat >= 10){
        X12_lq_effect_X1 = results[[j]]$stage2_X12_lq$coefficients[2,1]
        X12_lq_effect_X2 = results[[j]]$stage2_X12_lq$coefficients[3,1]
        X12_lq_intercept = results[[j]]$stage2_X12_lq$coefficients[1,1]
        X12_lq_global_pval = pf(results[[j]]$stage2_X12_lq$fstatistic[1],
                                 results[[j]]$stage2_X12_lq$fstatistic[2],
                                 results[[j]]$stage2_X12_lq$fstatistic[3],lower.tail=F)
        X12_lq_nonlinear_pval = results[[j]]$stage2_X12_lq$coefficients[3,4]
        
        X12_lq_results = rbind(X12_lq_results,
                             c(results[[j]]$gene,results[[j]]$chr,
                               # results[[j]]$rsq_stage1_X12,
                               # results[[j]]$adjrsq_stage1_X12,
                               X12_lq_effect_X1,
                               X12_lq_effect_X2,
                               X12_lq_intercept,
                               X12_lq_global_pval,
                               X12_lq_nonlinear_pval,
                               results[[j]]$stage1_X12_fstat[1],
                               results[[j]]$stage1_X12_pval,
                               results[[j]]$num_snps_X12))
    }
   
    if (!is.null(results[[j]]$stage2_X12_q) && 
        !is.null(results[[j]]$stage1_X12_fstat) &&
        results[[j]]$stage1_X12_fstat >= 10){
        X12_q_effect = results[[j]]$stage2_X12_q$coefficients[2,1]
        X12_q_intercept = results[[j]]$stage2_X12_q$coefficients[1,1]
        X12_q_pval = pf(results[[j]]$stage2_X12_q$fstatistic[1],
                        results[[j]]$stage2_X12_q$fstatistic[2],
                        results[[j]]$stage2_X12_q$fstatistic[3],lower.tail=F)
        
        X12_q_results = rbind(X12_q_results,
                              c(results[[j]]$gene,results[[j]]$chr,
                                # results[[j]]$rsq_stage1_X12,
                                # results[[j]]$adjrsq_stage1_X12,
                                X12_q_effect,
                                X12_q_intercept,
                                X12_q_pval,
                                results[[j]]$stage1_X12_fstat[1],
                                results[[j]]$stage1_X12_pval,
                                results[[j]]$num_snps_X12))
    }
  }
  lq_results = as.data.frame(lq_results)
  if (length(lq_results) > 0){
    colnames(lq_results) = lq_names
    fwrite(lq_results,"hdl_lq_results.txt", append = T,
           row.names = F, col.names = !file.exists('hdl_lq_results.txt'))
  }
  

  q_results = as.data.frame(q_results)
  if (length(q_results)>0){
    colnames(q_results) = q_names
    fwrite(q_results,"hdl_q_results.txt", append = T,
           row.names = F, col.names = !file.exists('hdl_q_results.txt'))
  }
 

  poly_results = as.data.frame(poly_results)
  if (length(poly_results) > 0){
    colnames(poly_results) = poly_names
    fwrite(poly_results,"hdl_poly_results.txt", append = T,
           row.names = F, col.names = !file.exists('hdl_poly_results.txt'))
  }
 
  X12_lq_results = as.data.frame(X12_lq_results)
  if (length(X12_lq_results) > 0){
    colnames(X12_lq_results) = X12_lq_names
    fwrite(X12_lq_results,"hdl_X12_lq_results.txt", append = T,
           row.names = F, col.names = !file.exists('hdl_X12_lq_results.txt'))
  }
  
  X12_q_results = as.data.frame(X12_q_results)
  if (length(X12_q_results) > 0){
    colnames(X12_q_results) = X12_q_names
    fwrite(X12_q_results,"hdl_X12_q_results.txt", append = T,
           row.names = F, col.names = !file.exists('hdl_X12_q_results.txt'))
  }
  
  rm(list=ls())
}




