for (i in 1:100){
  setwd("~/deepRIV/UKB/results/TWAS_LQ_Q/ldl")
  library(data.table)
  lq_names = c('gene_name','chr',
               'rsq_stage1_X1',
               'adjrsq_stage1_X1',
               'rsq_stage1_X2',
               'adjrsq_stage1_X2',
               'lq_effect_X1',
               'lq_effect_X2',
               'lq_intercept',
               'lq_pval',
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
  file_name = paste0('ldl_',i,'.RData')
  load(file_name)
  n = length(results)
  lq_results = NULL
  q_results = NULL
  for (j in 1:n){
    print(paste('i, j:',i,j))
    if (is.null(results[[j]])) next()
    if ((results[[j]]$num_snps_X2 >= 2)&&(results[[j]]$num_snps_X1 >= 2)&&
        (results[[j]]$stage1_X1_fstat >=10)&&(results[[j]]$stage1_X2_fstat >=10)){
        lq_pval = pf(results[[j]]$stage2_lq$fstatistic[1],
                     results[[j]]$stage2_lq$fstatistic[2],
                     results[[j]]$stage2_lq$fstatistic[3],lower.tail=F)
        lq_intercept = results[[j]]$stage2_lq$coefficients[1,1]
        lq_effect_X1 = results[[j]]$stage2_lq$coefficients[2,1]
        lq_effect_X2 = results[[j]]$stage2_lq$coefficients[3,1]
        lq_results  = rbind(lq_results,
                           c(results[[j]]$gene,results[[j]]$chr,
                             results[[j]]$rsq_stage1_X1,
                             results[[j]]$adjrsq_stage1_X1,
                             results[[j]]$rsq_stage1_X2,
                             results[[j]]$adjrsq_stage1_X2,
                             lq_effect_X1,
                             lq_effect_X2,
                             lq_intercept,
                             lq_pval,
                             results[[j]]$stage1_X1_fstat,
                             results[[j]]$stage1_X1_pval,
                             results[[j]]$stage1_X2_fstat,
                             results[[j]]$stage1_X2_pval,
                             results[[j]]$num_snps_X1,
                             results[[j]]$num_snps_X2))
    }else if((results[[j]]$num_snps_X2 >= 2)&&(results[[j]]$stage1_X2_fstat >=10)){
      q_pval = pf(results[[j]]$stage2_q$fstatistic[1],
                   results[[j]]$stage2_q$fstatistic[2],
                   results[[j]]$stage2_q$fstatistic[3],lower.tail=F)
      q_intercept = results[[j]]$stage2_q$coefficients[1,1]
      q_effect = results[[j]]$stage2_q$coefficients[2,1]
      q_results  = rbind(q_results,
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
  }
  lq_results = as.data.frame(lq_results)
  if (ncol(lq_results)==length(lq_names)){
    colnames(lq_results) = lq_names
    fwrite(lq_results,"ldl_lq_results.txt", append = T,
           row.names = F, col.names = !file.exists('ldl_lq_results.txt'))
  }
  
  q_results = as.data.frame(q_results)
  if (ncol(q_results)==length(q_names)){
    colnames(q_results) = q_names
    fwrite(q_results,"ldl_q_results.txt", append = T,
           row.names = F, col.names = !file.exists('ldl_q_results.txt'))
    
  }
  
  rm(list=ls())
}



