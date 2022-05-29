# R/3.6.0
library(plink2R)
library(dplyr)
library(tidyr)
library(data.table)
library(TwoSampleMR)
library(ieugwasr)
library(glmnet)

source("~/deepRIV/UKB/code/allele_qc.R")

Mode <- function(x){
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

ModeImpute <- function(x){
  n = ncol(x)
  for(i in 1:n){
    na_idx = which(is.na(x[,i]))
    tmp_mode = Mode(x[-na_idx,1])
    x[na_idx,i] = tmp_mode
  }
  return(x)
}

white_unrelated_keep = fread('~/deepRIV/UKB/data/white_unrelated_keep_ind.txt',
                             header = F)
exclude_id = fread('~/deepRIV/UKB/data/exclude_ID.csv',
                   header = F)
white_unrelated_keep = white_unrelated_keep %>% filter(!(V1 %in% exclude_id$V1))
ukb_pheno_all = fread('~/deepRIV/UKB/data/ukb_hdl.txt')
ukb_pheno_all = na.omit(ukb_pheno_all)
keep_idx = sort(na.omit(match(white_unrelated_keep$V1, ukb_pheno_all$f.eid)))
ukb_pheno_all = ukb_pheno_all[keep_idx,]
rm(keep_idx)
covariates = fread('~/deepRIV/UKB/data/covariates.csv')

stage1 <- function(trait, 
                   save_models = FALSE, 
                   z_path = NULL, 
                   beta_path = NULL,
                   theta_path = NULL,
                   y_path = NULL){
  
  bed_dir = '/home/panwei/shared/UKBiobankIndiv/imputed/pgen/'
  ukb_pheno = ukb_pheno_all
  
  if (trait=='bmi'){
    id = 'ukb-a-248'
    exposure = fread('~/deepRIV/UKB/data/exposures.csv')[,c('f.eid','bmi')]
  }
  
  #exp_dat = extract_instruments(outcomes=id, clump=F)
  exp_dat = extract_instruments(outcomes=id)
  chrs = unique(exp_dat$chr.exposure)
  ukb_snp_all = NULL
  current_eids = ukb_pheno$f.eid
  
  for (i in 1:length(chrs)){
    print(i)
    chr = chrs[i]
    chr_idx = which(exp_dat$chr.exposure == chr)
    rsids = exp_dat$SNP[chr_idx]
    fwrite(data.frame(rsids), paste0('~/deepRIV/UKB/data/ForUKB/',trait,'_rs.txt'),
           row.names=F, col.names=F)
    prefix = paste0(bed_dir,'ukbb_chr',chr,'_1')
    ukb_out = paste0('~/deepRIV/UKB/data/ForUKB/',trait, '_chr', chr)
    ukb_plink_command = paste0('module load plink/2.00-alpha-091019; plink2 --pfile vzs ',prefix, ' --chr ',chr,
                               #' --keep /home/panwei/he000176/deepRIV/UKB/data/white_unrelated_keep_ind.txt',
                               ' --extract ','~/deepRIV/UKB/data/ForUKB/',trait,'_rs.txt',
                               #' --mind 0.1',
                               ' --make-bed',
                               ' --out ',ukb_out)
    ukb_plink_msg = system(ukb_plink_command,intern=TRUE)
    
    ukb_snp = tryCatch(read_plink(paste0("~/deepRIV/UKB/data/ForUKB/",trait,'_chr', chr)),
                       error=function(e){cat("ERROR :",
                                             conditionMessage(e),
                                             "\n")})
    if(is.null(ukb_snp))
    {
      next
    }
    remove_command = paste0("rm ~/deepRIV/UKB/data/ForUKB/",
                            trait,"_chr*")
    system(remove_command)
    remove_command = paste0("rm ~/deepRIV/UKB/data/ForUKB/",
                            trait,"_rs*")
    system(remove_command)
    
    row.names(ukb_snp$bed) = ukb_snp$fam$V1
    tmp_eids = intersect(current_eids, ukb_snp$fam$V1)
    if (length(tmp_eids) <= 35000){
      next
    }
    current_eids = tmp_eids
    ukb_snp$fam = ukb_snp$fam %>% filter(ukb_snp$fam$V1 %in% current_eids)
    ukb_snp$bed = as.data.frame(ukb_snp$bed) %>% filter(row.names(ukb_snp$bed) %in% current_eids)
    
    if(is.null(ukb_snp_all$bed)){
      ukb_snp_all$bed = ukb_snp$bed
      ukb_snp_all$bim = ukb_snp$bim
      # ukb_snp_all$fam = ukb_snp$fam
    }else{
      ukb_snp_all$bed = ukb_snp_all$bed %>% filter(row.names(ukb_snp_all$bed) %in% current_eids)
      ukb_snp$bed = ukb_snp$bed[match(rownames(ukb_snp_all$bed),rownames(ukb_snp$bed)),]
      print(c("are eids aligned? ", all(rownames(ukb_snp_all$bed)==rownames(ukb_snp$bed))))
      ukb_snp_all$bed = cbind(ukb_snp_all$bed, ukb_snp$bed)
      ukb_snp_all$bim = rbind(ukb_snp_all$bim, ukb_snp$bim)
    }
    
    print(length(current_eids))
  }
  
  ukb_snp_all$bed = ModeImpute(ukb_snp_all$bed)
  
  # flip allele
  tmp_exp_dat = exp_dat[,c(8,9,10,4,5)]
  colnames(ukb_snp_all$bim)[2] = 'SNP'
  tmp_exp_dat = merge(ukb_snp_all$bim, tmp_exp_dat, by = 'SNP', sort = F)
  remove_flip = allele.qc(tmp_exp_dat$V5,tmp_exp_dat$V6,
                          tmp_exp_dat$effect_allele.exposure,tmp_exp_dat$other_allele.exposure)
  tmp_exp_dat$beta.exposure[remove_flip$flip] = -1 * tmp_exp_dat$beta.exposure[remove_flip$flip]
  
  keep_ukb_indiv = intersect(ukb_pheno$f.eid,row.names(ukb_snp_all$bed))
  ukb_pheno = ukb_pheno %>% filter(f.eid %in% keep_ukb_indiv)
  ukb_snp_all$bed = ukb_snp_all$bed[which(rownames(ukb_snp_all$bed)%in%keep_ukb_indiv),,drop=FALSE]
  ukb_pheno = ukb_pheno[match(rownames(ukb_snp_all$bed),ukb_pheno$f.eid),]
  exposure = exposure[match(rownames(ukb_snp_all$bed),exposure$f.eid),]
  covariates = covariates[match(rownames(ukb_snp_all$bed),covariates$f.eid),]
  
  # check if f.eid's are aligned
  print(ukb_pheno$f.eid==exposure$f.eid && ukb_pheno$f.eid==covariates$f.eid
        && ukb_pheno$f.eid==rownames(ukb_snp_all$bed))
  
  exposure[,2] = scale(exposure[,2])
  ukb_pheno[,2] = scale(ukb_pheno[,2])
  ukb_snp_all$bed = scale(ukb_snp_all$bed)
  
  # regress out covariates
  exposure$bmi[is.na(exposure$bmi)] = mean(exposure$bmi, na.rm = T)
  lm_RegCov = lm(exposure.bmi ~.,data = data.frame(exposure$bmi,covariates[,-c(1)]))
  X1 = summary(lm_RegCov)$residuals
  X1 = scale(X1)
  
  # lasso
  lambdas = c(0, exp(seq(-5,5, length.out = 9)))
  lasso_stage1 = cv.glmnet(ukb_snp_all$bed, X1, lambda = lambdas,
                           nfolds = 10, alpha=1)
  lasso_stage1_X2 = cv.glmnet(ukb_snp_all$bed, X1^2, lambda = lambdas,
                              nfolds = 10, alpha=1)
  
  X_hat = predict(lasso_stage1, newx = ukb_snp_all$bed, s = 'lambda.min')
  X2_hat = predict(lasso_stage1_X2, newx = ukb_snp_all$bed, s = 'lambda.min')
  y_ukb = as.numeric(unlist(ukb_pheno[,2]))
  stage2_test = lm(y_ukb~X_hat)
  stage2_test = summary(stage2_test)
  
  stage2_polyq_test = lm(y_ukb~X_hat+I(X_hat^2))
  stage2_polyq_test = summary(stage2_polyq_test)
  
  stage2_lq_test = lm(y_ukb~X_hat+X2_hat)
  stage2_lq_test = summary(stage2_lq_test)
  
  FINAL_RESULT = list(Xhat_ukb_total = X_hat,
                      y_ukb_total = y_ukb,
                      stage2_pvalue = stage2_test$coefficients[2,4],
                      stage2_tvalue = stage2_test$coefficients[2,3],
                      stage2_polyq_pvalue = stage2_polyq_test$coefficients[3,4],
                      stage2_polyq_tvalue = stage2_polyq_test$coefficients[3,3],
                      stage2_lq_global = pf(stage2_lq_test$fstatistic[1], stage2_lq_test$fstatistic[2],
                                            stage2_lq_test$fstatistic[3], lower.tail = F),
                      stage2_lq_nonlinear = stage2_lq_test$coefficients[3,4],
                      num_snps = ncol(ukb_snp_all$bed),
                      sample_size = nrow(ukb_snp_all$bed)
  )
          
}



