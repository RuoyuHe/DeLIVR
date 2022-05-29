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
covariates = fread('~/deepRIV/UKB/data/covariates.csv')

stage1 <- function(trait, 
                   save_models = FALSE, 
                   z_path = NULL, 
                   beta_path = NULL,
                   theta_path = NULL,
                   y_path = NULL){
  
  bed_dir = '/home/panwei/shared/UKBiobankIndiv/imputed/pgen/'
  
  if (trait=='bmi'){
    id = 'ukb-a-248'
    exposure = fread('~/deepRIV/UKB/data/exposures.csv')[,c('f.eid','bmi')]
  }
  
  # read outcomes
  hdl = fread('~/deepRIV/UKB/data/ukb_hdl.txt')
  hdl = na.omit(hdl)
  ldl = fread('~/deepRIV/UKB/data/ukb_ldl.txt')
  ldl = na.omit(ldl)
  keep_idx = sort(na.omit(match(white_unrelated_keep$V1, hdl$f.eid)))
  hdl = hdl[keep_idx,]
  rm(keep_idx)
  outcomes = fread('~/deepRIV/UKB/data/outcomes.csv')
  outcomes = merge(outcomes, hdl, by = 'f.eid')
  outcomes = merge(outcomes, ldl, by = 'f.eid')
  colnames(outcomes)[5] = 'hdl'
  outcomes = as.data.frame(na.omit(outcomes))
  #exp_dat = extract_instruments(outcomes=id, clump=F)
  exp_dat = extract_instruments(outcomes=id)
  chrs = unique(exp_dat$chr.exposure)
  ukb_snp_all = NULL
  current_eids = outcomes$f.eid
  
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
  
  keep_ukb_indiv = intersect(outcomes$f.eid,row.names(ukb_snp_all$bed))
  outcomes = outcomes %>% filter(f.eid %in% keep_ukb_indiv)
  ukb_snp_all$bed = ukb_snp_all$bed[which(rownames(ukb_snp_all$bed)%in%keep_ukb_indiv),,drop=FALSE]
  outcomes = outcomes[match(rownames(ukb_snp_all$bed),outcomes$f.eid),]
  exposure = exposure[match(rownames(ukb_snp_all$bed),exposure$f.eid),]
  covariates = covariates[match(rownames(ukb_snp_all$bed),covariates$f.eid),]
  
  # check if f.eid's are aligned
  print(outcomes$f.eid==exposure$f.eid && outcomes$f.eid==covariates$f.eid
        && outcomes$f.eid==rownames(ukb_snp_all$bed))
  
  exposure$bmi[is.na(exposure$bmi)] = mean(exposure$bmi, na.rm = T)
  exposure[,2] = scale(exposure[,2])
  outcomes[c(2:ncol(outcomes))] = scale(outcomes[c(2:ncol(outcomes))])
  ukb_snp_all$bed = scale(ukb_snp_all$bed)
  
  fwrite(outcomes,paste0('~/deepRIV/UKB/data/MR_tmp/normalized_outcomes_',trait,'.csv'))
  fwrite(exposure,paste0('~/deepRIV/UKB/data/MR_tmp/normalized_',trait,'.csv'))
  
  # regress out covariates
  lm_RegCov = lm(exposure.bmi ~.,data = data.frame(exposure$bmi,covariates[,-c(1)]))
  X1 = summary(lm_RegCov)$residuals
  X1 = scale(X1)
  X2 = X1^2

  fwrite(as.data.frame(X1),
         paste0('~/deepRIV/UKB/data/MR_tmp/normalized_',trait,'_adjusted_cov.csv'))
  
  if (F){
    # lasso
    lambdas = c(0, exp(seq(-5,5, length.out = 9)))
    lasso_stage1 = cv.glmnet(ukb_snp_all$bed, X1, lambda = lambdas,
                             nfolds = 10, alpha=1)
    lasso_stage1_X2 = cv.glmnet(ukb_snp_all$bed, X1^2, lambda = lambdas,
                                nfolds = 10, alpha=1)
    
    X_hat = predict(lasso_stage1, newx = ukb_snp_all$bed, s = 'lambda.min')
    fwrite(as.data.frame(X_hat),
           paste0('~/deepRIV/UKB/data/MR_tmp/lasso_',trait,'_adjusted_cov.csv'))
  }
  
  # LR
  lm_stage1 = lm(X1 ~ ., data.frame(X1,ukb_snp_all$bed))
  X_hat_lm = lm_stage1$fitted.values
  fwrite(as.data.frame(X_hat_lm),
         paste0('~/deepRIV/UKB/data/MR_tmp/lm_',trait,'_adjusted_cov.csv'))
  
  lm_stage1_X2 = lm(X2 ~ ., data.frame(X2,ukb_snp_all$bed))
  X2_hat_lm = lm_stage1_X2$fitted.values
  fwrite(as.data.frame(X_hat_lm),
         paste0('~/deepRIV/UKB/data/MR_tmp/lm_X2_',trait,'_adjusted_cov.csv'))
}

stage1('bmi')

