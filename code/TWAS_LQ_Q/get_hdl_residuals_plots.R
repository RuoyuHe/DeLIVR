library(plink2R)
library(dplyr)
library(tidyr)
library(data.table)

source("~/deepRIV/UKB/code/allele_qc.R")

##  # Covariates
COV5TABLE = fread("~/deepRIV/UKB/data/GTEx/Whole_Blood.v8.covariates.txt",header = T)
### Gene Expression
GENE_EXP = fread("~/deepRIV/UKB/data/GTEx/Whole_Blood.v8.normalized_expression.bed.gz")
colnames(GENE_EXP)[1] = 'chr'
whole_blood_buddy = '/home/panwei/he000176/deepRIV/UKB/data/Whole_Blood_buddy.txt'


### ukb phenotype
white_unrelated_keep = fread('~/deepRIV/UKB/data/white_unrelated_keep_ind.txt',
                             header = F)
ukb_pheno_all = fread('~/deepRIV/UKB/data/ukb_hdl.txt')
ukb_pheno_all = na.omit(ukb_pheno_all)
keep_idx = sort(na.omit(match(white_unrelated_keep$V1, ukb_pheno_all$f.eid)))
ukb_pheno_all = ukb_pheno_all[keep_idx,]
rm(keep_idx)

repeat_genes = fread("~/deepRIV/UKB/code/hdl/combine_p/deepIV/hdl_deepiv_genes.csv")
hdl = fread("/home/panwei/he000176/deepRIV/UKB/results/hdl/test40p_unrelated2/hdl_combined_results.csv")
repeat_genes$gene_ind = hdl %>% filter(hdl$gene_name %in% repeat_genes$x) %>% select(gene_ind)

path = '~/deepRIV/UKB/results/TWAS_LQ_Q/hdl/'

stage1 <- function(gene_ind, 
                   save_models = FALSE, 
                   z_path = NULL, 
                   beta_path = NULL,
                   theta_path = NULL,
                   y_path = NULL){
  # get covariates
  cov5table = COV5TABLE
  # get gene expression
  gene_exp = GENE_EXP[gene_ind,]
  chr = gene_exp$chr[1]
  chr = as.numeric(gsub("[^0-9.-]", "", chr))
  
  ###UKB
  bed_dir = '/home/panwei/shared/UKBiobankIndiv/genetic_data/ukb_cal_chr'
  prefix = paste0(bed_dir,chr,'_v2')
  
  # start analysis ----------------------------------------------------------
  ukb_pheno = ukb_pheno_all
  gene_name = as.character(gene_exp$gene_id[1])
  start = floor(gene_exp$start[1]/1000) - 100
  end = ceiling(gene_exp$end[1]/1000) + 100
  
  plink_command = paste0("module load plink; \n","plink --bfile /home/panwei/he000176/deepRIV/UKB/data/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_MAF01",
                         " --chr ",chr," --from-kb ",
                         start," --to-kb ",end,
                         " --keep ",whole_blood_buddy,
                         " --geno 0 --maf 0.05 --hwe 0.001  ",
                         " --make-bed --out ~/deepRIV/UKB/data/For_Individual_Genes/",
                         gene_name,sep="")
  plink_msg=system(plink_command,intern=TRUE)
  
  snp = tryCatch(read_plink(paste("~/deepRIV/UKB/data/For_Individual_Genes/",gene_name,sep="")),
                 error=function(e){cat("ERROR :",
                                       conditionMessage(e),
                                       "\n")})
  if(is.null(snp))
  {
    FINAL_RESULT = NULL
    return(FINAL_RESULT)
  }
  remove_command = paste("rm ~/deepRIV/UKB/data/For_Individual_Genes/",
                         gene_name,".*",sep="")
  system(remove_command)
  
  # match individuals with the covariates
  # cov5table = na.omit(cov5table[match(snp$fam$V2,cov5table$ID),])
  # snp$bed = na.omit(snp$bed[match(cov5table$ID,snp$fam$V2),])
  # snp$fam = na.omit(snp$fam[match(cov5table$ID,snp$fam$V2),])
  
  write.table(snp$bim$V2,paste0('~/deepRIV/UKB/data/ForUKB/',gene_name,'_rs.txt'),row.names=F,quote=F,col.names=F)
  
  grep_command = paste0("zgrep -F -m ",length(snp$bim$V2)," -f ~/deepRIV/UKB/data/ForUKB/",gene_name,"_rs.txt /home/panwei/he000176/deepRIV/UKB/data/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt  > ~/deepRIV/UKB/data/ForUKB/",gene_name,"_rs1.txt")
  system(grep_command)
  gtex_dic = read.table(paste0('~/deepRIV/UKB/data/ForUKB/',gene_name,'_rs1.txt'))
  write.table(gtex_dic$V7,paste0('~/deepRIV/UKB/data/ForUKB/',gene_name,'_rs.txt'),row.names=F,quote=F,col.names=F)
  
  ukb_out = paste0('~/deepRIV/UKB/data/ForUKB/',gene_name)
  ukb_plink_command = paste0('module load plink; plink --bfile ',prefix,' --chr ',chr,
                             ' --keep /home/panwei/he000176/deepRIV/UKB/data/white_unrelated_keep_ind.txt',
                             ' --extract ','~/deepRIV/UKB/data/ForUKB/',gene_name,'_rs.txt',
                             ' --mind 0 --make-bed',
                             ' --out ',ukb_out)
  ukb_plink_msg = system(ukb_plink_command,intern=TRUE)
  
  ukb_snp = tryCatch(read_plink(paste("~/deepRIV/UKB/data/ForUKB/",gene_name,sep="")),
                     error=function(e){cat("ERROR :",
                                           conditionMessage(e),
                                           "\n")})
  if(is.null(ukb_snp))
  {
    FINAL_RESULT = NULL
    return(FINAL_RESULT)
  }
  remove_command = paste0("rm ~/deepRIV/UKB/data/ForUKB/",
                          gene_name,".*")
  system(remove_command)
  remove_command = paste0("rm ~/deepRIV/UKB/data/ForUKB/",
                          gene_name,"_rs*")
  system(remove_command)
  
  ### overlap between GTEx SNP and UKB
  #colnames(snp$bed) = snp$bim$V4
  rownames(ukb_snp$bed) = ukb_snp$fam[,2]
  
  #colnames(ukb_snp$bim)[4:6] = c("Position","A1","A2")
  gtex_bim_tmp = merge(snp$bim,gtex_dic,by.x='V2',by.y='V1',sort=F)
  snp$bim[,2] = gtex_bim_tmp$V7
  colnames(snp$bed) = snp$bim[,2]
  rownames(snp$bed) = snp$fam[,2]
  snp_bim_ukb = merge(snp$bim,ukb_snp$bim,by=c("V1","V2"))
  remove_flip = allele.qc(snp_bim_ukb$V5.x,snp_bim_ukb$V6.x,
                          snp_bim_ukb$V5.y,snp_bim_ukb$V6.y)
  flip_snp = snp_bim_ukb$V2[remove_flip$flip]
  snp_bim_ukb = snp_bim_ukb[remove_flip$keep,]
  
  ukb_snp_bed = ukb_snp$bed[,which(colnames(ukb_snp$bed) %in% snp_bim_ukb$V2),drop=FALSE] 
  ukb_snp_bed[,which(colnames(ukb_snp_bed) %in% flip_snp)] = 2 - ukb_snp_bed[,which(colnames(ukb_snp_bed) %in% flip_snp)]
  
  snp_bed_ind = NULL
  snp_bed_colname = colnames(snp$bed)
  for(i in 1:nrow(snp_bim_ukb))
  {
    snp_bed_ind = c(snp_bed_ind,which(snp_bed_colname == colnames(ukb_snp_bed)[i]))
  }
  SNP_BED = snp$bed[,snp_bed_ind,drop=FALSE]
  cov5table = cov5table[match(rownames(SNP_BED),cov5table$ID),]
  
  keep_ukb_indiv = intersect(ukb_pheno$f.eid,ukb_snp$fam[,2])
  ukb_pheno = ukb_pheno %>% filter(f.eid %in% keep_ukb_indiv)
  ukb_snp_bed = ukb_snp_bed[which(rownames(ukb_snp_bed)%in%keep_ukb_indiv),,drop=FALSE]
  ukb_pheno = ukb_pheno[match(rownames(ukb_snp_bed),ukb_pheno$f.eid),]
  
  if(nrow(ukb_pheno)< 3){
    FINAL_RESULT = NULL
    return(FINAL_RESULT)
  }
  if(length(snp_bed_ind) < 2)
  {
    FINAL_RESULT = NULL
    return(FINAL_RESULT)
  }
  ### prune
  cor_cutoff = 0.8
  cor_bed = abs(cor(SNP_BED))
  cor_bed = (cor_bed < cor_cutoff)^2
  diag(cor_bed) = 1
  i = 1
  while(i < nrow(cor_bed) )
  {
    ind = which(cor_bed[i,] == 1)
    cor_bed = as.matrix(cor_bed[ind,ind])
    i = i + 1
  }
  if(nrow(cor_bed) == 1)
  {
    FINAL_RESULT = NULL
    return(FINAL_RESULT)
  }
  ind = which(is.element(colnames(SNP_BED),colnames(cor_bed)))
  SNP_BED = SNP_BED[,ind,drop=F]
  ukb_snp_bed = ukb_snp_bed[,ind]
  snp_bim_ukb = snp_bim_ukb[ind,,drop=F]
  
  
  ### regress gene exp on covariates
  X = as.matrix(gene_exp[1,-(1:4)])
  X = X[,match(rownames(SNP_BED),colnames(X)),drop=F]
  X1 = scale(as.numeric(X))
  X2 = X1^2
  
  lm_RegCov = lm(X1 ~.,data = data.frame(X1,cov5table[,-c(1)]))
  X1 = lm_RegCov$residuals
  X1 = scale(X1)
  
  lm_RegCov = lm(X2 ~.,data = data.frame(X2,cov5table[,-c(1,2)]))
  X2 = lm_RegCov$residuals
  X2 = scale(X2)
  
  if(ncol(SNP_BED) > 50)
  {
    ind1 = order(abs(cor(X1,SNP_BED)),decreasing = T)[1:50]
    ind2 = order(abs(cor(X2,SNP_BED)),decreasing = T)[1:50]
    ind = union(ind1,ind2)
    
    ind = sort(ind)
    SNP_BED = SNP_BED[,ind]
    ukb_snp_bed = ukb_snp_bed[,ind]
    snp_bim_ukb = snp_bim_ukb[ind,]
  }
  
  SNP_BED = scale(SNP_BED)
  ukb_snp_bed = scale(ukb_snp_bed)
  
  ###
  lm_stage1_X1 = lm(X1 ~ SNP_BED)
  hatbetaX1 = lm_stage1_X1$coefficients[-1]
  na_ind = which(!is.na(hatbetaX1))
  check_ukb = apply(ukb_snp_bed,2,sd)
  na_ind_ukb = which(!is.na(check_ukb))
  na_ind = intersect(na_ind,na_ind_ukb)
  SNP_BED = SNP_BED[,na_ind,drop=FALSE]
  ukb_snp_bed = ukb_snp_bed[,na_ind,drop=FALSE]
  
  ## get E(X1 | Z)
  lm_stage1_X1 = 
    step(lm(X1~.,data = data.frame(X1,SNP_BED)),direction = "backward",trace=FALSE)
  AIC_stage1_X1 = AIC(lm_stage1_X1)
  BIC_stage1_X1 = BIC(lm_stage1_X1)
  lm_stage1_X1 = summary(lm_stage1_X1)
  stage1_sigma = lm_stage1_X1$sigma
  rsq_stage1_X1 = lm_stage1_X1$r.squared
  adjrsq_stage1_X1 = lm_stage1_X1$adj.r.squared
  se_betaX1 = hatbetaX1 = rep(0,ncol(SNP_BED))
  coef_X1 = lm_stage1_X1$coefficients
  name_coef_X1 = substr(rownames(coef_X1),1,30)
  name_SNP_BED = colnames(SNP_BED)
  for(beta_ind in 1:nrow(coef_X1))
  {
    ii = which(name_SNP_BED == name_coef_X1[beta_ind])
    hatbetaX1[ii] = coef_X1[beta_ind,1]
    se_betaX1[ii] = coef_X1[beta_ind,2]
  }
  
  ## get E(X2 | Z)
  lm_stage1_X2 = 
    step(lm(X2~.,data = data.frame(X2,SNP_BED)),direction = "backward",trace=FALSE)
  AIC_stage1_X2 = AIC(lm_stage1_X2)
  BIC_stage1_X2 = BIC(lm_stage1_X2)
  lm_stage1_X2 = summary(lm_stage1_X2)
  rsq_stage1_X2 = lm_stage1_X2$r.squared
  adjrsq_stage1_X2 = lm_stage1_X2$adj.r.squared
  se_betaX2 = hatbetaX2 = rep(0,ncol(SNP_BED))
  coef_X2 = lm_stage1_X2$coefficients
  name_coef_X2 = substr(rownames(coef_X2),1,30)
  name_SNP_BED = colnames(SNP_BED)
  for(beta_ind in 1:nrow(coef_X2))
  {
    ii = which(name_SNP_BED == name_coef_X2[beta_ind])
    hatbetaX2[ii] = coef_X2[beta_ind,1]
    se_betaX2[ii] = coef_X2[beta_ind,2]
  }
  
  ukb_pheno[,2] = scale(ukb_pheno[,2])
  y_ukb = as.numeric(unlist(ukb_pheno[,2]))
  
  num_snps_X1 = sum(hatbetaX1 != 0)
  num_snps_X2 = sum(hatbetaX2 != 0)
  
  Xhat_ukb = ukb_snp_bed %*% hatbetaX1 + coef_X1[1,1]
  X2hat_ukb = ukb_snp_bed %*% hatbetaX2 + coef_X2[1,1]
  stage2_l = lm(y_ukb ~ Xhat_ukb)
  tmpl = summary(stage2_l)
  stage2_lq = lm(y_ukb ~ Xhat_ukb + X2hat_ukb)
  tmplq = summary(stage2_lq)
  lq_nonlinear_p = anova(stage2_l, stage2_lq)$`Pr(>F)`[-1]
  
  stage2_polyq = lm(y_ukb ~ Xhat_ukb + I(Xhat_ukb^2))
  tmppolyq = summary(stage2_polyq)
  
  qt = quantile(Xhat_ukb, c(0.025, 0.975))
  unq1 = unique(Xhat_ukb)
  unq1 = (as.data.frame(unq1) %>% filter(unq1 >= qt[1] & unq1 <= qt[2]))$V1
  unq = c()
  mean_res_lq = c()
  mean_res_polyq = c()
  for(i in 1:length(unq1)){
    idx = which(Xhat_ukb == unq1[i])
    # print(length(idx))
    if (length(idx) > 1000){
      unq = c(unq, unq1[i])
      mean_res_lq = c(mean_res_lq, mean(tmplq$residuals[idx]))
      mean_res_polyq = c(mean_res_polyq, mean(tmppolyq$residuals[idx]))
    }
  }
  unq = sort(unq, index.return = T)
  idx = unq$ix
  unq = unq$x
  mean_res_lq = mean_res_lq[idx]
  mean_res_polyq = mean_res_polyq[idx]
  
  gene_name = repeat_genes$x[which(repeat_genes$gene_ind==gene_ind)]
  pdf(paste0(path,'lq_', gene_name,'.pdf')) 
  plot(x = unq, y = mean_res_lq, pch = '.',
       xlab = "mu", ylab = "res")
  dev.off() 
  
  pdf(paste0(path,'polyq_', gene_name,'.pdf')) 
  plot(x = unq, y = mean_res_polyq, pch = '.',
       xlab = "mu", ylab = "res")
  dev.off() 
}

for (i in 1:nrow(repeat_genes)){
  stage1(repeat_genes$gene_ind[i])
}