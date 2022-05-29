library(plink2R)
library(dplyr)
library(tidyr)
library(data.table)
library(jsonlite)

setwd("~/deepRIV/UKB/data")
source("~/deepRIV/UKB/code/allele_qc.R")

chr = 1
chr0 = paste0('chr',chr)

### train, val, test
train_id = fread("hdl_train_id.txt", header = T)
val_id = fread("hdl_val_id.txt", header = T)
test_id = fread("hdl_test_id.txt", header = T)

### ukb phenotype
ukb_pheno_all = fread('~/deepRIV/UKB/data/ukb_hdl.txt')
ukb_pheno_all = na.omit(ukb_pheno_all)

###UKB
bed_dir = '/home/panwei/shared/UKBiobankIndiv/genetic_data/ukb_cal_chr'
prefix = paste0(bed_dir,chr,'_v2')

### Covariates
cov5table = fread("~/deepRIV/UKB/data/Whole_Blood.v8.covariates.txt",header = T)
### Gene Expression
gene_exp = fread("/home/panwei/shared/GTEx_v8/GTEx_Analysis_v8_eQTL_expression_matrices/Whole_Blood.v8.normalized_expression.bed.gz")
colnames(gene_exp)[1] = 'chr'
gene_exp = gene_exp %>% filter(chr==chr0)

cat(nrow(gene_exp),"\n")
real_data_result = list()

# start analysis ----------------------------------------------------------
for(gene_ind in 1:nrow(gene_exp)){
  cat(gene_ind,"\n")
  ukb_pheno = ukb_pheno_all
  gene_name = as.character(gene_exp$gene_id[gene_ind])
  start = floor(gene_exp$start[gene_ind]/1000) - 100
  end = ceiling(gene_exp$end[gene_ind]/1000) + 100
  
  plink_command = paste0("module load plink; \n","plink --bfile /home/panwei/shared/GTEx_v8/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_MAF01",
                         " --chr ",chr," --from-kb ",
                         start," --to-kb ",end,
                         " --geno 0 --maf 0.05 --hwe 0.001  ",
                         " --make-bed --out For_Individual_Genes/",
                         gene_name,sep="")
  plink_msg=system(plink_command,intern=TRUE)
  
  snp = tryCatch(read_plink(paste("For_Individual_Genes/",gene_name,sep="")),
                 error=function(e){cat("ERROR :",
                                       conditionMessage(e),
                                       "\n")})
  
  if(is.null(snp))
  {
    FINAL_RESULT = NULL
    # real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
    next()
  }
  remove_command = paste("rm For_Individual_Genes/",
                         gene_name,".*",sep="")
  system(remove_command)
  write.table(snp$bim$V2,paste0('ForUKB/',gene_name,'_rs.txt'),row.names=F,quote=F,col.names=F)
  
  grep_command = paste0("zgrep -F -m ",length(snp$bim$V2)," -f ForUKB/",gene_name,"_rs.txt /home/panwei/shared/GTEx_v8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz  > ForUKB/",gene_name,"_rs1.txt")
  system(grep_command)
  gtex_dic = read.table(paste0('ForUKB/',gene_name,'_rs1.txt'))
  write.table(gtex_dic$V7,paste0('ForUKB/',gene_name,'_rs.txt'),row.names=F,quote=F,col.names=F)
  
  ukb_out = paste0('ForUKB/',gene_name)
  ukb_plink_command = paste0('module load plink; plink --bfile ',prefix,' --chr ',chr,
                             ' --extract ','ForUKB/',gene_name,'_rs.txt',
                             ' --mind 0 --make-bed',
                             ' --out ',ukb_out)
  ukb_plink_msg = system(ukb_plink_command,intern=TRUE)
  
  ukb_snp = tryCatch(read_plink(paste("ForUKB/",gene_name,sep="")),
                     error=function(e){cat("ERROR :",
                                           conditionMessage(e),
                                           "\n")})
  if(is.null(ukb_snp))
  {
    FINAL_RESULT = NULL
    # real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
    next()
  }
  remove_command = paste("rm ForUKB/",
                         gene_name,".*",sep="")
  system(remove_command)
  remove_command = paste("rm ForUKB/",
                         gene_name,"_rs*",sep="")
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
  
  ukb_snp_bed = ukb_snp$bed[,which(colnames(ukb_snp$bed)%in% snp_bim_ukb$V2),drop=FALSE] 
  ukb_snp_bed[,which(colnames(ukb_snp_bed) %in% flip_snp)] = 2 - ukb_snp_bed[,which(colnames(ukb_snp_bed) %in% flip_snp)]
  
  snp_bed_ind = NULL
  snp_bed_colname = colnames(snp$bed)
  for(i in 1:nrow(snp_bim_ukb))
  {
    snp_bed_ind = c(snp_bed_ind,which(snp_bed_colname == colnames(ukb_snp_bed)[i]))
  }
  SNP_BED = snp$bed[,snp_bed_ind,drop=FALSE]
  cov5table = na.omit(cov5table[match(rownames(SNP_BED),cov5table$ID),])
  SNP_BED = na.omit(SNP_BED[match(cov5table$ID,rownames(SNP_BED)),])
  
  
  keep_ukb_indiv = intersect(ukb_pheno$f.eid,ukb_snp$fam[,2])
  ukb_pheno = ukb_pheno %>% filter(f.eid %in% keep_ukb_indiv)
  ukb_snp_bed = ukb_snp_bed[which(rownames(ukb_snp_bed)%in%keep_ukb_indiv),,drop=FALSE]
  ukb_pheno = ukb_pheno[match(rownames(ukb_snp_bed),ukb_pheno$f.eid),]
  
  if(nrow(ukb_pheno)< 3){
    FINAL_RESULT = NULL
    # real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
    next()
  }
  
  if(length(snp_bed_ind) < 2)
  {
    FINAL_RESULT = NULL
    # real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
    next()
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
    # real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
    next()
  }
  ind = which(is.element(colnames(SNP_BED),colnames(cor_bed)))
  SNP_BED = SNP_BED[,ind,drop=F]
  ukb_snp_bed = ukb_snp_bed[,ind]
  snp_bim_ukb = snp_bim_ukb[ind,,drop=F]
  
  
  ### regress gene exp on covariates
  X = as.matrix(gene_exp[gene_ind,-(1:4)])
  X = X[,match(rownames(SNP_BED),colnames(X)),drop=F]
  # X = scale(as.numeric(X))
  # X = as.numeric(gene_exp[gene_ind,-(1:4)]) #A1,A3
  X1 = as.numeric(X)
  
  lm_RegCov = lm(X1 ~.,data = data.frame(X1,cov5table[,-c(1,2)]))
  X1 = lm_RegCov$residuals
  # X1 = scale(X1)
  if(ncol(SNP_BED) > 50)
  {
    ind = order(abs(cor(X1,SNP_BED)),decreasing = T)[1:50]
    ind = sort(ind)
    SNP_BED = SNP_BED[,ind]
    ukb_snp_bed = ukb_snp_bed[,ind]
    snp_bim_ukb = snp_bim_ukb[ind,]
  }
  # SNP_BED = scale(SNP_BED)
  # ukb_snp_bed = scale(ukb_snp_bed)
  
  
  ###
  lm_stage1_X1 = lm(X1 ~ SNP_BED)
  hatbetaX1 = lm_stage1_X1$coefficients
  na_ind = which(!is.na(hatbetaX1))
  check_ukb = apply(ukb_snp_bed,2,sd)
  na_ind_ukb = which(!is.na(check_ukb))
  na_ind = intersect(na_ind,na_ind_ukb)
  SNP_BED = SNP_BED[,na_ind,drop=FALSE]
  ukb_snp_bed = ukb_snp_bed[,na_ind,drop=FALSE]
  
  
  lm_stage1_X1 = 
    step(lm(X1~.,data = data.frame(X1,SNP_BED)),direction = "backward",trace=FALSE)
  AIC_stage1_X1 = AIC(lm_stage1_X1)
  BIC_stage1_X1 = BIC(lm_stage1_X1)
  lm_stage1_X1 = summary(lm_stage1_X1)
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
  
  if(sum(abs(hatbetaX1)) == 0 )
  {
    FINAL_RESULT = NULL
    # real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
    next()
  }
  
  ### UKB
  train_idx = sort(na.omit(match(train_id$f.eid, ukb_pheno$f.eid)))
  val_idx = sort(na.omit(match(val_id$f.eid, ukb_pheno$f.eid)))
  test_idx = sort(na.omit(match(test_id$f.eid, ukb_pheno$f.eid)))
  
  Xhat_ukb_train = ukb_snp_bed[train_idx,] %*% hatbetaX1 + coef_X1[1,1]
  Xhat_ukb_val = ukb_snp_bed[val_idx,] %*% hatbetaX1 + coef_X1[1,1]
  Xhat_ukb_test = ukb_snp_bed[test_idx,] %*% hatbetaX1 + coef_X1[1,1]
  y_ukb_train = as.numeric(unlist(ukb_pheno[train_idx,2]))
  y_ukb_val = as.numeric(unlist(ukb_pheno[val_idx,2]))
  y_ukb_test = as.numeric(unlist(ukb_pheno[test_idx,2]))
  
  ### Stage 2
  stage2_test = lm(y_ukb_test~Xhat_ukb_test)
  stage2_test = summary(stage2_test)
  
  FINAL_RESULT = list(gene=gene_name,
                      chr=chr,
                      Xhat_ukb_train = Xhat_ukb_train,
                      Xhat_ukb_val = Xhat_ukb_val,
                      Xhat_ukb_test = Xhat_ukb_test,
                      y_ukb_train = y_ukb_train,
                      y_ukb_val = y_ukb_val,
                      y_ukb_test = y_ukb_test,
                      rsq_stage1 = rsq_stage1_X1,
                      adjrsq_stage1 = adjrsq_stage1_X1,
                      stage2_pvalue = stage2_test$coefficients[2,4]
  )
  real_data_result[[gene_name]] = FINAL_RESULT
}
real_data_result = toJSON(real_data_result)
write(real_data_result, 
      paste0("~/deepRIV/UKB/results/stage1_results/hdl_chr",chr,".json"))
