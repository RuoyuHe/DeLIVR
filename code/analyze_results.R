library(dplyr)
library(tidyr)
library(data.table)
setwd("~/deepRIV/UKB/results/hdl/test40p_unrelated2/")

hdl_results = fread("hdl_combined_results.csv", header = T)

gene_gtf = rtracklayer::import('/home/panwei/shared/GTEx-v8/gencode.v26.GRCh38.genes.gtf')
gene_gtf = as.data.frame(gene_gtf) %>% select(gene_id,gene_name) %>% distinct()
hdl_results = merge(hdl_results, gene_gtf, by.x='gene_names', by.y='gene_id',
                    all.x=T, all.y=F, sort=F)
GENE_EXP = fread("/home/panwei/he000176/deepRIV/UKB/data/GTEx/Whole_Blood.v8.normalized_expression.bed.gz")
colnames(GENE_EXP)[1] = 'chr'
GENE_EXP$chr = as.numeric(gsub("[^0-9.-]", "", GENE_EXP$chr))
chr_idx = na.omit(match(hdl_results$gene_names, GENE_EXP$gene_id))
hdl_results$gene_ind = chr_idx
hdl_results$chr = GENE_EXP$chr[chr_idx]
hdl_results$bp = 1:nrow(hdl_results)
fwrite(hdl_results, "hdl_combined_results.csv", row.names = F)

    
library(qqman)
library(data.table)

setwd("~/Dropbox/UMN/Research_Projects/deepRIV/UKB/results/hdl")

hdl_results = fread("hdl_combined_results.csv")
# hdl_results$bp = 1:nrow(hdl_results)
# hdl_results$chr = as.numeric(gsub("[^0-9.-]", "", hdl_results$chr))
# hdl_results$IVa_pval_nonlinear = hdl_results$IVa_pval_nonlinear + runif(nrow(hdl_results))*1e-15
# hdl_results[,c(2,6)] = hdl_results[,c(2,6)] + 1e-15
# hdl_results$IVa_pval_global = jitter(hdl_results$IVa_pval_global, factor = 1e-10)

alpha = 0.05
adj_alpha = 0.05/nrow(hdl_results)
# valid_idx = which(hdl_results$adj_rsq >= 0.01)

manhattan(hdl_results, chr = "chr", bp = "bp", 
          p = "IVa_pval_global", snp = "gene_name", 
          suggestiveline = F, genomewideline = -log10(adj_alpha), 
          ylim = c(0,18), main = "DeLIVER Global Test", 
          annotatePval = adj_alpha)

manhattan(hdl_results, chr = "chr", bp = "bp", 
          p = "IVa_pval_nonlinear", snp = "gene_name", 
          suggestiveline = F, genomewideline = -log10(adj_alpha), 
          ylim = c(0,18), main = "DeLIVER Nonlinearity Test", 
          annotatePval = adj_alpha)

manhattan(hdl_results, chr = "chr", bp = "bp", 
          p = "LR_pval", snp = "gene_name", 
          suggestiveline = F, genomewideline = -log10(adj_alpha), 
          ylim = c(0,18), main = "Linear Regression", 
          annotatePval = adj_alpha)
