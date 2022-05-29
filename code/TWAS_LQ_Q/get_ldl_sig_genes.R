library(data.table)
setwd('~/deepRIV/UKB/code/TWAS_LQ_Q/')
source('~/deepRIV/UKB/code/TWAS_LQ_Q/ldl_sig_genes.R')

args = commandArgs(trailingOnly = T)
file_name = args[1]
gene_ind = fread(file_name, header = F)$V1

print(gene_ind)

results = list()

k=1
for (i in gene_ind){
  print(paste0('ldl: ',i))
  FINAL_RESULT = stage1(i)
  results[[k]] = FINAL_RESULT
  k = k + 1
}

out_path = paste0('/home/panwei/he000176/deepRIV/UKB/results/TWAS_LQ_Q/ldl/ldl_lq_sig_genes.RData')
save(results, file = out_path)



