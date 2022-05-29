library(data.table)
source("~/deepRIV/UKB/code/TWAS_LQ_Q/ldl_twaslq_q.R")

args = commandArgs(trailingOnly = T)
file_num = as.integer(args[1])

num_genes = 20315
num_files = 60
file_size = ceiling(num_genes/num_files)
begin = (file_num-1)*file_size + 1
finish = min(file_num*file_size, num_genes)

results = list()
k=1
for (i in begin:finish){
  print(paste0('ldl: ',i))
  FINAL_RESULT = stage1(i)
  results[[k]] = FINAL_RESULT
  k = k + 1
}

out_path = paste0('/home/panwei/he000176/deepRIV/UKB/results/TWAS_LQ_Q/ldl/ldl_',
                  file_num,'.RData')
save(results, file = out_path)

