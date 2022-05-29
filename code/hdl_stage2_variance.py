import sys
print(sys.version)
import pdb
import json
import pickle
import numpy as np
import pandas as pd
import tensorflow as tf
import tensorflow.keras as K
import statsmodels.api as sm
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from tensorflow.keras import Model
from tensorflow.keras.layers import Activation, Add, Input, Dense, ReLU, Reshape
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.optimizers import Adam
from sklearn.model_selection import StratifiedShuffleSplit
from tensorflow.keras.backend import expand_dims, squeeze, clear_session
from sklearn.model_selection import train_test_split
from wrappers import create_stage1, create_stage2, fit_stage2, tune_l2, stage2Tests
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt


# white_unrelated_keep = fread('~/deepRIV/UKB/data/white_unrelated_keep_ind.txt',
# 							 header = F)
# ukb_pheno_all = fread('~/deepRIV/UKB/data/ukb_hdl.txt')
# ukb_pheno_all = na.omit(ukb_pheno_all)
# keep_idx = sort(na.omit(match(white_unrelated_keep$V1, ukb_pheno_all$f.eid)))
# ukb_pheno_all = ukb_pheno_all[keep_idx,]
# rm(keep_idx)



# Defining the R script and loading the instance in Python
r = robjects.r
r['source']('hdl_stage1_FUN.R')

# Loading the function we have defined in R.
stage1_r = robjects.globalenv['stage1']

file_num = 1
gene_ind = 15037
# runs = len(sig_gene_ind)

runs = 30


gene_names = []
num_snps = np.zeros(runs)
IVa_pval_global = np.zeros(runs)
IVa_pval_nonlinear = np.zeros(runs)
IVa_tval_global = np.zeros(runs)
IVa_tval_nonlinear = np.zeros(runs)
LR_pval = np.zeros(runs)
LR_tval = np.zeros(runs)
rsq = np.zeros(runs)
adj_rsq = np.zeros(runs)
stage1_fstat = np.zeros(runs)
stage1_pval = np.zeros(runs)
stage1_sigma = np.zeros(runs)
nonlinear_coef = np.zeros(runs)
nonlinear_se = np.zeros(runs)
mu_max = np.zeros(runs)
mu_min = np.zeros(runs)

i = 0
while i < runs:
# for gene_ind in sig_gene_ind:
	gene_ind = int(gene_ind)
	print("gene_ind: {}".format(gene_ind))
	#converting it into r object for passing into r function
	with localconverter(robjects.default_converter + pandas2ri.converter):
		stage1_results = stage1_r(gene_ind)
		gene_name = robjects.conversion.rpy2py(stage1_results[0])
		#
		if gene_name[0] == "None":
			gene_names.append("None")
			i += 1
			continue
		#
		gene_names.append(gene_name[0])
		#
		mean_x_train = robjects.conversion.rpy2py(stage1_results[2])
		mean_x_val = robjects.conversion.rpy2py(stage1_results[3])
		mean_x_test = robjects.conversion.rpy2py(stage1_results[4])
		y_train = robjects.conversion.rpy2py(stage1_results[5])
		y_val = robjects.conversion.rpy2py(stage1_results[6])
		y_test = robjects.conversion.rpy2py(stage1_results[7])
		rsq[i] = robjects.conversion.rpy2py(stage1_results[8])
		adj_rsq[i] = robjects.conversion.rpy2py(stage1_results[9])
		LR_pval[i] = robjects.conversion.rpy2py(stage1_results[10])
		LR_tval[i] = robjects.conversion.rpy2py(stage1_results[11])
		stage1_fstat[i] = robjects.conversion.rpy2py(stage1_results[12])
		stage1_pval[i] = robjects.conversion.rpy2py(stage1_results[13])
		num_snps[i] = robjects.conversion.rpy2py(stage1_results[14])
		stage1_sigma[i] = robjects.conversion.rpy2py(stage1_results[15])
		# if the test set is too5.3 qt small, skip
		if mean_x_test.shape[0] < 150:
			IVa_pval_global[i] = 1
			IVa_pval_nonlinear[i] = 1
			continue

	#l2 = tune_l2(mean_x_train, y_train, (1,), 1, log_l2 = np.linspace(-5,4,5))
	l2 = 0
	rsp1 = create_stage2((1,),1, l2 = l2)
	_ = fit_stage2(rsp1, mean_x_train, y_train, mean_x_val, y_val,
				learning_rate = 0.00001, 
				batch_size = 64, 
				training_steps = 8000,
				patience = 300)
	# h = fit_stage2(rsp1, mean_x_train, y_train, mean_x_val, y_val)
	# json.dump(h.history, open('../results/stage2_results/history_{}_{}'.format(gene_ind,i), 'w'))
	# h = json.load(open('../results/stage2_results/history7468_{}'.format(i), 'r'))
	# np.savetxt('../results/stage2_results/l2_{}_{}.txt'.format(gene_ind,i), np.array([l2]))

	#tests
	tmp = rsp1.predict(mean_x_test).squeeze()
	IVa_pval_global[i], IVa_tval_global[i], IVa_pval_nonlinear[i], IVa_tval_nonlinear[i], nonlinear_coef[i], nonlinear_se[i] = stage2Tests(mean_x_test, y_test, tmp)
	#
	mu_max[i] = np.quantile(mean_x_test, 0.975)
	mu_min[i] = np.quantile(mean_x_test, 0.025)
	# save model
	if True:
		rsp1.save("../models/hdl/nfatc3_test/{}_{}".format(gene_ind, i))
	#
	i += 1

print("Length of gene_names: {}".format(len(gene_names)))
results = pd.DataFrame({"gene_names":gene_names, 
						"IVa_pval_global":IVa_pval_global, 
						"IVa_pval_nonlinear":IVa_pval_nonlinear,
						"IVa_tval_global":IVa_tval_global, 
						"IVa_tval_nonlinear":IVa_tval_nonlinear,
						"LR_pval":LR_pval,
						"LR_tval":LR_tval,
						"rsq":rsq,
						"adj_rsq":adj_rsq,
						"stage1_fstat":stage1_fstat,
						"stage1_pval":stage1_pval,
						"stage1_sigma":stage1_sigma,
						"nonlinear_coef":nonlinear_coef,
						"nonlinear_se":nonlinear_se,
						"mu_max":mu_max,
						"mu_min":mu_min,
						"num_snps":num_snps
						}
						)
results.to_csv("~/deepRIV/UKB/results/stage2_results/nfat3_test_{}.txt".format(gene_ind), index = False, sep = " ")


