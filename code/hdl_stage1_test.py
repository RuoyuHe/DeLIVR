import sys
print(sys.version)
import pickle
import pdb
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
from wrappers import create_stage1, create_stage2, fit_stage2, tune_l2
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt



# Defining the R script and loading the instance in Python
r = robjects.r
r['source']('hdl_stage1_FUN.R')

# Loading the function we have defined in R.
stage1_r = robjects.globalenv['stage1']

file_num = 1
runs = 150
start = (file_num-1)*runs+1
end = min(file_num*runs+1, 19696)
runs = end - start


# start = 22
# end = 23


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
nonlinear_coef = np.zeros(runs)
nonlinear_se = np.zeros(runs)
mu_max = np.zeros(runs)
mu_min = np.zeros(runs)

i=0
for gene_ind in range(start, end):
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
		# if the test set is too small, skip
		if mean_x_test.shape[0] < 150:
			IVa_pval_global[i] = 1
			IVa_pval_nonlinear[i] = 1
			continue

	# l2 = tune_l2(mean_x_train, y_train, (1,), 1, log_l2 = np.linspace(-5,4,5))
	# rsp1 = create_stage2((1,),1, l2 = l2)
	# _ = fit_stage2(rsp1, mean_x_train, y_train, mean_x_val, y_val)
	# #tests
	# tmp = rsp1.predict(mean_x_test).squeeze()
	# design = sm.add_constant(tmp.reshape(-1,1))
	# testM = sm.OLS(y_test, design)
	# testM_results = testM.fit()
	# IVa_pval_global[i] = testM_results.pvalues[-1]
	# IVa_tval_global[i] = testM_results.tvalues[-1]
	# # y~ E(X|Z) + E(g(X)|Z)
	# testM = sm.OLS(y_test, sm.add_constant(mean_x_test))
	# testM_results = testM.fit()
	# residuals = y_test - testM_results.predict(sm.add_constant(mean_x_test))
	# testM = sm.OLS(residuals, sm.add_constant(tmp.reshape(-1,1)))
	# testM_results = testM.fit()	
	# IVa_pval_nonlinear[i] = testM_results.pvalues[-1]
	# IVa_tval_nonlinear[i] = testM_results.tvalues[-1]
	# nonlinear_coef[i] = testM_results.params[-1]
	# nonlinear_se[i] = testM_results.bse[-1]
	# #
	# mu_max[i] = np.quantile(mean_x_test, 0.975)
	# mu_min[i] = np.quantile(mean_x_test, 0.025)
	# # save model
	# if IVa_pval_nonlinear[i] <= 0.05/100:
	# 	rsp1.save("../models/hdl/model_{}".format(gene_ind))
	# #
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
						"nonlinear_coef":nonlinear_coef,
						"nonlinear_se":nonlinear_se,
						"mu_max":mu_max,
						"mu_min":mu_min,
						"num_snps":num_snps
						}
						)
results.to_csv("~/deepRIV/UKB/results/stage1_results/p_vals_{}.txt".format(file_num), index = False, sep = " ")
