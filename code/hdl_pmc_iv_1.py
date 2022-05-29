# hdl_pmc_iv
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
from tensorflow.keras.layers import Activation, Add, Input, Dense, LeakyReLU, ReLU, Reshape
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.optimizers import Adam
from sklearn.model_selection import StratifiedShuffleSplit
from tensorflow.keras.backend import expand_dims, squeeze, clear_session
from sklearn.model_selection import train_test_split
from neuralnet import stage2Loss, MCsample, ResNet, directM, MC2sls
from fit_model import deeprivOPT, fit_deepRIV, unbiased_grad, fit_deepIV
from wrappers import create_stage1, create_stage2, fit_stage2, tune_l2, stage2Tests
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt


# Defining the R script and loading the instance in Python
r = robjects.r
r['source']('hdl_stage1_FUN2.R')

# Loading the function we have defined in R.
stage1_r = robjects.globalenv['stage1']

file_num = 3
sig_gene_ind = [4639, 12030, 12441, 15043]
# gene_ind = 15043
runs = len(sig_gene_ind)

# runs = 10


gene_names = []
num_snps = np.zeros(runs)
PMC_pval_global = np.zeros(runs)
PMC_pval_nonlinear = np.zeros(runs)
PMC_tval_global = np.zeros(runs)
PMC_tval_nonlinear = np.zeros(runs)
PMC_nonlinear_coef = np.zeros(runs)
PMC_nonlinear_se = np.zeros(runs)
IV_pval_global = np.zeros(runs)
IV_pval_nonlinear = np.zeros(runs)
IV_tval_global = np.zeros(runs)
IV_tval_nonlinear = np.zeros(runs)
IV_nonlinear_coef = np.zeros(runs)
IV_nonlinear_se = np.zeros(runs)
LR_pval = np.zeros(runs)
LR_tval = np.zeros(runs)
rsq = np.zeros(runs)
adj_rsq = np.zeros(runs)
stage1_fstat = np.zeros(runs)
stage1_pval = np.zeros(runs)
stage1_sigma = np.zeros(runs)
mu_max = np.zeros(runs)
mu_min = np.zeros(runs)

i = 0
# while i < runs:
for gene_ind in sig_gene_ind:
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
	#
	mu_max[i] = np.quantile(mean_x_test, 0.975)
	mu_min[i] = np.quantile(mean_x_test, 0.025)
	# DeepIV
	n_draws = 60
	learning_rate = 0.00033
	batch_size = 64
	training_steps = 15000
	display_step = 1
	resample_prop = 1

	pred_x_train = np.random.normal(mean_x_train, stage1_sigma[i], (mean_x_train.shape[0],n_draws))
	pred_x_train2 = np.random.normal(mean_x_train, stage1_sigma[i], (mean_x_train.shape[0],n_draws))
	pred_x_val = np.random.normal(mean_x_val, stage1_sigma[i], (mean_x_val.shape[0],n_draws))

	deepiv = ResNet([8,8,8], LeakyReLU(alpha = 0.05), n_draws, l2 = 2)
	lr_schedule = tf.keras.optimizers.schedules.ExponentialDecay(
		learning_rate,
		decay_steps=8000,
		decay_rate=0.7,
		staircase=True)
	optimizer = tf.keras.optimizers.Adam(lr_schedule)
	fit_deepIV(deepiv, optimizer, unbiased_grad, stage2Loss, mean_x_train, stage1_sigma[i], 
					pred_x_train, pred_x_train2, pred_x_val, y_train, y_val,
					n_draws, batch_size, training_steps,
					display_step = 1, resample_prop = 1, early_stopping = 500)

	# save model
	deepiv.save("/home/panwei/he000176/deepRIV/UKB/models/hdl/deepiv_{}_{}".format(gene_ind, gene_name[0]))
	#tests
	pred_x_test = np.random.normal(mean_x_test, stage1_sigma[i], (mean_x_test.shape[0], n_draws))
	# deepIV
	tmp = deepiv(pred_x_test.reshape(-1,1), training = False).numpy()
	tmp = np.mean(tmp.reshape(-1,n_draws), axis=1)
	IV_pval_global[i], IV_tval_global[i], IV_pval_nonlinear[i], IV_tval_nonlinear[i], IV_nonlinear_coef[i], IV_nonlinear_se[i] = stage2Tests(mean_x_test, y_test, tmp)
	#
	# PMC
	# train parametric (2SLS)
	learning_rate = 0.001
	training_steps = 120000
	m_p = MC2sls(1)
	lr_schedule = tf.keras.optimizers.schedules.ExponentialDecay(
		learning_rate,
		decay_steps=50000,
		decay_rate=0.7,
		staircase=True)
	optimizer = tf.keras.optimizers.Adam(lr_schedule)
	fit_deepIV(m_p, optimizer, unbiased_grad, stage2Loss, mean_x_train, stage1_sigma[i], 
					pred_x_train, pred_x_train2, pred_x_val, y_train, y_val,
					n_draws, batch_size, training_steps,
					display_step = 1, resample_prop = 1, early_stopping = 500)
	# save model
	m_p_weights = pd.DataFrame({"w0":m_p.w0.numpy().squeeze(), "w":m_p.w.numpy().squeeze(), 
								"b":m_p.b.numpy().squeeze()}, index=[0])
	m_p_weights.to_csv("~/deepRIV/UKB/models/hdl/pmc_{}_{}.txt".format(gene_ind, gene_name[0]),
						index = False, sep = " ")
	# tests
	tmp = m_p(pred_x_test.reshape(-1,1), training = False).numpy()
	tmp = np.mean(tmp.reshape(-1,n_draws), axis=1)
	PMC_pval_global[i], PMC_tval_global[i], PMC_pval_nonlinear[i], PMC_tval_nonlinear[i], PMC_nonlinear_coef[i], PMC_nonlinear_se[i] = stage2Tests(mean_x_test, y_test, tmp)
	#
	i += 1

print("Length of gene_names: {}".format(len(gene_names)))
results = pd.DataFrame({"gene_names":gene_names, 
						"IV_pval_global":IV_pval_global, 
						"IV_pval_nonlinear":IV_pval_nonlinear,
						"IV_tval_global":IV_tval_global, 
						"IV_tval_nonlinear":IV_tval_nonlinear,
						"IV_nonlinear_coef":IV_nonlinear_coef,
						"IV_nonlinear_se":IV_nonlinear_se,
						"PMC_pval_global":PMC_pval_global, 
						"PMC_pval_nonlinear":PMC_pval_nonlinear,
						"PMC_tval_global":PMC_tval_global, 
						"PMC_tval_nonlinear":PMC_tval_nonlinear,
						"PMC_nonlinear_coef":PMC_nonlinear_coef,
						"PMC_nonlinear_se":PMC_nonlinear_se,
						"LR_pval":LR_pval,
						"LR_tval":LR_tval,
						"rsq":rsq,
						"adj_rsq":adj_rsq,
						"stage1_fstat":stage1_fstat,
						"stage1_pval":stage1_pval,
						"stage1_sigma":stage1_sigma,
						"mu_max":mu_max,
						"mu_min":mu_min,
						"num_snps":num_snps
						}
						)
results.to_csv("~/deepRIV/UKB/results/stage2_results/iv_pmc_{}.txt".format(file_num), index = False, sep = " ")


