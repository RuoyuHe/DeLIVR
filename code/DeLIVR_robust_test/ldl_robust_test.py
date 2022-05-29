import os
import sys
print(sys.version)
import pdb
import json
import pickle
import itertools
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
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from wrappers import create_stage2, fit_stage2, stage2Tests
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt


file_num = int(sys.argv[-1])

# Defining the R script and loading the instance in Python
r = robjects.r
r['source']('/home/panwei/he000176/deepRIV/UKB/code/robust_DeLIVR/ldl_stage1_FUN_unnormalized.R')

# Loading the function we have defined in R.
stage1_r = robjects.globalenv['stage1']

repeat_genes = pd.read_csv("/home/panwei/he000176/deepRIV/UKB/results/ldl/combine_p/cauchy_NL_results.csv", index_col = 0)
hdl = pd.read_csv("/home/panwei/he000176/deepRIV/UKB/results/hdl/test40p_unrelated2/hdl_combined_results.csv")
repeat_genes = repeat_genes.merge(hdl, how = 'left', on = 'gene_name')
gene_names = repeat_genes.gene_name.values
repeat_genes = repeat_genes.gene_ind.values


# split_seeds = (ord('u'), ord('k'), ord('b'))
split_ratio = ({'test_ratio':0.4, 
				'val_ratio':0.1},
				{'test_ratio':0.5, 
				'val_ratio':0.1},
				{'test_ratio':0.5, 
				'val_ratio':0.1})

training_params = ({'learning_rate':0.00001,
					'decay_steps':150,
					'decay_rate':1,
					'patience':80},
					{'learning_rate':0.00001,
					'training_steps':500,
					'decay_rate':1,
					'patience':60},
					{'learning_rate':0.00001,
					'training_steps':500,
					'decay_rate':1,
					'patience':30})

num_params = len(training_params)
num_repeats = 7

gene_ind = repeat_genes[file_num]
gene_ind = int(gene_ind)
print("gene_ind: {}".format(gene_ind))

out_path = '/home/panwei/he000176/deepRIV/UKB/'
IVa_path = out_path + 'results/DeLIVR_robust_test/ldl_robust_{}.txt'.format(file_num)
debug_path = out_path + 'results/debug/DeLIVR_robust_test/ldl_{}.txt'.format(file_num)

#converting it into r object for passing into r function
with localconverter(robjects.default_converter + pandas2ri.converter):
	stage1_results = stage1_r(gene_ind)
	mean_X = robjects.conversion.rpy2py(stage1_results[2])
	Y = robjects.conversion.rpy2py(stage1_results[3])
	Z = robjects.conversion.rpy2py(stage1_results[14])
	beta = robjects.conversion.rpy2py(stage1_results[15])

nonzero_idx = np.where(beta.squeeze()!=0)[0]
Z = Z[:, nonzero_idx]
# add a column of ones
Z = np.concatenate((np.ones(Z.shape[0]).reshape(-1,1),Z),axis=-1)

for j in range(num_params):
	for k in range(j*num_repeats, (j+1)*num_repeats):
		# train test split
		test_ratio = split_ratio[j]['test_ratio']
		val_ratio = split_ratio[j]['val_ratio']/(1-test_ratio)
		z_train, z_test, mean_x_train, mean_x_test, y_train, y_test = train_test_split(Z, mean_X, Y, test_size=test_ratio, random_state=k)
		z_train, z_val, mean_x_train, mean_x_val, y_train, y_val = train_test_split(z_train, mean_x_train, y_train, test_size=val_ratio, random_state=k) 
		# normalize data
		scaler_x = StandardScaler().fit(mean_x_train)
		scaler_y = StandardScaler().fit(y_train.reshape(-1,1))
		mean_x_train = scaler_x.transform(mean_x_train)
		mean_x_val = scaler_x.transform(mean_x_val)
		mean_x_test = scaler_x.transform(mean_x_test)
		y_train = scaler_y.transform(y_train.reshape(-1,1)).squeeze()
		y_val = scaler_y.transform(y_val.reshape(-1,1)).squeeze()
		y_test = scaler_y.transform(y_test.reshape(-1,1)).squeeze()
		# train
		W_train = z_train @ np.linalg.inv(z_train.T @ z_train)
		W_val = z_val @ np.linalg.inv(z_val.T @ z_val)
		W_test = z_test @ np.linalg.inv(z_test.T @ z_test)
		QzY_train = y_train - W_train @ (z_train.T @ y_train)
		QzY_val = y_val - W_val @ (z_val.T @ y_val)
		QzY_test = y_test - W_test @ (z_test.T @ y_test)
		# train DeLIVR
		l2 = 0
		rsp1 = create_stage2((1,),1, l2 = l2)
		_ = fit_stage2(rsp1, mean_x_train, y_train, mean_x_val, y_val,
						**training_params[j])
		# save model
		model_path = out_path + 'models/ldl/combine_p/model_{}_{}_{}'.format(gene_ind, j, k)
		rsp1.save(model_path)
		# tests normal
		pred = rsp1.predict(mean_x_test).squeeze()
		IVa_pval_global, IVa_tval_global, IVa_pval_nonlinear, IVa_tval_nonlinear, _, _ = stage2Tests(mean_x_test, y_test, pred)
		# robust
		design = sm.add_constant(pred - W_test@(z_test.T@pred))
		testM = sm.OLS(QzY_test, design)
		testM_results = testM.fit()
		p_robust = testM_results.pvalues[-1]
		t_robust = testM_results.tvalues[-1]

		with open(debug_path, "a") as f:
			f.write("{}: Done\n".format(k))

		#
		IVa_results = pd.DataFrame({"gene_name":gene_names[file_num],
									"IVa_pval_global":IVa_pval_global, 
									"IVa_pval_nonlinear":IVa_pval_nonlinear,
									"IVa_tval_global":IVa_tval_global, 
									"IVa_tval_nonlinear":IVa_tval_nonlinear, 
									"p_robust":p_robust,
									"t_robust":t_robust
									}, index = [0])

		IVa_results.to_csv(IVa_path, mode = 'a', index = False, 
								sep = " ", header = not os.path.exists(IVa_path))