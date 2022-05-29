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
from wrappers import create_stage1, create_stage2, fit_stage2, tune_l2, stage2Tests
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

def train_val_test_split(data, train_ratio, val_ratio, random_state = None):
	'''
		data: a pd.DataFrame of size n x p, first column is Y, the rest are X
	'''
	n = data.shape[0]
	train, val, test = np.split(data.sample(frac = 1, random_state = random_state), 
											[int(train_ratio * n), int((val_ratio + train_ratio) * n)])
	if train.shape[0]<=150 or val.shape[0]<=150 or test.shape[0] <= 150:
		return None
	y_train = train.iloc[:,0].to_numpy()
	y_val = val.iloc[:,0].to_numpy()
	y_test = test.iloc[:,0].to_numpy()
	x_train = train.iloc[:,1:].to_numpy()
	x_val = val.iloc[:,1:].to_numpy()
	x_test = test.iloc[:,1:].to_numpy()
	return {'y_train':y_train,
			'y_val':y_val,
			'y_test':y_test,
			'mean_x_train':x_train,
			'mean_x_val':x_val,
			'mean_x_test':x_test}


file_num = 0

# Defining the R script and loading the instance in Python
r = robjects.r
r['source']('/home/panwei/he000176/deepRIV/UKB/code/hdl/combine_p/hdl_stage1_FUN_unnormalized.R')

# Loading the function we have defined in R.
stage1_r = robjects.globalenv['stage1']

repeat_genes = pd.read_csv("repeat_genes_ind.csv", index_col = 0)
IVa_results = pd.read_csv("~/deepRIV/UKB/results/hdl/combine_p/ind_test/IVa_combined.csv")
idx = repeat_genes.reset_index().set_index('gene_name').loc[IVa_results.gene_name, 'index'].values
idx.sort()
gene_names = repeat_genes.loc[idx, 'gene_name'].values
repeat_genes = repeat_genes.loc[idx, 'gene_ind'].values
idx = None

# split_seeds = (ord('u'), ord('k'), ord('b'))
split_ratio = ({'train_ratio':0.5, 
				'val_ratio':0.1},
				{'train_ratio':0.4, 
				'val_ratio':0.1},
				{'train_ratio':0.3, 
				'val_ratio':0.1})

training_params = ({'learning_rate':0.00001,
					'decay_steps':150,
					'decay_rate':1,
					'patience':80},
					{'learning_rate':0.00001,
					'training_steps':500,
					'decay_rate':1,
					'patience':30},
					{'learning_rate':0.00001,
					'training_steps':500,
					'batch_size':32,
					'decay_rate':1,
					'patience':30})

num_params = len(training_params)
num_repeats = 7

gene_ind = repeat_genes[file_num]
gene_ind = int(gene_ind)
print("gene_ind: {}".format(gene_ind))

out_path = '/home/panwei/he000176/deepRIV/UKB/'
IVa_path = out_path + 'results/hdl/combine_p/IVa_{}.txt'.format(file_num)
debug_path = out_path + 'results/debug/combine_p/{}.txt'.format(file_num)

#converting it into r object for passing into r function
with localconverter(robjects.default_converter + pandas2ri.converter):
	stage1_results = stage1_r(gene_ind)
	# gene_name = robjects.conversion.rpy2py(stage1_results[0])
	#
	# gene_names = gene_name[0]
	#
	mean_X = robjects.conversion.rpy2py(stage1_results[2])
	Y = robjects.conversion.rpy2py(stage1_results[3])
	rsq = robjects.conversion.rpy2py(stage1_results[4])
	adj_rsq = robjects.conversion.rpy2py(stage1_results[5])
	LR_pval = robjects.conversion.rpy2py(stage1_results[6])
	LR_tval = robjects.conversion.rpy2py(stage1_results[7])
	stage1_fstat = robjects.conversion.rpy2py(stage1_results[8])
	stage1_pval = robjects.conversion.rpy2py(stage1_results[9])
	num_snps = robjects.conversion.rpy2py(stage1_results[10])
	stage1_sigma = robjects.conversion.rpy2py(stage1_results[11])
	stage2_theta = robjects.conversion.rpy2py(stage1_results[12])
	stage2_intercept = robjects.conversion.rpy2py(stage1_results[13])

#
with open(debug_path, "a") as f:
	f.write("{}\n".format(gene_ind))
	f.write("Common Done\n")


for j in range(num_params):
	# train test split
	tmp_data = pd.DataFrame(np.append(Y.reshape(-1,1), mean_X, axis = 1))

	for k in range(j*num_repeats, (j+1)*num_repeats):
		data = train_val_test_split(tmp_data, **split_ratio[j], random_state = k)
		if data is None or np.unique(data['mean_x_test']).shape[0] < 2:
			IVa_pval_global = 100
			IVa_pval_nonlinear = 100
			IVa_tval_global = 100
			IVa_tval_nonlinear = 100
			continue
		mu_max = np.quantile(data['mean_x_test'], 0.975)
		mu_min = np.quantile(data['mean_x_test'], 0.025)
		# normalize data
		scaler_x = StandardScaler().fit(data['mean_x_train'])
		scaler_y = StandardScaler().fit(data['y_train'].reshape(-1,1))
		data['mean_x_train'] = scaler_x.transform(data['mean_x_train'])
		data['mean_x_val'] = scaler_x.transform(data['mean_x_val'])
		data['mean_x_test'] = scaler_x.transform(data['mean_x_test'])
		data['y_train'] = scaler_y.transform(data['y_train'].reshape(-1,1)).squeeze()
		data['y_val'] = scaler_y.transform(data['y_val'].reshape(-1,1)).squeeze()
		data['y_test'] = scaler_y.transform(data['y_test'].reshape(-1,1)).squeeze()
		# train
		l2 = 0
		rsp1 = create_stage2((1,),1, l2 = l2)
		_ = fit_stage2(rsp1, data['mean_x_train'], data['y_train'], data['mean_x_val'], data['y_val'],
						**training_params[j])
		# tests
		tmp = rsp1.predict(data['mean_x_test']).squeeze()
		IVa_pval_global, IVa_tval_global, IVa_pval_nonlinear, IVa_tval_nonlinear, nonlinear_coef, nonlinear_se = stage2Tests(data['mean_x_test'], data['y_test'], tmp)
		# save model
		if True:
			model_path = out_path + 'models/hdl/combine_p/model_{}_{}_{}'.format(gene_ind, j, k)
			rsp1.save(model_path)
		with open(debug_path, "a") as f:
			f.write("{}: Done\n".format(k))
		#
		IVa_results = pd.DataFrame({"gene_name":gene_names[file_num],
									"IVa_pval_global":IVa_pval_global, 
									"IVa_pval_nonlinear":IVa_pval_nonlinear,
									"IVa_tval_global":IVa_tval_global, 
									"IVa_tval_nonlinear":IVa_tval_nonlinear, 
									"nonlinear_coef":nonlinear_coef, 
									"nonlinear_se":nonlinear_se,
									"mu_max":mu_max, 
									"mu_min":mu_min
									}, index = [0])

		IVa_results.to_csv(IVa_path, mode = 'a', index = False, 
								sep = " ", header = not os.path.exists(IVa_path))