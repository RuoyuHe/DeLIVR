import os
import sys
print(sys.version)
import pdb
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
from sklearn.model_selection import train_test_split
from neuralnet import stage2Loss
from fit_model import deeprivOPT, fit_deepRIV, unbiased_grad, fit_deepIV, evaluate
from wrappers import create_stage1, create_stage2, stage2Tests
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

tf.config.run_functions_eagerly(True)

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


file_num = int(sys.argv[-1])

# Defining the R script and loading the instance in Python
r = robjects.r
r['source']('/home/panwei/he000176/deepRIV/UKB/code/hdl/combine_p/hdl_stage1_FUN.R')

# Loading the function we have defined in R.
stage1_r = robjects.globalenv['stage1']

repeat_genes = pd.read_csv("~/deepRIV/UKB/code/hdl/combine_p/deepIV/hdl_deepiv_genes.csv")
hdl = pd.read_csv("/home/panwei/he000176/deepRIV/UKB/results/hdl/test40p_unrelated2/hdl_combined_results.csv")
repeat_genes['gene_ind'] = hdl.reset_index().set_index('gene_name').loc[repeat_genes.x, 'gene_ind'].values

# split_seeds = (ord('u'), ord('k'), ord('b'))
split_ratio = ({'train_ratio':0.5, 
				'val_ratio':0.1},
				{'train_ratio':0.4, 
				'val_ratio':0.1})

training_params = ({'n_draws':16,
					'learning_rate':0.00033,
					'batch_size':640,
					'training_steps':92503,
					'display_step':1,
					'resample_prop':1,
					'early_stopping':5000},
					{'n_draws':8,
					'learning_rate':0.00033,
					'batch_size':640,
					'training_steps':92503,
					'display_step':1,
					'resample_prop':1,
					'early_stopping':3000})

num_params = len(training_params)
num_repeats = 10

gene_ind = repeat_genes.gene_ind[file_num]
gene_ind = int(gene_ind)
print("gene_ind: {}".format(gene_ind))
gene_name = hdl.loc[hdl.gene_ind==gene_ind,'gene_name'].values[0]

out_path = '/home/panwei/he000176/deepRIV/UKB/'
IVa_path = out_path + 'results/hdl/combine_p/repeat_deepiv/IV_{}.txt'.format(file_num)
debug_path = out_path + "results/debug/hdl_deepiv/{}.txt".format(file_num)

#converting it into r object for passing into r function
with localconverter(robjects.default_converter + pandas2ri.converter):
	stage1_results = stage1_r(gene_ind)
	# gene_name = robjects.conversion.rpy2py(stage1_results[0])
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
		loss_path = out_path + 'results/hdl/combine_p/repeat_deepiv/loss_{}_{}_{}.txt'.format(file_num,j,k)
		data = train_val_test_split(tmp_data, **split_ratio[j], random_state = k)
		if data is None or np.unique(data['mean_x_test']).shape[0] < 2:
			IVa_pval_global = 100
			IVa_pval_nonlinear = 100
			IVa_tval_global = 100
			IVa_tval_nonlinear = 100
			continue
		mu_max = np.quantile(data['mean_x_test'], 0.975)
		mu_min = np.quantile(data['mean_x_test'], 0.025)
		# train
		std_x_train = stage1_sigma
		pred_x_train = np.random.normal(data['mean_x_train'], std_x_train, 
										(data['mean_x_train'].shape[0],training_params[j]['n_draws']))
		pred_x_train2 = np.random.normal(data['mean_x_train'], std_x_train, 
										(data['mean_x_train'].shape[0],training_params[j]['n_draws']))
		pred_x_val = np.random.normal(data['mean_x_val'], std_x_train, 
									(data['mean_x_val'].shape[0],training_params[j]['n_draws']))
		pred_x_test = np.random.normal(data['mean_x_test'], std_x_train, 
									(data['mean_x_test'].shape[0],training_params[j]['n_draws']))

		stage2 = create_stage2((1,),1, l2 = 2)
		lr_schedule = tf.keras.optimizers.schedules.ExponentialDecay(
			training_params[j]['learning_rate'],
			decay_steps=8000,
			decay_rate=0.9,
			staircase=True)
		optimizer = tf.keras.optimizers.Adam(lr_schedule)
		fit_deepIV(stage2, optimizer, unbiased_grad, stage2Loss, data['mean_x_train'], std_x_train, 
							pred_x_train, pred_x_train2, pred_x_val, data['y_train'], data['y_val'],
							training_params[j]['n_draws'], training_params[j]['batch_size'], 
							training_params[j]['training_steps'], 
							display_step = 1, resample_prop = 1, 
							early_stopping = training_params[j]['early_stopping'],
							saving_path = loss_path)
		tmp = evaluate(stage2, inputs = pred_x_test.reshape(-1,1))
		tmp = tf.reduce_mean(tf.reshape(tmp, (-1, training_params[j]['n_draws'])) ,axis = 1).numpy()
		IV_pval_global, IV_tval_global, IV_pval_nonlinear, IV_tval_nonlinear, nonlinear_coef, nonlinear_se = stage2Tests(data['mean_x_test'], data['y_test'], tmp)

		# save model
		if True:
			stage2.save("/home/panwei/he000176/deepRIV/UKB/models/hdl/repeat_deepiv/model_{}_{}_{}".format(gene_ind, j, k))
		with open(debug_path, "a") as f:
			f.write("{}: Done\n".format(k))
		#
		IV_results = pd.DataFrame({"gene_name":gene_name,
									"IV_pval_global":IV_pval_global, 
									"IV_pval_nonlinear":IV_pval_nonlinear,
									"IV_tval_global":IV_tval_global, 
									"IV_tval_nonlinear":IV_tval_nonlinear, 
									"nonlinear_coef":nonlinear_coef, 
									"nonlinear_se":nonlinear_se,
									"mu_max":mu_max, 
									"mu_min":mu_min}, index = [0])

		IV_results.to_csv(IVa_path, mode = 'a', index = False, 
								sep = " ", header = not os.path.exists(IVa_path))
