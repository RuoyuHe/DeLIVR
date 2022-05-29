# included independent tests for Y_hat and residuals
import os
import sys
print(sys.version)
import pdb
import dcor
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
from wrappers import create_stage1, create_stage2, fit_stage2, tune_l2, stage2Tests
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

def train_val_test_split(data, train_ratio, val_ratio, random_state):
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


# Defining the R script and loading the instance in Python
r = robjects.r
r['source']('/home/panwei/he000176/deepRIV/UKB/code/hdl/combine_p/hdl_stage1_FUN.R')

# Loading the function we have defined in R.
stage1_r = robjects.globalenv['stage1']

file_num = 1
runs = 30
start = (file_num-1)*runs
end = min(file_num*runs, 4701)
runs = end - start

hdl_path = "/home/panwei/he000176/deepRIV/UKB/results/hdl/test40p_unrelated2/"
hdl = pd.read_csv(hdl_path + "hdl_combined_results.csv")
valid_gene_inds = hdl['gene_ind'].iloc[start:end].to_numpy()


# split_seeds = (ord('u'), ord('k'), ord('b'))
split_ratio = ({'train_ratio':0.5, 
				'val_ratio':0.1, 
				'random_state':ord('u')},
				{'train_ratio':0.4, 
				'val_ratio':0.1, 
				'random_state':ord('k')},
				{'train_ratio':0.3, 
				'val_ratio':0.1, 
				'random_state':ord('b')})

training_params = ({'learning_rate':0.00001,
					'training_steps':300,
					'decay_steps':150,
					'decay_rate':1,
					'patience':60},
					{'learning_rate':0.00001,
					'training_steps':300,
					'decay_rate':1,
					'patience':30},
					{'learning_rate':0.00001,
					'batch_size':32,
					'decay_rate':1,
					'patience':30})

num_repeats = len(training_params)

out_path = '/home/panwei/he000176/deepRIV/UKB/'
common_path = out_path + 'results/hdl/combine_p/common_{}.txt'.format(file_num)
debug_path = out_path + "results/debug/combine_p/{}.txt".format(file_num)

for gene_ind in valid_gene_inds:
	gene_ind = int(gene_ind)
	print("gene_ind: {}".format(gene_ind))
	#converting it into r object for passing into r function
	with localconverter(robjects.default_converter + pandas2ri.converter):
		stage1_results = stage1_r(gene_ind)
		gene_name = robjects.conversion.rpy2py(stage1_results[0])
		#
		if gene_name[0] == "None":
			print(1)
			continue
		gene_names = gene_name[0]
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
		# if the test set is too small, skip
		if mean_X.shape[0] < 500:
			print(1)
			continue
	# independent tests (screening)
	Y_hat = mean_X @ stage2_theta + stage2_intercept
	stage2_residuals = Y - Y_hat
	common_dcor_test = dcor.independence.distance_correlation_t_test(Y_hat, stage2_residuals)
	#
	common_results = pd.DataFrame({"gene_names":gene_names,
						"LR_pval":LR_pval,
						"LR_tval":LR_tval,
						"dcor_pval":common_dcor_test.p_value,
						"rsq":rsq,
						"adj_rsq":adj_rsq,
						"stage1_fstat":stage1_fstat,
						"stage1_pval":stage1_pval,
						"stage1_sigma":stage1_sigma,
						"num_snps":num_snps,
						"gene_ind":gene_ind
						})
	common_results.to_csv(common_path, mode = 'a', index = False, 
							sep = " ", header = not os.path.exists(common_path))
	with open(debug_path, "a") as f:
		f.write("{}\n".format(gene_ind))
		f.write("Common Done\n")

	if common_dcor_test.p_value <= 0.05/4701:
		IVa_pval_global = np.zeros(num_repeats)
		IVa_pval_nonlinear = np.zeros(num_repeats)
		IVa_tval_global = np.zeros(num_repeats)
		IVa_tval_nonlinear = np.zeros(num_repeats)
		IVa_dcor_pval = np.zeros(num_repeats)
		nonlinear_coef = np.zeros(num_repeats)
		nonlinear_se = np.zeros(num_repeats)
		mu_max = np.zeros(num_repeats)
		mu_min = np.zeros(num_repeats)
		for j in range(num_repeats):
			# train test split
			tmp_data = pd.DataFrame(np.append(Y.reshape(-1,1), mean_X, axis = 1))
			data = train_val_test_split(tmp_data, **split_ratio[j])
			if data is None or np.unique(data['mean_x_test']).shape[0] < 2:
				IVa_pval_global[j] = 100
				IVa_pval_nonlinear[j] = 100
				IVa_tval_global[j] = 100
				IVa_tval_nonlinear[j] = 100
				continue
			mu_max[j] = np.quantile(data['mean_x_test'], 0.975)
			mu_min[j] = np.quantile(data['mean_x_test'], 0.025)
			# train
			l2 = 0
			rsp1 = create_stage2((1,),1, l2 = l2)
			_ = fit_stage2(rsp1, data['mean_x_train'], data['y_train'], data['mean_x_val'], data['y_val'],
							**training_params[j])
			# tests
			tmp = rsp1.predict(data['mean_x_test']).squeeze()
			IVa_pval_global[j], IVa_tval_global[j], IVa_pval_nonlinear[j], IVa_tval_nonlinear[j], nonlinear_coef[j], nonlinear_se[j] = stage2Tests(data['mean_x_test'], data['y_test'], tmp)
			# dcor tests
			residuals_test = data['y_test'] - tmp
			IVa_dcor_pval[j] = dcor.independence.distance_correlation_t_test(tmp.astype(np.float64), 
				residuals_test.astype(np.float64)).p_value
			# save model
			if IVa_pval_nonlinear[j] <= 0.05/100:
				rsp1.save("/home/panwei/he000176/deepRIV/UKB/models/hdl/combine_p/model_{}_{}".format(gene_ind, j))
			with open(debug_path, "a") as f:
				f.write("{}: Done\n".format(j))
		#
		gene_names = np.array(gene_names).reshape(-1,1)
		IVa_results = pd.DataFrame(gene_names)
		IVa_results = pd.concat([IVa_results, pd.DataFrame(np.concatenate((IVa_pval_global, IVa_pval_nonlinear,
													IVa_tval_global, IVa_tval_nonlinear, 
													IVa_dcor_pval,
													nonlinear_coef, nonlinear_se,
													mu_max, mu_min)).reshape(1,-1))], axis = 1)
		colnames = [['gene_names'], ['global_p']*num_repeats, ['nonlinear_p']*num_repeats,
					['global_t']*num_repeats, ['nonlinear_t']*num_repeats, 
					['dcor_p']*num_repeats,
					['nonlinear_coef']*num_repeats, ['nonlinear_se']*num_repeats,
					['mu_max']*num_repeats, ['mu_min']*num_repeats]
		IVa_results.columns = list(itertools.chain(*colnames))
		IVa_path = out_path + 'results/hdl/combine_p/IVa_{}.txt'.format(file_num)
		IVa_results.to_csv(IVa_path, mode = 'a', index = False, 
								sep = " ", header = not os.path.exists(IVa_path))

