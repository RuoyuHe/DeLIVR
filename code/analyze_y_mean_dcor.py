import os
import sys
import dcor
import itertools
import numpy as np
import pandas as pd
import statsmodels.api as sm
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from sklearn.model_selection import train_test_split
from wrappers import stage2Tests
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
r['source']('/home/panwei/he000176/deepRIV/UKB/code/hdl/combine_p/hdl_stage1_FUN.R')

# Loading the function we have defined in R.
stage1_r = robjects.globalenv['stage1']

path = "/home/panwei/he000176/deepRIV/UKB/results/hdl/combine_p/repeat_ind/"
cauchy_results = pd.read_csv(path+"cauchy_results.csv")
repeat_genes = cauchy_results.loc[cauchy_results.IVa_pval_nonlinear <= 0.05/4701, 'gene_ind'].values

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

num_genes = 9
cutoff = 0.05/4701

cauchy_results = pd.DataFrame(np.zeros((repeat_genes.shape[0],3)))
cauchy_results.columns = ['gene_names', 'y_mean_dcor_pval', 'y_mean_nonlinear_pval']

l=0
for gene_ind in repeat_genes:
	gene_ind = int(gene_ind)
	print(gene_ind)
	#converting it into r object for passing into r function
	with localconverter(robjects.default_converter + pandas2ri.converter):
		stage1_results = stage1_r(gene_ind)
		gene_name = robjects.conversion.rpy2py(stage1_results[0])
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
	#
	cauchy_results.iloc[l,0] = gene_names
	IVa = np.zeros(21)
	IVa_pval_nonlinear = np.zeros(21)
	for j in range(num_params):
		# train test split
		tmp_data = pd.DataFrame(np.append(Y.reshape(-1,1), mean_X, axis = 1))
		for k in range(j*num_repeats, (j+1)*num_repeats):
			data = train_val_test_split(tmp_data, **split_ratio[j], random_state = k)
			unq1 = np.unique(data['mean_x_test'])
			unq = []
			mean_y = []
			for i in range(unq1.shape[0]):
				idx = np.where(data['mean_x_test'][:,0] == unq1[i])[0]
				if idx.shape[0] > 1000:
					unq.append(unq1[i])
					mean_y.append(np.mean(data['y_test'][idx]))
			unq = np.array(unq)
			mean_y = np.array(mean_y)
			# unique train
			unq1 = np.unique(data['mean_x_train'])
			unq_train = []
			mean_y_train = []
			for i in range(unq1.shape[0]):
				idx = np.where(data['mean_x_train'][:,0] == unq1[i])[0]
				if idx.shape[0] > 1000:
					unq_train.append(unq1[i])
					mean_y_train.append(np.mean(data['y_train'][idx]))
			unq_train = np.array(unq_train)
			mean_y_train = np.array(mean_y_train)
			#
			mean_x_test = np.array([])
			y_test = np.array([])
			mean_y_test = np.array([])
			# tmp_dict = dict(zip(unq, mean_y))
			tmp_dict = dict(zip(unq_train, mean_y_train))
			for i in range(data['mean_x_test'].shape[0]):
				if data['mean_x_test'][i,0] in tmp_dict:
					mean_x_test = np.append(mean_x_test, data['mean_x_test'][i,0])
					y_test = np.append(y_test, data['y_test'][i])
					mean_y_test = np.append(mean_y_test, tmp_dict[data['mean_x_test'][i,0]])
			mean_x_test = mean_x_test.reshape(-1,1)
			residuals_test = y_test - mean_y_test
			IVa[k] = dcor.independence.distance_correlation_t_test(mean_y_test.astype(np.float64), 
				residuals_test.astype(np.float64)).p_value
			_, _, IVa_pval_nonlinear[k], _, _, _ = stage2Tests(mean_x_test, y_test, mean_y_test)
		
	cauchy_dcor_T = np.tan((0.5-IVa)*np.pi).sum() / IVa.shape[0]
	cauchy_nl_T = np.tan((0.5-IVa_pval_nonlinear)*np.pi).sum() / IVa_pval_nonlinear.shape[0]
	cauchy_results.iloc[l,1] = 1/2 - np.arctan(cauchy_dcor_T)/np.pi
	cauchy_results.iloc[l,2] = 1/2 - np.arctan(cauchy_nl_T)/np.pi
	l += 1

cauchy_results.to_csv(path + "y_mean_cauchy.csv", index = False)
