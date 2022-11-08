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
from utils import create_stage1, create_stage2, fit_stage2, tune_l2, stage2Tests
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
	x_train = train.iloc[:,1].to_numpy()
	x_val = val.iloc[:,1].to_numpy()
	x_test = test.iloc[:,1].to_numpy()
	z_train = train.iloc[:,2:].to_numpy()
	z_val = val.iloc[:,2:].to_numpy()
	z_test = test.iloc[:,2:].to_numpy()
	return {'y_train':y_train,
			'y_val':y_val,
			'y_test':y_test,
			'x_train':x_train,
			'x_val':x_val,
			'x_test':x_test,
			'z_train':z_train,
			'z_val':z_val,
			'z_test':z_test}

file_num = int(sys.argv[-1])

# split_seeds = (ord('u'), ord('k'), ord('b'))
# split_ratio = [{'train_ratio':0.5, 
# 				'val_ratio':0.1},
# 				{'train_ratio':0.35, 
# 				'val_ratio':0.1}]

# training_params = [{'learning_rate':0.00001,
# 					'training_steps':500,
# 					'decay_steps':150,
# 					'decay_rate':1,
# 					'patience':60},
# 					{'learning_rate':0.00001,
# 					'training_steps':500,
# 					'decay_steps':150,
# 					'decay_rate':1,
# 					'patience':30}]

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
runs = 10

gene_ind = 12441
print("gene_ind: {}".format(gene_ind))

true_model = 'typeI'

out_path = '/home/panwei/he000176/deepRIV/UKB/simulation/'
IVa_path = out_path + 'results/cdk2ap1/20repeats/delivr_' + true_model + '_{}.txt'.format(file_num)
debug_path = out_path + 'results/debug/cdk2ap1/' + true_model + '_{}.txt'.format(file_num)

z_path = out_path + 'data/cdk2ap1/ukb_snp_bed.csv'
beta_path = out_path + 'data/cdk2ap1/hatbetaX1.txt'
theta_path = out_path + 'data/cdk2ap1/hattheta_lq.txt'


z = pd.read_csv(z_path).to_numpy()
beta = pd.read_csv(beta_path, sep = ' ').to_numpy()
nonzero_idx = np.where(beta.squeeze()!=0)[0]
z = z[:, nonzero_idx]
beta = beta[nonzero_idx, :]
theta = pd.read_csv(theta_path, sep = ' ').to_numpy()

true_g = lambda x: 0
gamma = 0.7

for run in range(runs):
	e = np.random.multivariate_normal(mean = [0,0], cov = [[1,gamma],[gamma,1]], size = z.shape[0])
	x = (z@beta).squeeze() + e[:,0]
	y = true_g(x) + e[:,1]
	x = x.reshape(-1,1)
	tmp_data = pd.DataFrame(np.concatenate((y.reshape(-1,1), x, z), axis = 1))
	for j in range(num_params):
		# np.random.seed(k)
		for k in range(j*num_repeats, (j+1)*num_repeats):
			# train test split
			data = train_val_test_split(tmp_data, **split_ratio[j], random_state = k)

			design = sm.add_constant(data['z_train'])
			stage1 = sm.OLS(data['x_train'], design)
			stage1 = stage1.fit()
			data['mean_x_train'] = stage1.predict(design).reshape(-1,1)
			data['mean_x_val'] = stage1.predict(sm.add_constant(data['z_val'])).reshape(-1,1)
			data['mean_x_test'] = stage1.predict(sm.add_constant(data['z_test'])).reshape(-1,1)

			rsq = stage1.rsquared
			adj_rsq = stage1.rsquared_adj
			stage1_fstat = stage1.fvalue
			stage1_pval = stage1.f_pvalue
			stage1_sigma = stage1.ssr/ stage1.df_resid

			# TWAS-L, TWAS-LQ
			design = sm.add_constant(data['z_test'])
			stage1_l = sm.OLS(data['x_test'], design)
			stage1_l = stage1_l.fit()
			stage1_q = sm.OLS(data['x_test']**2, design)
			stage1_q = stage1_q.fit()
			mean_X = stage1_l.predict(design).reshape(-1,1)
			mean_X2 = stage1_q.predict(design).reshape(-1,1)

			design = sm.add_constant(mean_X)
			stage2_l = sm.OLS(data['y_test'], design)
			stage2_l = stage2_l.fit()
			L_pval = stage2_l.pvalues[-1]
			L_coef = stage2_l.params[-1]
			design = sm.add_constant(np.append(mean_X, mean_X2, axis = 1))
			stage2_lq = sm.OLS(data['y_test'], design)
			stage2_lq = stage2_lq.fit()
			LQ_pval = stage2_lq.f_pvalue

			if np.unique(data['mean_x_test']).shape[0] < 2:
				with open(debug_path, "a") as f:
					f.write("Data problem, skip\n")
				continue
			mu_max = np.quantile(data['mean_x_test'], 0.975)
			mu_min = np.quantile(data['mean_x_test'], 0.025)
			# train
			l2 = 0
			rsp1 = create_stage2((1,),1, l2 = l2)
			h = fit_stage2(rsp1, data['mean_x_train'], data['y_train'], data['mean_x_val'], data['y_val'],
										**training_params[j])
			# json.dump(h.history, open(out_path + 'results/cdk2ap1/history_{}'.format( k), 'w'))
			# h = json.load(open(out_path + 'results/cdk2ap1/history_{}'.format(k), 'r'))
			#tests
			tmp = rsp1.predict(data['mean_x_test']).squeeze()
			IVa_pval_global, IVa_tval_global, IVa_pval_nonlinear, IVa_tval_nonlinear, nonlinear_coef, nonlinear_se, linear_coef, linear_se = stage2Tests(data['mean_x_test'], data['y_test'], tmp)
			# save model
			if False:
				model_path = out_path + 'models/cdk2ap1/delivr_' + true_model + '_{}_{}'.format(k,j)
				rsp1.save(model_path)
			with open(debug_path, "a") as f:
				f.write("{}: Done\n".format(k))

			IVa_results = pd.DataFrame({"gene_names":'CDK2AP1',
										"twasl_pval":L_pval,
										"twaslq_pval":LQ_pval,
										"IVa_pval_global":IVa_pval_global, 
										"IVa_pval_nonlinear":IVa_pval_nonlinear,
										"IVa_tval_global":IVa_tval_global, 
										"IVa_tval_nonlinear":IVa_tval_nonlinear, 
										"nonlinear_coef":nonlinear_coef, 
										"nonlinear_se":nonlinear_se,
										"linear_coef":linear_coef, 
										"linear_se":linear_se,
										"mu_max":mu_max, 
										"mu_min":mu_min,
										"gene_ind":gene_ind,
										"num_snps":z.shape[1]
										}, index = [0])

			IVa_results.to_csv(IVa_path, mode = 'a', index = False, 
									sep = " ", header = not os.path.exists(IVa_path))

