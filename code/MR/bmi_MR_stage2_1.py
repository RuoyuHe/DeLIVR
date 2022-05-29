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


exposures = pd.read_csv('~/deepRIV/UKB/data/MR_tmp/lm_bmi_adjusted_cov.csv').iloc[:,0]
outcomes = pd.read_csv('~/deepRIV/UKB/data/MR_tmp/normalized_outcomes_bmi.csv')

file_num = int(sys.argv[-1])
exposure = 'bmi'
outcome = outcomes.columns[file_num]
mean_X = exposures.to_numpy().reshape(-1,1)
Y = outcomes.iloc[:,file_num].to_numpy().reshape(-1,1)

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

out_path = '/home/panwei/he000176/deepRIV/UKB/'
common_path = out_path + 'results/MR/' + exposure + '/' + outcome + '_common_{}.txt'.format(file_num)
IVa_path = out_path + 'results/MR/' + exposure + '/' + outcome + '_IVa_{}.txt'.format(file_num)
debug_path = out_path + "results/debug/MR/{}.txt".format(file_num)


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
		# train
		l2 = 0
		rsp1 = create_stage2((1,),1, l2 = l2)
		_ = fit_stage2(rsp1, data['mean_x_train'], data['y_train'], data['mean_x_val'], data['y_val'],
						**training_params[j])
		# remove extreme values
		# idx = np.where(data['mean_x_test'] > mu_min[k])[0]
		# idx = np.intersect1d(idx, np.where(data['mean_x_test'] < mu_max[k])[0])
		# tests
		tmp = rsp1.predict(data['mean_x_test']).squeeze()
		IVa_pval_global, IVa_tval_global, IVa_pval_nonlinear, IVa_tval_nonlinear, nonlinear_coef, nonlinear_se = stage2Tests(data['mean_x_test'], data['y_test'], tmp)
		# save model
		if True:
			rsp1.save("/home/panwei/he000176/deepRIV/UKB/models/MR/" + exposure + '/' + outcome + "_{}".format(k))
		with open(debug_path, "a") as f:
			f.write("{}: Done\n".format(k))
		#
		IVa_results = pd.DataFrame({"IVa_pval_global":IVa_pval_global, 
									"IVa_pval_nonlinear":IVa_pval_nonlinear,
									"IVa_tval_global":IVa_tval_global, 
									"IVa_tval_nonlinear":IVa_tval_nonlinear, 
									"nonlinear_coef":nonlinear_coef, 
									"nonlinear_se":nonlinear_se,
									"mu_max":mu_max, 
									"mu_min":mu_min}, index = [0])

		IVa_results.to_csv(IVa_path, mode = 'a', index = False, 
								sep = " ", header = not os.path.exists(IVa_path))

