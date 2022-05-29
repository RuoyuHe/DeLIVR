import os
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
from scipy.stats import binned_statistic
from matplotlib import pyplot as plt


exposures = pd.read_csv('~/deepRIV/UKB/data/MR_tmp/lm_bmi_adjusted_cov.csv').iloc[:,0]
outcomes = pd.read_csv('~/deepRIV/UKB/data/MR_tmp/normalized_outcomes_bmi.csv')


for l in range(1,6):
	exposure = 'bmi' 
	mean_X = exposures.to_numpy().reshape(-1,1)
	Y = outcomes.iloc[:,l].to_numpy().reshape(-1,1)
	#
	name = outcomes.columns[l]
	binned_mean, binned_edges, _ = binned_statistic(mean_X.squeeze(), Y.squeeze(), 
		statistic='mean', bins = 100)
	binned_edges = (binned_edges[:-1] + binned_edges[1:])/2
	#
	new_x = np.linspace(np.min(binned_edges), np.max(binned_edges), 500)
	pred = np.zeros((21,500))
	tmp_ind = 0
	for i in range(21):
		path = '/home/panwei/he000176/deepRIV/UKB/models/MR/bmi/' + name +'_{}'.format(i)
		if os.path.isdir(path):
			print(tmp_ind)
			rsp1 = tf.keras.models.load_model(path)
			pred[tmp_ind,:] = rsp1.predict(new_x).squeeze()
			tmp_ind += 1
	CI = np.quantile(pred, q = [0.025,0.975], axis = 0)
	mean_pred = np.mean(pred, axis = 0)
	#
	plt.xlim((np.min(binned_edges)-0.05, np.max(binned_edges)+0.05))
	plt.ylim((min(binned_mean)-0.1, max(binned_mean)+0.1))
	plt.xlabel(r"$\hat{\mu}_Z \ (BMI)$")
	plt.ylabel(name)
	plt.plot(new_x, mean_pred, color = "r", linestyle = "dashed", label = 'DeLIVR') #mean curve.
	plt.fill_between(new_x, CI[0,:], CI[1,:], color = 'r', alpha=.15) #std curves.
	plt.scatter(binned_edges, binned_mean, s=10, label = r"$\hat{\mu}_Z$ vs " + name)
	plt.title(exposure + ' vs ' + name)
	plt.legend(loc = "upper left")
	plt.savefig('/home/panwei/he000176/deepRIV/UKB/results/MR/' + exposure + '/' + name + '.png')
	plt.clf()

