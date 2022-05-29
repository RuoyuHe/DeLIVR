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
from tensorflow.keras.layers import Activation, Add, Input, Dense, ReLU, Reshape
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.optimizers import Adam
from sklearn.model_selection import StratifiedShuffleSplit
from tensorflow.keras.backend import expand_dims, squeeze, clear_session
from sklearn.model_selection import train_test_split
from wrappers import create_stage1, create_stage2, fit_stage2, tune_l2, stage2Tests
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt


r = robjects.r
r['source']('/home/panwei/he000176/deepRIV/UKB/code/hdl_stage1_FUN.R')

# Loading the function we have defined in R.
stage1_r = robjects.globalenv['stage1']

path = "/home/panwei/he000176/deepRIV/UKB/results/hdl/test40p_unrelated2/"
hdl = pd.read_csv(path + "hdl_combined_results.csv")
sig_idx = np.where(hdl['IVa_pval_nonlinear'] <= 0.05/hdl.shape[0])[0]
indices = hdl['gene_ind'].iloc[sig_idx].to_numpy()
gene_names = hdl['gene_name'].iloc[sig_idx].to_numpy()

k = 5
gene_ind = 15037
ind = 21

gene_ind = int(gene_ind)
print("gene_ind: {}".format(gene_ind))
#converting it into r object for passing into r function
with localconverter(robjects.default_converter + pandas2ri.converter):
	stage1_results = stage1_r(gene_ind)
	#
	mean_x_train = robjects.conversion.rpy2py(stage1_results[2])
	mean_x_val = robjects.conversion.rpy2py(stage1_results[3])
	mean_x_test = robjects.conversion.rpy2py(stage1_results[4])
	y_train = robjects.conversion.rpy2py(stage1_results[5])
	y_val = robjects.conversion.rpy2py(stage1_results[6])
	y_test = robjects.conversion.rpy2py(stage1_results[7])

lower_p = hdl['mu_min'].iloc[sig_idx[k]]
upper_p = hdl['mu_max'].iloc[sig_idx[k]]
#
unq1 = np.unique(mean_x_test)
unq1 = np.array([v for v in unq1 if v>=(lower_p - 1e12) and v<=(upper_p + 1e12)])
unq = []
mean_y = []
for i in range(unq1.shape[0]):
	idx = np.where(mean_x_test[:,0] == unq1[i])[0]
	if idx.shape[0] > 1000:
		unq.append(unq1[i])
		mean_y.append(np.mean(y_test[idx]))


unq = np.array(unq)
mean_y = np.array(mean_y)
#
new_x = np.linspace(np.min(unq), np.max(unq), 500)

l2 = 0
rsp1 = create_stage2((1,),1, l2 = l2)
h = fit_stage2(rsp1, mean_x_train, y_train, mean_x_val, y_val)
json.dump(h.history, open('../results/stage2_results/history_{}_default'.format(gene_ind), 'w'))
# h = json.load(open('../results/stage2_results/history7468_{}'.format(i), 'r'))

rsp2 = create_stage2((1,),1, l2 = l2)
h = fit_stage2(rsp2, mean_x_train, y_train, mean_x_val, y_val,
					learning_rate = 0.00001, 
					batch_size = 64, 
					training_steps = 8000,
					decay_steps = 150,
					decay_rate = 0.3,
					patience = 80)
json.dump(h.history, open('../results/stage2_results/history_{}_patience20_decay50'.format(gene_ind), 'w'))

tmp = rsp2.predict(mean_x_test).squeeze()
_, _, p_val, _, _, _ = stage2Tests(mean_x_test, y_test, tmp)
np.savetxt("../results/stage2_results/{}_patience20_decay50_pval.txt".format(gene_ind), [p_val])

rsp3 = create_stage2((1,),1, l2 = l2)
h = fit_stage2(rsp3, mean_x_train, y_train, mean_x_val, y_val,
					learning_rate = 0.00001, 
					batch_size = 64, 
					training_steps = 8000,
					patience = 300)
json.dump(h.history, open('../results/stage2_results/history_{}_patience300'.format(gene_ind), 'w'))


# comment out for new tests
if False:
	rsp1 = tf.keras.models.load_model("/home/panwei/he000176/deepRIV/UKB/models/hdl/sig_genes/{}_21".format(gene_ind))

	# before results
	pred_before = rsp1.predict(new_x).squeeze()
	tmp = rsp1.predict(mean_x_test).squeeze()
	mse_before = np.mean((y_test - tmp)**2)

	_ = fit_stage2(rsp1, mean_x_train, y_train, mean_x_val, y_val,
					learning_rate = 0.00001, 
					batch_size = 64, 
					training_steps = 8000,
					patience = 300)

	pred_after = rsp1.predict(new_x).squeeze()
	tmp = rsp1.predict(mean_x_test).squeeze()
	mse_after = np.mean((y_test - tmp)**2)
	_, _, p_val_after, _, _, _ = stage2Tests(mean_x_test, y_test, tmp)


	plt.xlim((np.min(unq)-0.05, np.max(unq)+0.05))
	plt.ylim((min(mean_y)-0.1, max(mean_y)+0.1))
	plt.plot(new_x, pred_before, color = "r", linestyle = "dashed", label = "before")
	plt.plot(new_x, pred_after, linestyle = "dashed", label = "after")
	plt.plot(unq, mean_y, color = "black", linestyle = "dashdot", label = r"$\bar{Y}_{\hat{\mu}}$")
	plt.scatter(mean_x_test, y_test, s=1, color = "lightgrey", label = r"$\hat{\mu}$ vs Y")
	plt.title("NFATC3")
	plt.legend(loc = "upper right")
	plt.savefig('/home/panwei/he000176/deepRIV/UKB/code/NFATC3_test.png')
	plt.clf()


	results = pd.DataFrame({"mse_before":mse_before, 
							"mse_after":mse_after, 
							"p_val_after":p_val_after
							}, index = [0]
							)

	results.to_csv("~/deepRIV/UKB/results/stage2_results/nfatc3_test.txt", index = False, sep = " ")

