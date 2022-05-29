# CI
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
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from matplotlib import pyplot as plt

r = robjects.r
r['source']('/home/panwei/he000176/deepRIV/UKB/code/hdl/combine_p/hdl_stage1_FUN.R')
# Loading the function we have defined in R.
stage1_r = robjects.globalenv['stage1']


path = "/home/panwei/he000176/deepRIV/UKB/results/hdl/test40p_unrelated2/"
hdl = pd.read_csv(path + "hdl_combined_results.csv")
# sig_idx = np.where(hdl['IVa_pval_nonlinear'] <= 0.05/hdl.shape[0])[0]
# indices = hdl['gene_ind'].iloc[sig_idx].to_numpy()
# gene_names = hdl['gene_name'].iloc[sig_idx].to_numpy()

gene_name = 'NT5DC2'
indices = hdl.loc[hdl.gene_name==gene_name,'gene_ind'].to_numpy()
gene_names = hdl['gene_name'].iloc[sig_idx].to_numpy()

runs = 30
num_params = 3
num_repeats = 7

for k, gene_ind in enumerate(indices):
	gene_ind = int(gene_ind)
	print("gene_ind: {}".format(gene_ind))
	#converting it into r object for passing into r function
	with localconverter(robjects.default_converter + pandas2ri.converter):
		stage1_results = stage1_r(gene_ind)
		#
		mean_X = robjects.conversion.rpy2py(stage1_results[2])
		Y = robjects.conversion.rpy2py(stage1_results[3])
	#
	lower_p = hdl.loc[hdl.gene_name==gene_name,'mu_min'].values[0]
	upper_p = hdl.loc[hdl.gene_name==gene_name,'mu_max'].values[0]
	#
	unq1 = np.unique(mean_X)
	unq1 = np.array([v for v in unq1 if v>=(lower_p - 1e12) and v<=(upper_p + 1e12)])
	unq = []
	mean_y = []
	for i in range(unq1.shape[0]):
		idx = np.where(mean_X[:,0] == unq1[i])[0]
		if idx.shape[0] > 1000:
			unq.append(unq1[i])
			mean_y.append(np.mean(Y[idx]))
	unq = np.array(unq)
	mean_y = np.array(mean_y)
	#
	new_x = np.linspace(np.min(unq), np.max(unq), 500)
	#
	tmp = np.zeros((runs,mean_X.shape[0]))
	pred = np.zeros((runs,500))
	for j in range(num_params):
		for k in range(j*num_repeats, (j+1)*num_repeats):
			rsp1 = tf.keras.models.load_model("/home/panwei/he000176/deepRIV/UKB/models/hdl/combine_p/model_{}_{}_{}".format(gene_ind,j,k))
			# tmp[i,:] = rsp1.predict(mean_x_test).squeeze()
			pred[i,:] = rsp1.predict(new_x).squeeze()
			# design = sm.add_constant(tmp[i,:].reshape(-1,1))
			# eig_min = np.min(np.linalg.eig(design.T @ design)[0])
			# if eig_min < 1:
			# 	print("small eig value! {}, {}, {}".format(gene_ind, i, eig_min))
	CI = np.quantile(pred, q = [0.025,0.975], axis = 0)
	mean_pred = np.mean(pred, axis = 0)
	#
	plt.xlim((np.min(unq)-0.05, np.max(unq)+0.05))
	plt.ylim((min(mean_y)-0.1, max(mean_y)+0.1))
	plt.plot(new_x, mean_pred, color = "r", linestyle = "dashed", label = gene_names[k]) #mean curve.
	plt.fill_between(new_x, CI[0,:], CI[1,:], color = 'r', alpha=.15) #std curves.
	plt.plot(unq, mean_y, color = "black", linestyle = "dashdot", label = r"$\bar{Y}_{\hat{\mu}}$")
	plt.scatter(mean_x_test, y_test, s=1, color = "lightgrey", label = r"$\hat{\mu}$ vs Y")
	plt.title(gene_names[k])
	plt.legend(loc = "upper right")
	plt.savefig('/home/panwei/he000176/deepRIV/UKB/code/' + gene_names[k] + '.png')
	plt.clf()







# CI delivr
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
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from matplotlib import pyplot as plt

r = robjects.r
r['source']('/home/panwei/he000176/deepRIV/UKB/code/hdl/combine_p/hdl_stage1_FUN.R')
# r['source']('/home/panwei/he000176/deepRIV/UKB/code/ldl/combine_p/ldl_stage1_FUN.R')

# Loading the function we have defined in R.
stage1_r = robjects.globalenv['stage1']


path = "/home/panwei/he000176/deepRIV/UKB/results/hdl/test40p_unrelated2/"
hdl = pd.read_csv(path + "hdl_combined_results.csv")
sig_idx = np.where(hdl['IVa_pval_nonlinear'] <= 0.05/hdl.shape[0])[0]
indices = hdl['gene_ind'].iloc[sig_idx].to_numpy()

runs = 21
gene_inds = [12441, 3834, 11276, 15100, 11286, 
				17467, 17471, 17475, 17199] # 11286:AP000892.6 17467: MAU2 17471:YJEFN3 17475: GMIP 17199: SLC44A2

gene_ind = 11276

gene_ind = int(gene_ind)
print("gene_ind: {}".format(gene_ind))
#converting it into r object for passing into r function
with localconverter(robjects.default_converter + pandas2ri.converter):
	stage1_results = stage1_r(gene_ind)
	gene_name = robjects.conversion.rpy2py(stage1_results[0])
	mean_X = robjects.conversion.rpy2py(stage1_results[2])
	Y = robjects.conversion.rpy2py(stage1_results[3])
#

r = robjects.r
r['source']('/home/panwei/he000176/deepRIV/UKB/code/TWAS_LQ_Q/hdl_twaslq_q.R')
# r['source']('/home/panwei/he000176/deepRIV/UKB/code/TWAS_LQ_Q/ldl_twaslq_q.R')

# Loading the function we have defined in R.
twaslq_r = robjects.globalenv['stage1']
with localconverter(robjects.default_converter + pandas2ri.converter):
	twaslq_results = twaslq_r(gene_ind)
	twaslq_coefs = robjects.conversion.rpy2py(twaslq_results[21])
#

# for AP000892.6
# twaslq_coefs = [-0.01057596, 0.06615957]
# for MAU2
# twaslq_coefs = [0.106062, -0.06130475]

tmp = np.zeros((runs,mean_X.shape[0]))
lower_p = hdl.loc[hdl.gene_ind==gene_ind,'mu_min'].values
upper_p = hdl.loc[hdl.gene_ind==gene_ind,'mu_max'].values
#
unq1 = np.unique(mean_X)
unq1 = np.array([v for v in unq1 if v>=(lower_p - 1e12) and v<=(upper_p + 1e12)])
unq = []
mean_y = []
for i in range(unq1.shape[0]):
	idx = np.where(mean_X[:,0] == unq1[i])[0]
	if idx.shape[0] > 1000:
		unq.append(unq1[i])
		mean_y.append(np.mean(Y[idx]))

unq = np.array(unq)
mean_y = np.array(mean_y)
#
new_x = np.linspace(np.min(unq), np.max(unq), 500)
num_repeats = 7
#
k = 0
pred = np.zeros((runs,500))
pred_robust = np.zeros((runs,500))
for i in range(3):
	for j in range(i*num_repeats, (i+1)*num_repeats):
		# path = "/home/panwei/he000176/deepRIV/UKB/models/hdl/combine_p/tmp_for_models/model_test_{}_{}_{}".format(gene_ind,i,j)
		path = "/home/panwei/he000176/deepRIV/UKB/models/hdl/combine_p//model_{}_{}_{}".format(gene_ind,i,j)
		robust_path = "/home/panwei/he000176/deepRIV/UKB/models/ldl/robust_DeLIVR/model_{}_{}_{}".format(gene_ind,i,j)
		if os.path.isdir(path):
			rsp1 = tf.keras.models.load_model(path)
			# pred[k,:] = rsp1.predict(new_x).squeeze()
			pred[k,:] = rsp1.predict((new_x-np.mean(mean_X))/np.std(mean_X)).squeeze()
			# rsp_robust = tf.keras.models.load_model(robust_path)
			# pred_robust[k,:] = rsp_robust.predict((new_x-np.mean(mean_X))/np.std(mean_X)).squeeze()
			k += 1

CI = np.quantile(pred, q = [0.025,0.975], axis = 0)
mean_pred = np.mean(pred, axis = 0)
CI_robust = np.quantile(pred_robust, q = [0.025,0.975], axis = 0)
mean_pred_robust = np.mean(pred_robust, axis = 0)
twaslq_pred = twaslq_coefs[-2]*new_x + twaslq_coefs[-1]*new_x**2
#
gene_name = hdl.loc[hdl.gene_ind==gene_ind,'gene_name'].values[0]
plt.xlim((np.min(unq)-0.05, np.max(unq)+0.05))
plt.ylim((min(mean_y)-0.1, max(mean_y)+0.1))
plt.xlabel(r"$\mu_Z$ or $X$")
plt.plot(new_x, mean_pred, color = "r", linestyle = "dashed", label = 'DeLIVR') #mean curve.
plt.fill_between(new_x, CI[0,:], CI[1,:], color = 'r', alpha=.15) #std curves.
# plt.plot(new_x, mean_pred_robust, color = "#549955", linestyle = "dashdot", label = 'rDeLIVR') #mean curve.
# plt.fill_between(new_x, CI_robust[0,:], CI_robust[1,:], color = '#549955', alpha=.2) #std curves.
plt.plot(new_x, twaslq_pred , color = "blue", linestyle = "dotted", label = 'TWAS-LQ') #mean curve.
# plt.plot(unq, mean_y, color = "black", linestyle = "dashdot", label = r"$\bar{Y}_{\hat{\mu}}$")
plt.scatter(mean_X, Y, s=1, color = "lightgrey", label = r"$\hat{\mu}$ vs Y")
plt.title(gene_name)
plt.legend(loc = "lower right")
plt.savefig('/home/panwei/he000176/deepRIV/UKB/' + gene_name + '_twas_delivr' + '.pdf')
plt.clf()




# CI simulation delivr
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
from matplotlib import pyplot as plt


path = "/home/panwei/he000176/deepRIV/UKB/results/hdl/test40p_unrelated2/"
hdl = pd.read_csv(path + "hdl_combined_results.csv")
sig_idx = np.where(hdl['IVa_pval_nonlinear'] <= 0.05/hdl.shape[0])[0]
indices = hdl['gene_ind'].iloc[sig_idx].to_numpy()

gene_ind = 3834
setup = ['Null', 'linear', 'quadratic', 'cubic']
gene_name = hdl.loc[hdl.gene_ind==gene_ind,'gene_name'].values[0]

for l in range(len(setup)-1):
	if l==0:
		true_model = 'typeI'
	else:
		true_model = setup[l]
	out_path = '/home/panwei/he000176/deepRIV/UKB/simulation/'
	#
	z_path = out_path + 'data/' + gene_name.lower() +'/ukb_snp_bed.csv'
	beta_path = out_path + 'data/' + gene_name.lower() +'/hatbetaX1.txt'
	theta_path = out_path + 'data/' + gene_name.lower() +'/hattheta_lq.txt'
	#
	z = pd.read_csv(z_path).to_numpy()
	beta = pd.read_csv(beta_path, sep = ' ').to_numpy()
	nonzero_idx = np.where(beta.squeeze()!=0)[0]
	z = z[:, nonzero_idx]
	beta = beta[nonzero_idx, :]
	theta = pd.read_csv(theta_path, sep = ' ').to_numpy()
	#
	gamma = 0.7
	if l==0:
		true_g = lambda x: 0
		true_g_mean = lambda x: 0*x
	elif l==1:
		true_g = lambda x: x*theta[0,0]
		true_g_mean = lambda x: x*theta[0,0]
	elif l==2:
		true_g = lambda x: x*theta[0,0] + x**2 * theta[1,0]
		true_g_mean = lambda x: x*theta[0,0] + (x**2 + 1) * theta[1,0]
	elif l==3:
		true_g = lambda x: -0.03*x**3 + 0.025*x**2 + 0.1*x
		true_g_mean = lambda x: -0.03*(x**3 + 3*x) + 0.025*(x**2 + 1) + 0.1*x
	#
	e = np.random.multivariate_normal(mean = [0,0], cov = [[1,gamma],[gamma,1]], size = z.shape[0])
	x = (z@beta).squeeze() + e[:,0]
	y = true_g(x) + e[:,1]
	x = x.reshape(-1,1)
	#
	design = sm.add_constant(z)
	stage1 = sm.OLS(x, design)
	stage1 = stage1.fit()
	mean_X = stage1.predict(design).reshape(-1,1)
	Y = y.copy()
	#
	lower_p = hdl.loc[hdl.gene_ind==gene_ind,'mu_min'].values
	upper_p = hdl.loc[hdl.gene_ind==gene_ind,'mu_max'].values
	#
	unq1 = np.unique(mean_X)
	unq1 = np.array([v for v in unq1 if v>=(lower_p - 1e12) and v<=(upper_p + 1e12)])
	unq = []
	mean_y = []
	for i in range(unq1.shape[0]):
		idx = np.where(mean_X[:,0] == unq1[i])[0]
		if idx.shape[0] > 1000:
			unq.append(unq1[i])
			mean_y.append(np.mean(Y[idx]))
	unq = np.array(unq)
	mean_y = np.array(mean_y)
	#
	new_x = np.linspace(np.min(unq), np.max(unq), 500)
	#
	tmp_ind = 0
	# pred = np.zeros((100,500))
	# for file_num in range(5):
	# 	for i in range(10):
	# 		for j in range(2):
	# 			path = '/home/panwei/he000176/deepRIV/UKB/simulation/models/' + gene_name.lower() +'/delivr_' + true_model + '_{}_{}_{}'.format(file_num,i,j)
	# 			if os.path.isdir(path):
	# 				print(tmp_ind)
	# 				rsp1 = tf.keras.models.load_model(path)
	# 				pred[tmp_ind,:] = rsp1.predict(new_x).squeeze()
	# 				tmp_ind += 1
	pred = np.zeros((20,500))
	for i in range(10):
		for j in range(2):
			path = '/home/panwei/he000176/deepRIV/UKB/simulation/models/' + gene_name.lower() +'/delivr_' + true_model + '_{}_{}'.format(i,j)
			if os.path.isdir(path):
				print(tmp_ind)
				rsp1 = tf.keras.models.load_model(path)
				pred[tmp_ind,:] = rsp1.predict(new_x).squeeze()
				tmp_ind += 1
	CI = np.quantile(pred, q = [0.025,0.975], axis = 0)
	mean_pred = np.mean(pred, axis = 0)
	#
	plt.xlim((np.min(unq)-0.05, np.max(unq)+0.05))
	plt.ylim((min(mean_y)-0.1, max(mean_y)+0.1))
	plt.xlabel(r"$\mu_Z$")
	plt.plot(new_x, mean_pred, color = "r", linestyle = "dashed", label = 'DeLIVR') #mean curve.
	plt.fill_between(new_x, CI[0,:], CI[1,:], color = 'r', alpha=.15) #std curves.
	# plt.plot(unq, mean_y, color = "blue", linestyle = "dashdot", label = r"$\bar{Y}_{\hat{\mu}}$")
	plt.plot(new_x, true_g_mean(new_x), color = "black", linestyle = "dashdot", label = r"$E(g(X)|Z)$")
	plt.scatter(mean_X, Y, s=1, color = "lightgrey", label = r"$\hat{\mu}$ vs Y")
	plt.title(gene_name + ' ' + setup[l])
	plt.legend(loc = "upper right")
	plt.savefig('/home/panwei/he000176/deepRIV/UKB/simulation/' + gene_name + '_delivr_' + true_model + '.png')
	plt.clf()






# CI simulation deepiv
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
from matplotlib import pyplot as plt

@tf.function
def evaluate(model, inputs, **kwargs):
	return model(inputs, **kwargs)


path = "/home/panwei/he000176/deepRIV/UKB/results/hdl/test40p_unrelated2/"
hdl = pd.read_csv(path + "hdl_combined_results.csv")
sig_idx = np.where(hdl['IVa_pval_nonlinear'] <= 0.05/hdl.shape[0])[0]
indices = hdl['gene_ind'].iloc[sig_idx].to_numpy()

gene_ind = 3834
setup = ['Null', 'linear', 'quadratic', 'cubic']
gene_name = hdl.loc[hdl.gene_ind==gene_ind,'gene_name'].values[0]

for l in range(1, len(setup)):
	if l==0:
		true_model = 'typeI'
	else:
		true_model = setup[l]
	out_path = '/home/panwei/he000176/deepRIV/UKB/simulation/'
	#
	z_path = out_path + 'data/' + gene_name.lower() +'/ukb_snp_bed.csv'
	beta_path = out_path + 'data/' + gene_name.lower() +'/hatbetaX1.txt'
	theta_path = out_path + 'data/' + gene_name.lower() +'/hattheta_lq.txt'
	#
	z = pd.read_csv(z_path).to_numpy()
	beta = pd.read_csv(beta_path, sep = ' ').to_numpy()
	nonzero_idx = np.where(beta.squeeze()!=0)[0]
	z = z[:, nonzero_idx]
	beta = beta[nonzero_idx, :]
	theta = pd.read_csv(theta_path, sep = ' ').to_numpy()
	#
	gamma = 0.7
	if l==0:
		true_g = lambda x: 0*x
		true_g_mean = lambda x: 0*x
	elif l==1:
		true_g = lambda x: x*theta[0,0]
		true_g_mean = lambda x: x*theta[0,0]
	elif l==2:
		true_g = lambda x: x*theta[0,0] + x**2 * theta[1,0]
		true_g_mean = lambda x: x*theta[0,0] + (x**2 + 1) * theta[1,0]
	elif l==3:
		true_g = lambda x: -0.03*x**3 + 0.025*x**2 + 0.1*x
		true_g_mean = lambda x: -0.03*(x**3 + 3*x) + 0.025*(x**2 + 1) + 0.1*x
	#
	e = np.random.multivariate_normal(mean = [0,0], cov = [[1,gamma],[gamma,1]], size = z.shape[0])
	x = (z@beta).squeeze() + e[:,0]
	y = true_g(x) + e[:,1]
	x = x.reshape(-1,1)
	#
	design = sm.add_constant(z)
	stage1 = sm.OLS(x, design)
	stage1 = stage1.fit()
	mean_X = stage1.predict(design).reshape(-1,1)
	Y = y.copy()
	#
	lower_p = hdl.loc[hdl.gene_ind==gene_ind,'mu_min'].values
	upper_p = hdl.loc[hdl.gene_ind==gene_ind,'mu_max'].values
	#
	unq1 = np.unique(mean_X)
	unq1 = np.array([v for v in unq1 if v>=(lower_p - 1e12) and v<=(upper_p + 1e12)])
	unq = []
	mean_y = []
	for i in range(unq1.shape[0]):
		idx = np.where(mean_X[:,0] == unq1[i])[0]
		if idx.shape[0] > 1000:
			unq.append(unq1[i])
			mean_y.append(np.mean(Y[idx]))
	unq = np.array(unq)
	mean_y = np.array(mean_y)
	#
	n_draws = 16
	new_x = np.linspace(np.min(unq), np.max(unq), 500)
	new_x2 = np.linspace(np.min(x), np.max(x), 500)
	pred_x = np.random.normal(new_x.reshape(-1,1), 1, (new_x.shape[0],n_draws))
	#
	tmp_ind = 0
	pred = np.zeros((100,500))
	pred_mean = np.zeros((100,500))
	for file_num in range(5):
		for i in range(10):
			for j in range(2):
				path = '/home/panwei/he000176/deepRIV/UKB/simulation/models/' + gene_name.lower() +'/deepiv_' + true_model + '_{}_{}_{}'.format(file_num,i,j)
				if os.path.isdir(path):
					print(tmp_ind)
					rsp1 = tf.keras.models.load_model(path)
					pred[tmp_ind,:] = evaluate(rsp1, new_x.reshape(-1,1)).numpy().squeeze()
					tmp_pred = evaluate(rsp1, pred_x.reshape(-1,1))
					pred_mean[tmp_ind,:] = tf.reduce_mean(tf.reshape(tmp_pred, (-1, n_draws)) ,axis = 1)
					tmp_ind += 1
	CI = np.quantile(pred_mean[~(pred_mean==0).all(1)], q = [0.025,0.975], axis = 0)
	mean_pred = np.mean(pred_mean[~(pred_mean==0).all(1)], axis = 0)
	#
	plt.xlim((np.min(unq)-0.05, np.max(unq)+0.05))
	plt.ylim((min(mean_y)-0.1, max(mean_y)+0.1))
	plt.xlabel(r"$\mu_Z$")
	# plt.ylabel(r"$y$")
	plt.plot(new_x, mean_pred, color = "r", linestyle = "dashed", label = 'DeepIV') #mean curve.
	plt.fill_between(new_x, CI[0,:], CI[1,:], color = 'r', alpha=.15) #std curves.
	# plt.plot(unq, mean_y, color = "black", linestyle = "dashdot", label = r"$\bar{Y}_{\hat{\mu}}$")
	plt.plot(new_x, true_g_mean(new_x), color = "black", linestyle = "dashdot", label = r"$E(g(X)|Z)$")
	plt.scatter(mean_X, Y, s=1, color = "lightgrey", label = r"$\hat{\mu}$ vs Y")
	plt.title(gene_name + ' ' + setup[l])
	plt.legend(loc = "upper right")
	plt.savefig('/home/panwei/he000176/deepRIV/UKB/simulation/' + gene_name + '_deepiv_' + true_model + '_mean.png')
	plt.clf()
	#
	CI = np.quantile(pred[~(pred==0).all(1)], q = [0.025,0.975], axis = 0)
	mean_pred = np.mean(pred[~(pred==0).all(1)], axis = 0)
	#
	plt.xlim((np.min(unq)-0.05, np.max(unq)+0.05))
	plt.ylim((min(mean_y)-0.1, max(mean_y)+0.1))
	plt.xlabel(r"$x$")
	# plt.ylabel(r"$y$")
	plt.plot(new_x, mean_pred, color = "r", linestyle = "dashed", label = 'DeepIV') #mean curve.
	plt.fill_between(new_x, CI[0,:], CI[1,:], color = 'r', alpha=.15) #std curves.
	# plt.plot(unq, mean_y, color = "black", linestyle = "dashdot", label = r"$\bar{Y}_{\hat{\mu}}$")
	plt.plot(new_x, true_g(new_x), color = "black", linestyle = "dashdot", label = r"$g(X)$")
	plt.scatter(mean_X, Y, s=1, color = "lightgrey", label = r"$\hat{\mu}$ vs Y")
	plt.title(gene_name + ' ' + setup[l])
	plt.legend(loc = "upper right")
	plt.savefig('/home/panwei/he000176/deepRIV/UKB/simulation/' + gene_name + '_deepiv_' + true_model + '.png')
	plt.clf()




