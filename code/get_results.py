import os
import numpy as np
import pandas as pd

path = "/home/panwei/he000176/deepRIV/UKB/results/stage2_results/"

p_vals = pd.read_csv(path + "p_vals_1.txt", sep = " ")
for i in range(2, 133):
	if not os.path.exists(path + "p_vals_{}.txt".format(i)):
		continue
	tmp = pd.read_csv(path + "p_vals_{}.txt".format(i), sep = " ")
	p_vals = pd.concat([p_vals, tmp])

p_vals = p_vals.mask(p_vals.eq('None')).dropna()

alpha = 0.05
adj_alpha = alpha / p_vals.shape[0]

nonlinear_idx = np.where(p_vals['IVa_pval_nonlinear'].to_numpy() <= adj_alpha)[0]
global_idx = np.where(p_vals['IVa_pval_global'].to_numpy() <= adj_alpha)[0]
lr_idx = np.where(p_vals['LR_pval'].to_numpy() <= adj_alpha)[0]


# combine p
path = "/home/panwei/he000176/deepRIV/UKB/results/ldl/combine_p/"
common = pd.read_csv(path + "common_1.txt", sep = " ")
for i in range(2, 158):
	if not os.path.exists(path + "common_{}.txt".format(i)):
		continue
	tmp = pd.read_csv(path + "common_{}.txt".format(i), sep = " ")
	common = pd.concat([common, tmp])

common.to_csv(path + "common_combined.csv", index = False)

# IVa = pd.read_csv(path + "IVa_1.txt", sep = " ")
IVa = None
for i in range(158):
	if not os.path.exists(path + "IVa_{}.txt".format(i)):
		continue
	tmp = pd.read_csv(path + "IVa_{}.txt".format(i), sep = " ")
	if IVa is not None:
		IVa = pd.concat([IVa, tmp])
	else:
		IVa = tmp.copy()

colnames = IVa.columns.to_numpy()
colnames[0] = 'gene_names'
IVa.columns = colnames
IVa.to_csv(path + "IVa_combined.csv", index = False)

hdl = pd.read_csv("/home/panwei/he000176/deepRIV/UKB/results/hdl/test40p_unrelated2/hdl_combined_results.csv")

common = pd.read_csv(path + "common_combined.csv")
IVa = pd.read_csv(path + "IVa_combined.csv")
common_missed = pd.read_csv(path + "common_missed.txt", sep = " ")
IVa_missed = pd.read_csv(path + "IVa_missed.txt", sep = " ")

tmp = np.in1d(common_missed['gene_names'], common['gene_names'])
idx = np.where(tmp==False)[0]
common_missed = common_missed.iloc[idx,:]
common = pd.concat([common, common_missed])
common = common.set_index(pd.Index(np.arange(4701)))
idx = common.reset_index().set_index('gene_names').loc[hdl.gene_names, 'index'].values
common = common.iloc[idx,:]
common = common.set_index(pd.Index(np.arange(4701)))
common['gene_name'] = hdl['gene_name']
common['gene_ind'] = hdl['gene_ind']
common['chr'] = hdl['chr']
common.to_csv(path + "common_combined.csv", index = False)

tmp = np.in1d(IVa_missed['gene_names'], IVa['gene_names'])
idx = np.where(tmp==False)[0]
IVa_missed = IVa_missed.iloc[idx,:]
IVa = pd.concat([IVa, IVa_missed])
IVa = IVa.set_index(pd.Index(np.arange(4701)))
idx = IVa.reset_index().set_index('gene_names').loc[hdl.gene_names, 'index'].values
IVa = IVa.iloc[idx,:]
IVa = IVa.set_index(pd.Index(np.arange(4701)))
IVa['gene_name'] = hdl['gene_name']
IVa['gene_ind'] = hdl['gene_ind']
IVa['chr'] = hdl['chr']
IVa.to_csv(path + "IVa_combined.csv", index = False)


# cv_test
path = "/home/panwei/he000176/deepRIV/UKB/results/hdl/combine_p/"
cv_test = pd.read_csv(path + "cv_test_0.txt", sep = " ")
for i in range(1, 10):
	if not os.path.exists(path + "cv_test_{}.txt".format(i)):
		continue
	tmp = pd.read_csv(path + "cv_test_{}.txt".format(i), sep = " ")
	cv_test = pd.concat([cv_test, tmp])

cv_test.to_csv(path + "cv_test_combined.csv", index = False)

cv_test_IVa = pd.read_csv(path + "cv_test_IVa_0.txt", sep = " ")
for i in range(1, 10):
	if not os.path.exists(path + "cv_test_IVa_{}.txt".format(i)):
		continue
	tmp = pd.read_csv(path + "cv_test_IVa_{}.txt".format(i), sep = " ")
	cv_test_IVa = pd.concat([cv_test_IVa, tmp])

cv_test_IVa.to_csv(path + "cv_test_IVa_combined.csv", index = False)





missed_debug = []
missed_common = []
missed_iva = []
for i in range(1,283):
	if not os.path.exists("common_{}.txt".format(i)):
		missed_common.append(i)
	if not os.path.exists("IVa_{}.txt".format(i)):
		missed_iva.append(i)
	if not os.path.exists("/home/panwei/he000176/deepRIV/UKB/results/debug/imputed_combine_p/{}.txt".format(i)):
		missed_debug.append(i)


import numpy as np, pandas as pd
a = np.arange(10)
a = np.append(a,10892)
tmp = pd.DataFrame()
num_repeats = 21
num_genes = 11

for i in a:
	print(tmp.shape, i)
	tmp = tmp.append(pd.read_csv('IV_{}.txt'.format(i), sep = ' '))


regrouped_global_pvals = tmp.IV_pval_global.to_numpy().reshape(-1,21)
regrouped_nl_pvals = tmp.IV_pval_nonlinear.to_numpy().reshape(-1,21)
cauchy_global_pvals = cauchy_p(regrouped_global_pvals, num_repeats)
cauchy_nl_pvals = cauchy_p(regrouped_nl_pvals, num_repeats)

gene_names = tmp.gene_names.to_numpy().reshape(-1,num_repeats)[:,0]

deepiv_hdl_cauchy = pd.DataFrame({'gene_names': gene_names,
									'IV_global_pvals':cauchy_global_pvals,
									'IV_nl_pvals':cauchy_nl_pvals})
deepiv_hdl_cauchy.to_csv('~/deepRIV/UKB/results/hdl/combine_p/repeat_deepiv/hdl_deepiv_cauchy_results.csv',
						index = False)



import numpy as np, pandas as pd

num_repeats = 2
setup = ['typeI', 'linear', 'quadratic', 'cos']
cutoff = 0.05

path = '~/deepRIV/UKB/simulation/results/cdk2ap1/'
cdk2ap1_results = pd.DataFrame()
for l in range(len(setup)):
	tmp = pd.DataFrame()
	for i in range(5):
		print(tmp.shape, i)
		tmp = tmp.append(pd.read_csv(path + 'delivr_' + setup[l] + '_{}.txt'.format(i), sep = ' '))
	regrouped_global_pvals = tmp.IVa_pval_global.to_numpy().reshape(-1,num_repeats)
	regrouped_nl_pvals = tmp.IVa_pval_nonlinear.to_numpy().reshape(-1,num_repeats)
	cauchy_global_pvals = cauchy_p(regrouped_global_pvals, num_repeats)
	cauchy_nl_pvals = cauchy_p(regrouped_nl_pvals, num_repeats)
	power_global = np.sum(cauchy_global_pvals <= cutoff)/50
	power_nl = np.sum(cauchy_nl_pvals <= cutoff)/50
	power_l = np.sum(tmp.twal_pval <= cutoff)/100
	power_lq = np.sum(tmp.twaslq_pval <= cutoff)/100
	cdk2ap1_results = cdk2ap1_results.append([[power_l, power_lq, power_global, power_nl]])

cdk2ap1_results.columns = ['TWAS-L', 'TWAS-LQ', 'Global', 'Nonlinear']

path = '~/deepRIV/UKB/simulation/results/nt5dc2/'
nt5dc2_results = pd.DataFrame()
for l in range(len(setup)):
	tmp = pd.DataFrame()
	for i in range(5):
		print(tmp.shape, i)
		tmp = tmp.append(pd.read_csv(path + 'delivr_' + setup[l] + '_{}.txt'.format(i), sep = ' '))
	regrouped_global_pvals = tmp.IVa_pval_global.to_numpy().reshape(-1,num_repeats)
	regrouped_nl_pvals = tmp.IVa_pval_nonlinear.to_numpy().reshape(-1,num_repeats)
	cauchy_global_pvals = cauchy_p(regrouped_global_pvals, num_repeats)
	cauchy_nl_pvals = cauchy_p(regrouped_nl_pvals, num_repeats)
	power_global = np.sum(cauchy_global_pvals <= cutoff)/50
	power_nl = np.sum(cauchy_nl_pvals <= cutoff)/50
	power_l = np.sum(tmp.twal_pval <= cutoff)/100
	power_lq = np.sum(tmp.twaslq_pval <= cutoff)/100
	nt5dc2_results = nt5dc2_results.append([[power_l, power_lq, power_global, power_nl]])

nt5dc2_results.columns = ['TWAS-L', 'TWAS-LQ', 'Global', 'Nonlinear']
