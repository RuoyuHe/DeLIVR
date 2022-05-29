import os
import itertools
import numpy as np
import pandas as pd

path = "/home/panwei/he000176/deepRIV/UKB/results/hdl/combine_p/no_ind_test/"

common = pd.read_csv(path + "common_combined.csv")
IVa = pd.read_csv(path + "IVa_combined.csv")
hdl = pd.read_csv("/home/panwei/he000176/deepRIV/UKB/results/hdl/test40p_unrelated2/hdl_combined_results.csv")




num_repeats = 3
colnames = [['gene_names'], ['global_p']*num_repeats, ['nonlinear_p']*num_repeats,
				['global_t']*num_repeats, ['nonlinear_t']*num_repeats,
				['nonlinear_coef']*num_repeats, ['nonlinear_se']*num_repeats,
				['mu_max']*num_repeats, ['mu_min']*num_repeats,
				['gene_name'], ['gene_ind']]
colnames = list(itertools.chain(*colnames))
#colnames.extend(list(IVa.columns[-num_repeats:]))
IVa.columns = colnames
IVa['rsq'] = common['rsq']
IVa['adj_rsq'] = common['adj_rsq']

cutoff = 0.05/common.shape[0]

IVa = IVa.loc[(IVa.nonlinear_p != 100.0).all(axis=1),:]

# Cauchy
cauchy_G_T = np.tan((0.5-IVa.global_p)*np.pi).sum(axis = 1) / num_repeats
cauchy_NL_T = np.tan((0.5-IVa.nonlinear_p)*np.pi).sum(axis = 1) / num_repeats

cauchy_G_P = 1/2 - np.arctan(cauchy_G_T)/np.pi
cauchy_NL_P = 1/2 - np.arctan(cauchy_NL_T)/np.pi

cauchy_G_results = pd.concat([IVa.loc[cauchy_G_P <= cutoff, 'gene_name'], 
	cauchy_G_P.loc[cauchy_G_P <= cutoff]], axis=1)
cauchy_NL_results = pd.concat([IVa.loc[cauchy_NL_P <= cutoff, 'gene_name'], 
	cauchy_NL_P.loc[cauchy_NL_P <= cutoff]], axis=1)

cauchy_LR_inter = pd.Series(list(set(hdl.loc[hdl.LR_pval<=cutoff, 'gene_name']).intersection(set(cauchy_G_results.gene_name))))
cauchy_G_NL_inter = pd.Series(list(set(cauchy_NL_results.gene_name).intersection(set(cauchy_G_results.gene_name))))


# Hommel's method
Cn = np.sum([1/i for i in range(1,num_repeats+1)])
hommel_cutoffs = np.array([k*cutoff / (Cn*num_repeats) for k in range(1, num_repeats+1)])

G_sorted = np.sort(IVa.global_p.to_numpy(), axis = 1)
NL_sorted = np.sort(IVa.nonlinear_p.to_numpy(), axis = 1)

G_idx = np.where(np.sum(G_sorted <= hommel_cutoffs, axis = 1))[0]
NL_idx = np.where(np.sum(NL_sorted <= hommel_cutoffs, axis = 1))[0]

hommel_G_results = IVa.gene_name.iloc[G_idx]
hommel_NL_results = IVa.gene_name.iloc[NL_idx]

hommel_LR_inter = pd.Series(list(set(hommel_G_results.values).intersection(set(hdl.loc[hdl.LR_pval<=cutoff, 'gene_name']))))
hommel_G_NL_inter = pd.Series(list(set(hommel_NL_results.values).intersection(set(hommel_G_results.values))))


# cv_test
from scipy.stats import chi2

df = 2*2

cv_test = pd.read_csv(path + "cv_test_combined.csv")
cv_test_IVa = pd.read_csv(path + "cv_test_IVa_combined.csv")

# cv cauchy
cv_original_p = cv_test.iloc[:,1:3].values
cv_cauchy_T = np.tan((0.5-cv_original_p)*np.pi).sum(axis = 1) / 2
cv_cauchy_P = 1/2 - np.arctan(cv_cauchy_T)/np.pi



# Fisher's method + cauchy
a = cv_test_IVa.global_p.values.reshape(-1,2)
b = cv_test_IVa['global_p.1'].values.reshape(-1,2)

a = -2*np.sum(np.log(a), axis = 1)
b = -2*np.sum(np.log(b), axis = 1)

a = 1 - chi2.cdf(a, df)
b = 1 - chi2.cdf(b, df)

cv_IVa_fcauchy_T = (np.tan((0.5-a)*np.pi) + np.tan((0.5-b)*np.pi)) / 2
cv_IVa_fcauchy_P = 1/2 - np.arctan(cv_IVa_fcauchy_T)/np.pi

# Cauchy
original_p = cv_test_IVa.iloc[:,1:3].values.reshape(-1,4)
cv_IVa_cauchy_T = np.tan((0.5-original_p)*np.pi).sum(axis = 1) / 4
cv_IVa_cauchy_P = 1/2 - np.arctan(cv_IVa_cauchy_T)/np.pi



#
# tmp = common.loc[common.LR_pval <= cutoff, ['gene_names','LR_pval','adj_rsq','gene_name','gene_ind']]
# to_exclude = IVa.loc[IVa.global_p.min(axis=1) <= cutoff, 'gene_names']
# tmp = tmp.loc[~tmp.index.isin(to_exclude.index)]

# close_call = IVa.loc[IVa.index.isin(tmp.index), ['gene_names', 'global_p', 'adj_rsq', 'gene_name', 'gene_ind']]
# close_call = close_call.sort_values(by = ['adj_rsq'])
# close_call.gene_ind.iloc[[:5,-5:]]

# target = pd.concat([close_call.head(5), close_call.tail(5)])
# target.to_csv("~/deepRIV/UKB/code/hdl/cv_p_test/target.csv")




repeat_NL = IVa.loc[cauchy_NL_P <= 1e-4, ['gene_name', 'gene_ind']]
repeat_NL2 = IVa.loc[(cauchy_NL_P <= 1e-3)&(IVa.adj_rsq >= 0.16), ['gene_name', 'gene_ind']]
idx = np.where(np.in1d(repeat_NL2.gene_name.values, repeat_NL.gene_name.values)==False)[0]

repeat_NL = pd.concat([repeat_NL, repeat_NL2.iloc[idx,:]])
repeat_NL = repeat_NL.sort_values(by = ['gene_ind'])


repeat_G = hdl.loc[hdl.LR_pval<=cutoff, ['gene_name','gene_ind']]
idx = np.where(np.in1d(repeat_G.gene_name.values, IVa.loc[cauchy_G_P <= cutoff, 'gene_name'].values)==False)[0]
repeat_G = repeat_G.iloc[idx,:]
repeat_G2 = hdl.loc[(hdl.LR_pval <= 1e-3)&(hdl.adj_rsq >= 0.5), ['gene_name', 'gene_ind']]
idx = np.where(np.in1d(repeat_G2.gene_name.values, repeat_G.gene_name.values)==False)[0]

repeat_G = pd.concat([repeat_G, repeat_G2.iloc[idx,:]])
repeat_G = repeat_G.sort_values(by = ['gene_ind'])

idx = np.where(np.in1d(repeat_NL.gene_name.values, repeat_G.gene_name.values)==False)[0]

repeat = pd.concat([repeat_G, repeat_NL.iloc[idx,:]])
repeat = repeat.sort_values(by = ['gene_ind'])
repeat.to_csv("~/deepRIV/UKB/code/hdl/combine_p/repeat_genes.csv")


# analyze repeat
path = "/home/panwei/he000176/deepRIV/UKB/results/hdl/combine_p/repeat_ind/"

hdl = pd.read_csv("/home/panwei/he000176/deepRIV/UKB/results/hdl/test40p_unrelated2/hdl_combined_results.csv")

num_repeats = 21
num_genes = 26
cutoff = 0.05/4701

cauchy_results = pd.DataFrame(np.zeros((num_genes,6)))
cauchy_results.columns = ['gene_names', 'IVa_pval_global', 'IVa_pval_nonlinear', 'IVa_dcor_pval',
		'mu_max', 'mu_min']
hommel_results = pd.DataFrame(np.zeros((num_genes,4)))
hommel_results.columns = ['gene_names', 'global_results', 'nonlinear_results', 'dcor_results']


for i in range(num_genes):
	IVa = pd.read_csv(path + "IVa_{}.txt".format(i), sep = " ")
	IVa = IVa.loc[(IVa.IVa_pval_global != 100.0),:]
	IVa = IVa.loc[(IVa.IVa_pval_nonlinear != 100.0),:]
	cauchy_G_T = np.tan((0.5-IVa.IVa_pval_global)*np.pi).sum() / IVa.shape[0]
	cauchy_NL_T = np.tan((0.5-IVa.IVa_pval_nonlinear)*np.pi).sum() / IVa.shape[0]
	cauchy_dcor_T = np.tan((0.5-IVa.IVa_dcor_pval)*np.pi).sum() / IVa.shape[0]
	cauchy_G_P = 1/2 - np.arctan(cauchy_G_T)/np.pi
	cauchy_NL_P = 1/2 - np.arctan(cauchy_NL_T)/np.pi
	cauchy_dcor_P = 1/2 - np.arctan(cauchy_dcor_T)/np.pi
	cauchy_results.loc[cauchy_results.index[i], 'gene_names'] = IVa.gene_names.iloc[0]
	cauchy_results.loc[cauchy_results.index[i], 'IVa_pval_global'] = cauchy_G_P
	cauchy_results.loc[cauchy_results.index[i], 'IVa_pval_nonlinear'] = cauchy_NL_P
	cauchy_results.loc[cauchy_results.index[i], 'IVa_dcor_pval'] = cauchy_dcor_P
	cauchy_results.loc[cauchy_results.index[i], 'mu_max'] = IVa.mu_max.iloc[0]
	cauchy_results.loc[cauchy_results.index[i], 'mu_min'] = IVa.mu_min.iloc[0]
	#
	Cn = np.sum([1/j for j in range(1,IVa.shape[0]+1)])
	hommel_cutoffs = np.array([k*cutoff / (Cn*IVa.shape[0]) for k in range(1, IVa.shape[0]+1)])
	#
	G_sorted = np.sort(IVa.IVa_pval_global.to_numpy())
	NL_sorted = np.sort(IVa.IVa_pval_nonlinear.to_numpy())
	dcor_sorted = np.sort(IVa.IVa_dcor_pval.to_numpy())
	#
	hommel_results.loc[hommel_results.index[i], 'gene_names'] = IVa.gene_names.iloc[0]
	hommel_results.loc[hommel_results.index[i], 'global_results'] = int(any(G_sorted <= hommel_cutoffs))
	hommel_results.loc[hommel_results.index[i], 'nonlinear_results'] = int(any(NL_sorted <= hommel_cutoffs))
	hommel_results.loc[hommel_results.index[i], 'dcor_results'] = int(any(dcor_sorted <= hommel_cutoffs))


tmp = hdl.reset_index().set_index('gene_names').loc[cauchy_results.gene_names, ['gene_name', 'gene_ind']]
cauchy_results['gene_name'] = tmp.gene_name.values
cauchy_results['gene_ind'] = tmp.gene_ind.values
hommel_results['gene_name'] = tmp.gene_name.values
hommel_results['gene_ind'] = tmp.gene_ind.values



