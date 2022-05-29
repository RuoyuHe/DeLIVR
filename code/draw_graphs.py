import numpy as np
import pandas as pd
import statsmodels.api as sm
from matplotlib import pyplot as plt


sig_idx = idx[0]

new_x = np.linspace(hdl['mu_min'].iloc[sig_idx], hdl['mu_max'].iloc[sig_idx], 500)
pred = rsp1.predict(new_x).squeeze()

plt.xlim((hdl['mu_min'].iloc[0]-0.1, hdl['mu_max'].iloc[0]+0.1))
plt.ylim((1.42,1.5))
#plt.scatter(mean_x_test, tmp, s=2, color = 'orange', label = r'$\hat{E}(g(X)|Z) vs Y$')
plt.plot(new_x, pred, color = "r", label = hdl['gene_name'].iloc[sig_idx])
plt.scatter(mean_x_test, y_test, s=1, label = r"$\hat{\mu}$ vs Y")
plt.legend()
plt.savefig('/home/panwei/he000176/deepRIV/UKB/code/' + hdl['gene_name'].iloc[sig_idx] + '.png')
plt.clf()


#import seaborn as sns
import seaborn.apionly as sns
import matplotlib.ticker as ticker
tmp_idx = np.where(mean_x_test.reshape(-1) >= hdl['mu_min'].iloc[sig_idx])[0]
tmp_idx = np.intersect1d(tmp_idx, np.where(mean_x_test.reshape(-1) <= hdl['mu_max'].iloc[sig_idx])[0])
ax = sns.boxplot(x = mean_x_test.reshape(-1)[tmp_idx], y = y_test[tmp_idx])
ax.set_title(hdl['gene_name'].iloc[sig_idx])
ax.minorticks_off()
ax.set_xticks(np.linspace(hdl['mu_min'].iloc[sig_idx], hdl['mu_max'].iloc[sig_idx], 5))
# ax.xaxis.set_major_locator(ticker.FixedLocator(np.linspace(hdl['mu_min'].iloc[sig_idx], hdl['mu_max'].iloc[sig_idx], 5)))
# ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.02f'))
plt.savefig('/home/panwei/he000176/deepRIV/UKB/code/' + hdl['gene_name'].iloc[sig_idx] + '_boxplot.png')
plt.clf()

unique_x = np.unique(mean_x_test.reshape(-1)[tmp_idx])
for i in range(len(unique_x)):
	v_idx = np.where(mean_x_test.reshape(-1)[tmp_idx] == unique_x[i])[0]
	plt.boxplot(y_test[tmp_idx][v_idx], positions = [unique_x[i]])


related = pd.read_csv("ukb35107_rel_s488265.dat", sep = ' ')


#
# plt.xlim((hdl['mu_min'].iloc[0]-0.1, hdl['mu_max'].iloc[0]+0.1))
# plt.ylim((-0.1,0.2))
# #plt.scatter(mean_x_test, tmp, s=2, color = 'orange', label = r'$\hat{E}(g(X)|Z) vs Y$')
# for i in range(29):
# 	plt.plot(new_x, pred[i,:], color = "r", linestyle = "dashed", alpha=.3)

# plt.plot(new_x, pred[29,:], color = "r", linestyle = "dashed", alpha=.3, label = "RP11-96D1.8")
# plt.plot(unq[:-1], mean_y[:-1], color = "black", linestyle = "dashdot", label = r"$\bar{Y}_{\hat{X_i}}$")
# plt.scatter(mean_x_test, y_test, s=1, color = "lightgrey", label = r"$\hat{\mu}$ vs Y")
# plt.legend()
# plt.savefig('/home/panwei/he000176/deepRIV/UKB/code/RP11-96D1.8_2.png')
# plt.clf()
#


# sig genes
path = "/home/panwei/he000176/deepRIV/UKB/results/hdl/test40p_unrelated2/"
hdl = pd.read_csv(path + "hdl_combined_results.csv")
sig_idx = np.where(hdl['IVa_pval_nonlinear'] <= 0.05/hdl.shape[0])[0]
indices = hdl['gene_ind'].iloc[sig_idx].to_numpy()
gene_names = hdl['gene_name'].iloc[sig_idx].to_numpy()

sig_path = "/home/panwei/he000176/deepRIV/UKB/results/hdl/sig_genes/v3/"
sig_p = pd.read_csv(sig_path + "sig_genes_{}.txt".format(indices[0]), sep = ' ')['IVa_pval_nonlinear']
for v in indices:
	if v == indices[0]:
		continue
	tmp = pd.read_csv(sig_path + "sig_genes_{}.txt".format(v), sep = ' ')['IVa_pval_nonlinear']
	sig_p = pd.concat([sig_p, tmp], axis = 1)

sig_p.columns = gene_names
sig_p.to_csv(sig_path + "sig_genes_combined.csv", index = False)

sig_p = pd.read_csv("sig_genes_combined.csv")
sig_p = -np.log10(sig_p)
sig_p = pd.melt(sig_p, value_vars=sig_p.columns, value_name='-log(p)')
sig_p.columns = ['gene', r'$-log_{10}(p)$']
sig_p.boxplot(by = 'gene', column = r'$-log_{10}(p)$', 
			positions = [1,4,7,10,13,16,19],
			grid = False)
plt.axhline(y = -np.log10(0.05/4701), 
			linestyle = "dashed", 
			color = "r",
			linewidth=1)

plt.show()

