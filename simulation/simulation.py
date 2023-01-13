import numpy as np
import pandas as pd
from DeLIVR import DeLIVR, generate_data

file_num = int(sys.argv[-1])

split_ratios = ({'test_ratio':0.4, 
				'val_ratio':0.1},)

training_params = ({'learning_rate':0.00001,
					'epochs':8000,
					'decay_steps':300,
					'decay_rate':0.95,
					'patience':30},)

n_repeats = 2
runs = 10

z = pd.read_csv('simulated_IV.txt', sep = ' ').to_numpy()
beta = pd.read_csv('s1_beta.txt', sep = ' ').to_numpy()
nonzero_idx = np.where(beta.squeeze()!=0)[0]
z = z[:, nonzero_idx]
beta = beta[nonzero_idx, :]

true_g = lambda x: 3*x**2
gamma = 0.7

IVa_path = '/path/to/results/'
for run in range(runs):
	s1_IV, s1_expr, s2_IV, s2_pheno = generate_data(z, beta, true_g, gamma)
	m = DeLIVR(s1_IV, s1_expr, s2_IV, s2_pheno)
	p_global, p_nl = m.fit(split_ratios, training_params, n_repeats)
	IVa_results = pd.DataFrame({"gene_names":'CDK2AP1',
								"IVa_pval_global":p_global, 
								"IVa_pval_nonlinear":p_nl,
								"num_snps":z.shape[1]
								}, index = [0])
	IVa_results.to_csv(IVa_path, mode = 'a', index = False, 
							sep = " ", header = not os.path.exists(IVa_path))

