import sys
import numpy as np
import statsmodels.api as sm
from sklearn.model_selection import train_test_split
from utils import create_stage2, fit_stage2, stage2Tests, cauchy_p

class DeLIVR():
	def __init__(self, s1_IV, s1_expr, s2_IV, s2_pheno):
		self._s1_IV = s1_IV
		self._s1_expr = s1_expr
		self._IV = s2_IV
		self._pheno = s2_pheno
		self._expr = None
	#
	def sample_split(self, test_ratio, val_ratio, seed = None):
		if seed is None:
			seed = np.random.randint(0,2**32 - 1)
		if self._expr is None:
			self.fit_stage1()
		#
		self._expr_train, self._expr_val, self._pheno_train, self._pheno_val \
			= train_test_split(self._expr, self._pheno, test_size=val_ratio, random_state=seed)
		#
		self._expr_train, self._expr_test, self._pheno_train, self._pheno_test \
			= train_test_split(self._expr_train, self._pheno_train, 
				test_size=test_ratio/(1-val_ratio), random_state=seed)
	#
	def fit_stage1(self):
		design = sm.add_constant(self._s1_IV)
		stage1 = sm.OLS(self._s1_expr, design)
		stage1 = stage1.fit()
		s2_design = sm.add_constant(self._IV)
		self._expr = stage1.predict(s2_design).reshape(-1,1)
	#
	def fit(self, split_ratios, training_params, n_repeats):
		n_ratios = len(split_ratios)
		n_params = len(training_params)
		total_repeats = n_ratios*n_params*n_repeats
		p_global = np.ones(total_repeats)
		p_nl = np.ones(total_repeats)
		#
		counter = 0
		for i in range(n_ratios):
			for j in range(n_params):
				for k in range(n_repeats):
					self.sample_split(**split_ratios[i])
					s2_model = create_stage2((1,),1, l2 = 0)
					h = fit_stage2(s2_model, self._expr_train, self._pheno_train, 
						self._expr_val, self._pheno_val, **training_params[j])
					tmp = s2_model.predict(self._expr_test).squeeze()
					p_global[counter], _, p_nl[counter], _, _, _, _, _ \
						= stage2Tests(self._expr_test, self._pheno_test, tmp)
					counter += 1
		p_global = cauchy_p(p_global.reshape(-1,total_repeats), total_repeats)
		p_nl = cauchy_p(p_nl.reshape(-1,total_repeats), total_repeats)

		return p_global, p_nl


def generate_data(IV, beta, true_g, gamma, seed = 0):
	if seed is None:
		seed = np.random.randint(0,2**32 - 1)
	np.random.seed(seed)
	s1_size = IV.shape[0]//10
	s1_idx = np.random.choice(np.arange(IV.shape[0]), size = s1_size, replace = False)
	s2_idx = np.array(list(set(np.arange(IV.shape[0])).difference(s1_idx)))
	s1_IV = IV[s1_idx,:]
	s2_IV = IV[s2_idx,:]
	# stage 1 data
	s1_expr = (s1_IV@beta).squeeze() + np.random.normal(size = s1_size)
	s1_expr = s1_expr.reshape(-1,1)
	#
	e = np.random.multivariate_normal(mean = [0,0], cov = [[1,gamma],[gamma,1]], size = s2_IV.shape[0])
	s2_expr = (s2_IV@beta).squeeze() + e[:,0]
	s2_pheno = true_g(s2_expr) + e[:,1]
	s2_expr = s2_expr.reshape(-1,1)

	return s1_IV, s1_expr, s2_IV, s2_pheno



def main():
	pass

if __name__ == '__main__':
	import sys
	sys.exit(int(main() or 0))

