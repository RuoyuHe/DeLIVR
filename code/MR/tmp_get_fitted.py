# included independent tests for Y_hat and residuals
import os
import sys
print(sys.version)
import pdb
import dcor
import pickle
import itertools
import numpy as np
import pandas as pd
import statsmodels.api as sm
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt


# Defining the R script and loading the instance in Python
r = robjects.r
r['source']('/home/panwei/he000176/deepRIV/UKB/code/MR/hdl_MR_stage1_FUN.R')

# Loading the function we have defined in R.
stage1_r = robjects.globalenv['stage1']

exposure = 'bmi'
out_path = '/home/panwei/he000176/deepRIV/UKB/results/MR/hdl/' + exposure + '_stage2_data.csv'

#converting it into r object for passing into r function
with localconverter(robjects.default_converter + pandas2ri.converter):
	stage1_results = stage1_r('bmi')
	#
	mean_X = robjects.conversion.rpy2py(stage1_results[0])
	Y = robjects.conversion.rpy2py(stage1_results[1])
	LR_pval = robjects.conversion.rpy2py(stage1_results[2])
	LR_tval = robjects.conversion.rpy2py(stage1_results[3])
	polyq_pval = robjects.conversion.rpy2py(stage1_results[4])
	polyq_tval = robjects.conversion.rpy2py(stage1_results[5])
	lq_global = robjects.conversion.rpy2py(stage1_results[6])
	lq_nonlinear = robjects.conversion.rpy2py(stage1_results[7])
	num_snps = robjects.conversion.rpy2py(stage1_results[8])
	num_samples = robjects.conversion.rpy2py(stage1_results[9])

print(Y.shape)
print(mean_X.shape)
tmp = pd.DataFrame({'y':Y.squeeze(), 'mean_X':mean_X.squeeze()})
tmp.to_csv(out_path, index = False)


# Defining the R script and loading the instance in Python
r = robjects.r
r['source']('/home/panwei/he000176/deepRIV/UKB/code/MR/ldl_MR_stage1_FUN.R')

# Loading the function we have defined in R.
stage1_r = robjects.globalenv['stage1']


out_path = '/home/panwei/he000176/deepRIV/UKB/results/MR/ldl/' + exposure + '_stage2_data.csv'

#converting it into r object for passing into r function
with localconverter(robjects.default_converter + pandas2ri.converter):
	stage1_results = stage1_r('bmi')
	#
	mean_X = robjects.conversion.rpy2py(stage1_results[0])
	Y = robjects.conversion.rpy2py(stage1_results[1])
	LR_pval = robjects.conversion.rpy2py(stage1_results[2])
	LR_tval = robjects.conversion.rpy2py(stage1_results[3])
	polyq_pval = robjects.conversion.rpy2py(stage1_results[4])
	polyq_tval = robjects.conversion.rpy2py(stage1_results[5])
	lq_global = robjects.conversion.rpy2py(stage1_results[6])
	lq_nonlinear = robjects.conversion.rpy2py(stage1_results[7])
	num_snps = robjects.conversion.rpy2py(stage1_results[8])
	num_samples = robjects.conversion.rpy2py(stage1_results[9])

print(Y.shape)
print(mean_X.shape)
tmp = pd.DataFrame({'y':Y.squeeze(), 'mean_X':mean_X.squeeze()})
tmp.to_csv(out_path, index = False)
