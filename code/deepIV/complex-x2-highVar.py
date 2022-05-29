import sys
print(sys.version)
import pickle
import pdb
import numpy as np
import pandas as pd
import tensorflow as tf
import tensorflow.keras as K
import statsmodels.api as sm
from tensorflow.keras import Model
from tensorflow.keras.layers import Activation, Add, Input, Dense, Dropout, Flatten, Conv2D, MaxPooling2D
from tensorflow.keras.layers import BatchNormalization as BN
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.optimizers import Adam
from sklearn.model_selection import StratifiedShuffleSplit
from tensorflow.keras.layers import LeakyReLU, ReLU, Reshape
from tensorflow.keras.backend import expand_dims, squeeze, clear_session
from tensorflow.keras.preprocessing.sequence import pad_sequences
from tensorflow.keras.layers import ZeroPadding2D, Concatenate, GlobalAveragePooling1D, GlobalAveragePooling2D
from tensorflow.keras.constraints import max_norm
from sklearn.model_selection import train_test_split
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from neuralnet import gaussianLoss, UniNormalSampling, stage2Loss, MCsample, ResNet
from fit_model import deeprivOPT, fit_deepRIV, unbiased_grad, fit_deepIV
from wrappers import create_stage1, create_stage2, fit_stage2, tune_l2, stage2Tests



z = pd.read_csv("../data/simulated_SNP_July18.txt", sep = " ")
z = z.to_numpy().astype(np.float64)
beta = pd.read_csv("../data/simulated_beta_July18.txt", sep = " ")
beta = beta.to_numpy().squeeze()
beta_idx = np.where(beta!=0)[0]
z = z[:,beta_idx]
beta = beta[beta_idx]
true_x = pd.read_csv("../data/simulated_X_July18.txt", sep = " ")
true_x = true_x.to_numpy().squeeze()

save_path = "/home/panwei/he000176/deepRIV/deepIV/results/x2-highVar/"

# gamma is the correlation btw x and y
true_g = lambda x: 0.5*x**2
gamma = 0
sig = 1
max_true_x = np.max(np.abs(true_x))
new_x = np.linspace(-max_true_x, max_true_x,500)
new_y = true_g(new_x)


runs = 20
pred_total = np.zeros((runs,500))
prednn_total = np.zeros((runs,500))
for rep in range(runs):
	e = np.random.multivariate_normal(mean = [0,0], cov = [[sig,gamma],[gamma,1]], size = true_x.shape[0])
	U = np.random.normal(0,6,true_x.shape[0])
	#U = 0
	x = true_x + U + e[:,0]
	y = true_g(x) - 2*U**2 + e[:,1]
	x = x.reshape(-1,1)

	# traing, val, test split, .6/.2/.2
	z_train, z_test, x_train, x_test, true_x_train, true_x_test, y_train, y_test = train_test_split(z, x, true_x.reshape(-1,1), y, test_size=0.2, random_state=1)
	z_train, z_val, x_train, x_val, true_x_train, true_x_val, y_train, y_val = train_test_split(z_train, x_train, true_x_train, y_train, test_size=0.25, random_state=1) 

	### stage1
	design = sm.add_constant(z_train)
	trt1M = sm.OLS(x_train, design)
	stage1 = trt1M.fit()

	### complex version
	n_draws = 100
	learning_rate = 0.00033
	batch_size = 640
	training_steps = 32000
	display_step = 1
	resample_prop = 1

	mean_x_train = stage1.predict(sm.add_constant(z_train)).reshape(-1,1)
	mean_x_val = stage1.predict(sm.add_constant(z_val)).reshape(-1,1)
	rss = ((x_train - mean_x_train)**2).squeeze()
	rss = rss[np.where(rss<30)[0]]
	std_x_train = np.sqrt(np.sum(rss)/rss.shape[0])

	pred_x_train = np.random.normal(mean_x_train, std_x_train, (mean_x_train.shape[0],n_draws))
	pred_x_train2 = np.random.normal(mean_x_train, std_x_train, (mean_x_train.shape[0],n_draws))
	pred_x_val = np.random.normal(mean_x_val, std_x_train, (mean_x_val.shape[0],n_draws))

	stage2 = create_stage2((1,),1, l2 = 2)
	lr_schedule = tf.keras.optimizers.schedules.ExponentialDecay(
		learning_rate,
		decay_steps=8000,
		decay_rate=0.7,
		staircase=True)
	optimizer = tf.keras.optimizers.Adam(lr_schedule)
	fit_deepIV(stage2, optimizer, unbiased_grad, stage2Loss, mean_x_train, std_x_train, 
						pred_x_train, pred_x_train2, pred_x_val, z_train, y_train, y_val,
						n_draws, batch_size, training_steps,
						display_step = 1, resample_prop = 1, early_stopping = 500)
	pred_total[rep,:] = stage2(new_x.reshape(-1,1), training = False).numpy().squeeze()

treat_pred = stage1.predict(sm.add_constant(z_test))
plt.scatter(true_x_test, treat_pred, s=1)
plt.plot(new_x,new_x, color = "black")
plt.xlabel("true mean")
plt.ylabel("predicted mean")
plt.savefig(save_path + "treat_X10.png")
plt.clf()

for i in range(1,runs):
	plt.plot(new_x, pred_total[i,:], color = "r", alpha = 0.5)

plt.plot(new_x, pred_total[0,:], color = "r", alpha = 0.5,label = "deepIV")
plt.plot(new_x, new_y, color = "black", label = "True")
# plt.plot(new_x, prednn_total[0,:], color = "green", alpha = 1, label = "Direct NN")
# plt.scatter(x2,np.repeat(y, n_draws), color = "orange", s=1, label = "sampled X")
plt.legend()
plt.savefig(save_path + "compare_X10.png")
plt.clf()

np.savetxt(save_path + "pred_X10.txt", pred_total)
np.savetxt(save_path + "treat_pred_X10.txt", treat_pred)

# np.savez("../results/model/data_used.npz", x = x, y = y, z_train = z_train, z_val = z_val, z_test = z_test,
# 			x_train = x_train, x_val = x_val, x_test = x_test, true_x_train = true_x_train, true_x_val = true_x_val,
# 			true_x_test = true_x_test, y_train = y_train, y_val = y_val, y_test = y_test )
# stage1.save("../results/model/stage1")
# stage2.save("../results/model/stage2_888")
# dnn.save("../results/model/dnn")


# for i in range(1,runs):
# 	plt.plot(new_x, prednn_total[i,:], color = "green", alpha = 0.5)

# plt.plot(new_x, prednn_total[0,:], color = "green", alpha = 0.5, label = "Direct NN")
# plt.plot(new_x, new_y, color = "black", label = "True")
# plt.legend()
