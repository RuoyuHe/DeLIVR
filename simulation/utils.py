import sys
import pickle
import pdb
import numpy as np
import pandas as pd
import tensorflow as tf
import tensorflow.keras as K
import statsmodels.api as sm
from tensorflow.keras import Model
from tensorflow.keras.layers import Activation, Add, Input, Dense, Flatten
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.optimizers import Adam

def create_stage2(input_shape, 
	ouput_shape,
	layer_dims = [32,16,8,8,8],
	activation = tf.nn.relu,
	l2 = 0.05,
	use_bias = True):
	#
	input1 = Input(shape=input_shape)
	out = Dense(layer_dims[0], activation=activation, kernel_regularizer=K.regularizers.l2(l2=l2), 
				dtype = tf.float64, use_bias=use_bias)(input1)
	for i in range(1,len(layer_dims)):
		out = Dense(layer_dims[i], activation=activation, kernel_regularizer=K.regularizers.l2(l2=l2), 
					dtype = tf.float64, use_bias=use_bias)(out)
	#
	out = Dense(ouput_shape, activation=None, dtype = tf.float64, use_bias=True)(out)
	skp = Dense(ouput_shape, activation=None, dtype = tf.float64, use_bias=True)(input1)
	out = Add()([out,skp])
	rsp = Model(inputs=input1, outputs=out)
	return rsp


def fit_stage2(model, X, y, X_val, y_val,
	optimizer = "Adam", 
	loss = "MSE",
	learning_rate = 0.00033, 
	batch_size = 64, 
	epochs = 8000,
	decay_steps = 300,
	decay_rate = 0.95,
	patience = 20,
	shuffle = True):
	#
	lr_schedule = tf.keras.optimizers.schedules.ExponentialDecay(
					learning_rate,
					decay_steps=decay_steps,
					decay_rate=decay_rate,
					staircase=True)
	if optimizer=="Adam":
		opt = tf.keras.optimizers.Adam(learning_rate = lr_schedule)
	elif optimizer=="SGD":
		opt = tf.keras.optimizers.SGD(learning_rate = lr_schedule, momentum=0.9)

	if loss == "MSE":
		loss_f = tf.keras.losses.MSE
	model.compile(optimizer=opt,
					loss=loss_f,
					metrics=['MeanSquaredError'])
	erlystp = tf.keras.callbacks.EarlyStopping(monitor='val_mean_squared_error', 
		min_delta = 0.00001, patience = patience, restore_best_weights = True)
	history = model.fit(X, y, batch_size=batch_size, epochs=epochs, callbacks = [erlystp],
						validation_data=(X_val, y_val), shuffle = shuffle, verbose=1)
	return history


def stage2Tests(x_test, y_test, pred):
	# global tests
	design = sm.add_constant(pred.reshape(-1,1))
	testM = sm.OLS(y_test, design)
	testM_results = testM.fit()
	p_global = testM_results.pvalues[-1]
	t_global = testM_results.tvalues[-1]
	linear_coef = testM_results.params[-1]
	linear_se = testM_results.bse[-1]
	#
	# Nonlinearity test
	testM = sm.OLS(y_test, sm.add_constant(x_test))
	testM_results = testM.fit()
	residuals = y_test - testM_results.predict(sm.add_constant(x_test))
	testM = sm.OLS(residuals, sm.add_constant(pred.reshape(-1,1)))
	testM_results = testM.fit()
	p_nonlinear = testM_results.pvalues[-1]
	t_nonlinear = testM_results.tvalues[-1]
	nonlinear_coef = testM_results.params[-1]
	nonlinear_se = testM_results.bse[-1]
	#
	return p_global, t_global, p_nonlinear, t_nonlinear, nonlinear_coef, nonlinear_se, linear_coef, linear_se


def cauchy_p(p_val, num_repeats):
	'''
		p_val: matrix of p-values, shape = (# samples, # repeats)
	'''
	cauchy_stat = np.tan((0.5-p_val)*np.pi).sum(axis = 1) / num_repeats
	P = 1/2 - np.arctan(cauchy_stat)/np.pi
	return P



def main():
	pass

if __name__ == '__main__':
	import sys
	sys.exit(int(main() or 0))