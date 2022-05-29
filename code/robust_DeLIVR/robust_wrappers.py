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
from tensorflow.keras.layers import LeakyReLU, ReLU, Reshape, Concatenate


def slice_identity(n, column_idx):
	# let I be an n x n identity matrix
	# this function returns I[column_idx,:], the name column_idx is based on how tf.SparseTensor works
	indices_x = tf.reshape(tf.range(tf.shape(column_idx)[0], dtype = tf.int64), (-1,1))
	indices_y = tf.reshape(column_idx, (-1,1))
	st_a = tf.SparseTensor(indices=tf.concat([indices_x, indices_y], axis = -1),
					   values=tf.ones(tf.shape(column_idx)[0]), 
					   dense_shape=[tf.shape(column_idx)[0], n])
	return tf.cast(st_a, tf.float64)

class robust_DeLIVR(tf.keras.Model):
	def __init__(self, mean_x_train, z_train, mean_x_val, z_val, **kwargs):
		super(robust_DeLIVR, self).__init__(**kwargs)
		self.mean_x_train = tf.constant(mean_x_train)
		self.mean_x_val = tf.constant(mean_x_val)
		self.z_train = tf.constant(z_train)
		self.z_val = tf.constant(z_val)
	#
	def train_step(self, data):
		# unpack data
		x, y = data
		W = tf.cast(x[0], dtype = tf.float64) # W = Z(Z'Z)^-1
		I_idx = x[1]
		Identity = slice_identity(tf.shape(self.z_train)[0], I_idx)
		with tf.GradientTape() as tape:
			out = self(self.mean_x_train, training=True)  # Forward pass
			# compute y_pred = (I - W @ z_train.T) @ out)
			# I @ out
			Iout = tf.sparse.sparse_dense_matmul(Identity, out)
			# W @ z_train.T @ out
			tmp = tf.linalg.matmul(W, tf.linalg.matmul(tf.transpose(self.z_train), out))
			# y_pred = Iout - tmp
			y_pred = tf.math.subtract(Iout, tmp)

			# Compute the loss value.
			# The loss function is configured in `compile()`.
			loss = self.compiled_loss(
				y,
				y_pred,
				regularization_losses=self.losses,
			)

		# Compute gradients
		trainable_vars = self.trainable_variables
		gradients = tape.gradient(loss, trainable_vars)

		# Update weights
		self.optimizer.apply_gradients(zip(gradients, trainable_vars))

		# Update the metrics.
		# Metrics are configured in `compile()`.
		self.compiled_metrics.update_state(y, y_pred)

		# Return a dict mapping metric names to current value.
		# Note that it will include the loss (tracked in self.metrics).
		return {m.name: m.result() for m in self.metrics}

	def test_step(self, data):
		# Unpack the data
		x, y = data
		W = tf.cast(x[0], dtype = tf.float64)
		I_idx = x[1]
		Identity = slice_identity(tf.shape(self.z_val)[0], I_idx)

		# Compute predictions
		out = self(self.mean_x_val, training=False)  # Forward pass
		# compute y_pred = (I - W @ z_val.T) @ out)
		# I @ out
		Iout = tf.sparse.sparse_dense_matmul(Identity, out)
		# W @ z_val.T @ out
		tmp = tf.linalg.matmul(W, tf.linalg.matmul(tf.transpose(self.z_val), out))
		# y_pred = Iout - tmp
		y_pred = tf.math.subtract(Iout, tmp)

		# Updates the metrics tracking the loss
		self.compiled_loss(y, y_pred, regularization_losses=self.losses)
		# Update the metrics.
		self.compiled_metrics.update_state(y, y_pred)
		# Return a dict mapping metric names to current value.
		# Note that it will include the loss (tracked in self.metrics).
		return {m.name: m.result() for m in self.metrics}



class data_generator(tf.keras.utils.Sequence):
	def __init__(self, W, I_idx, outputs, batch_size, seed=None):
		self.rng = np.random.RandomState(seed)
		self.W = W.copy()
		self.outputs = outputs.copy()
		self.I_idx = I_idx
		self.batch_size = batch_size
		self.shuffle()

	def __len__(self):
		if isinstance(self.outputs, list):
			return self.outputs[0].shape[0] // self.batch_size
		else:
			return self.outputs.shape[0] // self.batch_size

	def shuffle(self):
		idx = self.rng.permutation(np.arange(self.outputs.shape[0]))
		self.I_idx = self.I_idx[idx]
		self.outputs = self.outputs[idx]
		self.W = self.W[idx,:]
	
	def __getitem__(self, idx):
		I_idx = self.I_idx[idx*self.batch_size:(idx+1)*self.batch_size]
		W = self.W[idx*self.batch_size:(idx+1)*self.batch_size,:]
		batch_x = [W, I_idx]
		batch_y = self.outputs[idx*self.batch_size:(idx+1)*self.batch_size]
		if idx == (len(self) - 1):
			self.shuffle()
		return batch_x, batch_y


def create_stage2(input_shape, 
	ouput_shape,
	mean_x_train, 
	z_train, 
	mean_x_val, 
	z_val,
	layer_dims = [32,16,8,8,8],
	activation = tf.nn.relu,
	l2 = 0.05,
	use_bias = True):
	#
	input1 = Input(shape=input_shape)
	# input_idx = Input(shape=input_shape)
	out = Dense(layer_dims[0], activation=activation, kernel_regularizer=K.regularizers.l2(l2=l2), 
				dtype = tf.float64, use_bias=use_bias)(input1)
	for i in range(1,len(layer_dims)):
		out = Dense(layer_dims[i], activation=activation, kernel_regularizer=K.regularizers.l2(l2=l2), 
					dtype = tf.float64, use_bias=use_bias)(out)
	#
	out = Dense(ouput_shape, activation=None, dtype = tf.float64, use_bias=True)(out)
	skp = Dense(ouput_shape, activation=None, dtype = tf.float64, use_bias=True)(input1)
	out = Add(dtype = tf.float64)([out,skp])		
	rsp = robust_DeLIVR(mean_x_train, 
					z_train, 
					mean_x_val, 
					z_val, inputs=input1, outputs=out)
	return rsp


def fit_stage2(model, X, y, W, X_val, y_val, W_val,
	optimizer = "Adam", 
	loss = "MSE",
	learning_rate = 0.00033, 
	batch_size = 64, 
	training_steps = 8000,
	decay_steps = 200,
	decay_rate = 0.03,
	patience = 20,
	seed = 0):
	# W = Z(Z'Z)^-1
	lr_schedule = tf.keras.optimizers.schedules.ExponentialDecay(
					learning_rate,
					decay_steps=decay_steps,
					decay_rate=decay_rate,
					staircase=True)
	if optimizer=="Adam":
		opt = tf.keras.optimizers.Adam(learning_rate = lr_schedule)
	elif optimizer=="SGD":
		opt = tf.keras.optimizers.SGD(learning_rate = lr_schedule, momentum=0.9)
	#
	if loss == "MSE":
		loss_f = tf.keras.losses.MSE
	model.compile(optimizer=opt,
					loss=loss_f,
					metrics=['MeanSquaredError'])
	erlystp = tf.keras.callbacks.EarlyStopping(monitor='val_mean_squared_error', min_delta = 0.00001, patience = patience)

	generator = data_generator(W, np.arange(y.shape[0]), y, batch_size, seed = seed)
	history = model.fit(generator, batch_size=batch_size, epochs=training_steps, callbacks = [erlystp],
						validation_data=([W_val, np.arange(y_val.shape[0])], y_val), verbose=1)
	return history


