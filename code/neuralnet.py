import pickle
import numpy as np
import pandas as pd
import tensorflow as tf
import tensorflow.keras as K
from tensorflow.keras.utils import to_categorical
from tensorflow.keras import Model
from tensorflow.keras.layers import Activation, Add, Input, Dense, Dropout, Flatten, MaxPooling2D
from tensorflow.keras.layers import BatchNormalization as BN
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.optimizers import Adam
from sklearn.model_selection import StratifiedShuffleSplit
from tensorflow.keras.layers import LeakyReLU, ReLU, Reshape
from tensorflow.keras.backend import expand_dims, squeeze, clear_session
from tensorflow.keras.preprocessing.sequence import pad_sequences
from tensorflow.keras.layers import ZeroPadding2D, Concatenate, GlobalAveragePooling1D, GlobalAveragePooling2D
from tensorflow.keras.constraints import max_norm


class ResnetBlock(Model):
	#
	def __init__(self, dim, activation, l1=0., l2=0., residual_path=False):
		super(ResnetBlock, self).__init__()
		#
		self.layer_dim = dim
		self.activation = activation
		self.residual_path = residual_path
		self.l1 = l1
		self.l2 = l2
		#
		self.layer1 = Dense(self.layer_dim, activation=None, 
							kernel_regularizer=K.regularizers.l2(l2=self.l2), 
							use_bias=True, dtype = tf.float64)
		self.bn1 = BN()
		self.layer2 = Dense(self.layer_dim, activation=None, 
							kernel_regularizer=K.regularizers.l2(l2=self.l2),
							use_bias=True, dtype = tf.float64)
		self.bn2 = BN()
		#
		if residual_path:
			self.down_layer = Dense(self.layer_dim, activation=None, 
									kernel_regularizer=K.regularizers.l2(l2=self.l2),
									use_bias=True, dtype = tf.float64)
			self.down_bn = BN()
	#
	#@tf.function
	def call(self, inputs, training=None):
		inputs = tf.cast(inputs, tf.float64)
		residual = inputs
		#
		x1 = self.bn1(inputs, training=training)
		x1 = self.activation(x1)
		x = self.layer1(x1)
		x = self.bn2(x, training=training)
		x = self.activation(x)
		x = self.layer2(x)
		#
		# this module can be added into self.
		# however, module in for can not be added.
		if self.residual_path:
			residual = self.down_layer(x1)
		#
		x = 0.4*x + residual
		return x


class ResNet(Model):
	#
	def __init__(self, block_list, activation, n_draws, l1 = 0., l2 = 0., **kwargs):
		super(ResNet, self).__init__(**kwargs)
		#
		self.num_blocks = len(block_list)
		self.block_list = block_list
		self.activation = activation
		self.n_draws = n_draws
		self.l1 = l1
		self.l2 = l2
		#
		self.layer_initial = Dense(block_list[0], activation=None, 
									kernel_regularizer=K.regularizers.l2(l2=self.l2),
									use_bias=True, dtype = tf.float64)
		#
		self.blocks = tf.keras.models.Sequential(name='dynamic-blocks')
		# build all the blocks
		for i in range(len(block_list)):
			if i>0 and block_list[i]!=block_list[i-1]:
				block = ResnetBlock(block_list[i], self.activation, l1 = self.l1,
									l2 = self.l2, residual_path=True)
			else:
				block = ResnetBlock(block_list[i], self.activation, l1 = self.l1,
									l2 = self.l2)
			self.blocks.add(block)
		#
		self.final_bn = BN()
		self.outlayer = Dense(1, activation=None, use_bias=True, dtype = tf.float64)
	#
	def reparametrize(self, inputs):
		self.mu = np.repeat(inputs[:,0], self.n_draws)
		self.sigma = np.repeat(inputs[:,1], self.n_draws)
		self.new_x = self.mu + self.sigma*np.random.normal(0,1,self.mu.shape[0])
		self.new_x = tf.reshape(self.new_x, (-1,1))
		return self.new_x
	#
	#@tf.function
	def call(self, inputs, training = True):
		# if training:
		# 	inputs = self.reparametrize(inputs)
		out = self.layer_initial(inputs)
		out = self.blocks(out, training=training)
		out = self.final_bn(out, training=training)
		out = self.activation(out)
		out = self.outlayer(out)
		# if training:
		# 	out = tf.reshape(out, (-1,self.n_draws))
		# 	out = tf.math.reduce_mean(out, axis = 1)
		#
		return out

class directM(K.layers.Layer):
	#
	def __init__(self, input_dim):
		super(directM, self).__init__()
		self.input_dim = input_dim
		self.w = self.add_weight('w', [self.input_dim, 1],
						initializer = tf.keras.initializers.HeNormal(), trainable = True, dtype = tf.float64)
		# [dim_out]
		self.b = self.add_weight('intercept', [1],
						initializer = tf.keras.initializers.HeNormal(), trainable = True, dtype = tf.float64)
		# # coef
		# self.coef = self.add_weight('coef', [1],
		# 				initializer = tf.keras.initializers.HeNormal(), trainable = True, dtype = tf.float64)
		# mean
		self.mu = self.add_weight('mu', [1], 
						initializer = tf.keras.initializers.HeNormal(), trainable = True, dtype = tf.float64)
		# l2 penalty
		self.sigma = self.add_weight('sigma', [1], 
						initializer = tf.keras.initializers.HeNormal(), trainable = True, dtype = tf.float64)
		# variance of y
		self.sigma_y = self.add_weight('sigma_y', [1], 
						initializer = tf.keras.initializers.HeNormal(), trainable = True, dtype = tf.float64)
		print(self.w.shape, self.b.shape)
		print(type(self.w), tf.is_tensor(self.w), self.w.name)
		print(type(self.b), tf.is_tensor(self.b), self.b.name)
	#
	def call(self, x):
		x = tf.cast(x, tf.float64)
		x = tf.matmul(x, self.w) + self.b
		penalty = 2*tf.nn.l2_loss(self.w-self.mu)/(self.sigma**2) + tf.math.log(self.sigma**2)
		#
		return x, penalty, self.sigma_y


def gaussianLoss(true, pred, pred_sig):
	pi = tf.constant(np.pi, dtype = tf.float64)
	kernel = tf.math.reduce_mean(tf.math.divide(tf.math.pow(true - pred, 2),
							2*tf.math.pow(pred_sig,2)))
	loss = 0.5*tf.math.log(2*pi) + kernel + tf.math.reduce_mean(tf.math.log(pred_sig + tf.constant(1e-06,
						dtype = tf.float64)))
	return loss

def UniNormalSampling(mu, sd, n):
	# mu and sd should be of shape (# samples, 1)
	return mu + sd*tf.random.normal((mu.shape[0],n))

def stage2Loss(y,x):
	return tf.math.reduce_mean(tf.math.pow(y - x, 2))


def deeprivLoss(y,pred_g,pred_d,var_y,pen,N):
	'''
		pred_g: prediction of E(g(x)|Z)
		pred_d: prediction of Z*alpha, the direct effect
		pen: lambda*l2_norm of weights
	'''
	return tf.math.reduce_sum(tf.math.pow(y - pred_g - pred_d, 2))/(var_y**2) + N*tf.math.log(var_y**2) + pen


def MCsample(model, input_train, input_test, output_train, n_draws):
	pred_x_train = model.predict(input_train)
	pred_x_test = model.predict(input_test)
	rss = ((output_train - pred_x_train)**2).squeeze()
	rss = rss[np.where(rss<30)[0]]
	std_x_train = np.sqrt(np.sum(rss)/rss.shape[0])
	#
	pred_x_train = np.repeat(pred_x_train, n_draws) + std_x_train*np.random.normal(0,1,pred_x_train.shape[0]*n_draws)
	pred_x_train = pred_x_train.reshape(-1,1)
	pred_x_test = np.repeat(pred_x_test, n_draws) + std_x_train*np.random.normal(0,1,pred_x_test.shape[0]*n_draws)
	pred_x_test = pred_x_test.reshape(-1,1)
	return pred_x_train, pred_x_test


class MC2sls(K.layers.Layer):
	def __init__(self, input_dim):
		super(MC2sls, self).__init__()
		self.input_dim = input_dim
		self.w0 = self.add_weight('w0', [self.input_dim, 1],
						initializer = tf.keras.initializers.HeNormal(), trainable = True, dtype = tf.float64)
		self.w = self.add_weight('w', [self.input_dim, 1],
						initializer = tf.keras.initializers.HeNormal(), trainable = True, dtype = tf.float64)
		# [dim_out]
		self.b = self.add_weight('intercept', [1],
						initializer = tf.keras.initializers.HeNormal(), trainable = True, dtype = tf.float64)
	#
	def call(self, x):
		x = tf.cast(x, tf.float64)
		x2 = x**2
		x = tf.matmul(x, self.w0) + tf.matmul(x2, self.w) + self.b
		return x

