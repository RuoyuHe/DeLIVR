import numpy as np
import pandas as pd
import tensorflow as tf
import tensorflow.keras as K
from neuralnet import deeprivLoss, stage2Loss

def deeprivOPT(model_g,model_d,x,x2,y,z,batch_size,n_draws,opt,lossFUN):
	'''
		model_g: the model estimating g(x)
		model_d: the model estimating the direct effect Z*alpha
	'''
	pred2 = model_g(tf.reshape(x2,(-1,1)))
	pred2 = tf.reduce_mean(tf.reshape(pred2, (-1, n_draws)), axis = 1)
	asdf = -2*(y-pred2)/batch_size
	with tf.GradientTape() as g:
		pred_g = model_g(x.reshape(-1,1))
		pred_g = tf.reduce_mean(tf.reshape(pred_g, (-1, n_draws)), axis = 1)
		pred_d,pen,var_y = model_d(z)
		trainable_variables = model_g.trainable_variables
		tmp = len(trainable_variables)
		trainable_variables.extend(model_d.trainable_variables)
	gradients = g.gradient(pred_g, trainable_variables[:tmp], output_gradients = asdf/(var_y**2))
	with tf.GradientTape() as g:
		pred_d,pen,var_y = model_d(z)
		pred_d = tf.squeeze(pred_d)
		loss = lossFUN(y,pred_g,pred_d,var_y,pen,batch_size)
	gradients.extend(g.gradient(loss, trainable_variables[tmp:]))
	opt.apply_gradients(zip(gradients, trainable_variables))


def fit_deepRIV(model_g, model_d, optimizer, optFUN, lossFUN, mean_x_train, std_x_train, 
				pred_x_train, pred_x_train2, pred_x_val, z_train, z_val, y_train, y_val,
				n_draws, batch_size, training_steps,
				display_step = 1, resample_prop = 1, early_stopping = 500):
	train_data = tf.data.Dataset.from_tensor_slices(np.arange(pred_x_train.shape[0]))
	train_data = train_data.repeat().shuffle(500).batch(batch_size).prefetch(50)
	# train model
	losses = np.zeros(early_stopping)
	best_loss = 1e9
	k=0
	for step, batch_idx in enumerate(train_data.take(training_steps), 1):
		# resample part of x_train
		n_resample = int(batch_size*resample_prop)
		resample_idx = np.random.choice(batch_idx, n_resample, replace = False)
		pred_x_train[resample_idx,:] = np.random.normal(mean_x_train[resample_idx,:], std_x_train, 
														(resample_idx.shape[0],n_draws))
		pred_x_train2[resample_idx,:] = np.random.normal(mean_x_train[resample_idx,:], std_x_train, 
														(resample_idx.shape[0],n_draws))
		batch_x = pred_x_train[batch_idx,:]
		batch_x2 = pred_x_train2[batch_idx,:]
		batch_z = z_train[batch_idx,:]
		batch_y = y_train[batch_idx]
		# Run the optimization to update W and b values.
		optFUN(model_g, model_d, batch_x, batch_x2, batch_y, batch_z, batch_size, n_draws, optimizer,
				lossFUN)
		pred = model_g(pred_x_val.reshape(-1,1), training = False)
		pred = tf.reduce_mean(tf.reshape(pred, (-1, n_draws)) ,axis = 1)
		pred_d,penalty,var_y = model_d(z_val)
		pred_d = tf.squeeze(pred_d)
		loss = lossFUN(y_val, pred, pred_d, var_y, penalty, y_val.shape[0]).numpy()
		if step % display_step == 0:
			print("step: {}, loss: {}, lambda: {}, l2_norm: {}".format(step, loss[0], 1/model_d.sigma[0].numpy()**2, 
					2*tf.nn.l2_loss(model_d.w)))
		if k >= early_stopping:
			mean_loss = np.mean(losses)
			if mean_loss < best_loss:
				best_loss = mean_loss
				losses = np.zeros(early_stopping)
				k=0
			else:
				break
		losses[k] = loss
		k+=1


@tf.function
def unbiased_grad(model, x, x2, y, batch_size, n_draws, optimizer):
	pred2 = model(tf.reshape(x2,(-1,1)))
	pred2 = tf.reduce_mean(tf.reshape(pred2, (-1, n_draws)), axis = 1)
	asdf = -2*(y-pred2)/batch_size
	with tf.GradientTape() as g:
		pred = model(tf.reshape(x,(-1,1)))
		pred = tf.reduce_mean(tf.reshape(pred, (-1, n_draws)), axis = 1)
		trainable_variables = model.trainable_variables
	grad = g.gradient(pred, trainable_variables, output_gradients = asdf)
	optimizer.apply_gradients(zip(grad, trainable_variables))


@tf.function
def evaluate(model, inputs, **kwargs):
	return model(inputs, **kwargs)


def fit_deepIV(model, optimizer, optFUN, lossFUN, mean_x_train, std_x_train, 
				pred_x_train, pred_x_train2, pred_x_val, y_train, y_val,
				n_draws, batch_size, training_steps,
				display_step = 1, resample_prop = 1, early_stopping = 500,
				saving_path = None):
	train_data = tf.data.Dataset.from_tensor_slices(np.arange(pred_x_train.shape[0]))
	train_data = train_data.repeat().shuffle(500).batch(batch_size).prefetch(50)
	# train stage1
	total_losses = np.zeros(training_steps)
	losses = np.zeros(early_stopping)
	best_loss = 1e9
	mean_loss = 1e9
	k=0
	for step, batch_idx in enumerate(train_data.take(training_steps), 1):
		# print("step: {}".format(step))
		# resample part of x_train
		n_resample = int(batch_size*resample_prop)
		resample_idx = np.random.choice(batch_idx, n_resample, replace = False)
		pred_x_train[resample_idx,:] = np.random.normal(mean_x_train[resample_idx,:], std_x_train, 
														(resample_idx.shape[0],n_draws))
		pred_x_train2[resample_idx,:] = np.random.normal(mean_x_train[resample_idx,:], std_x_train, 
														(resample_idx.shape[0],n_draws))
		batch_x = pred_x_train[batch_idx,:]
		batch_x2 = pred_x_train2[batch_idx,:]
		batch_y = y_train[batch_idx]
		# Run the optimization to update W and b values.
		optFUN(model, batch_x, batch_x2, batch_y, batch_size, n_draws, optimizer)
		pred = evaluate(model, inputs = pred_x_val.reshape(-1,1), training = False)
		# pred = model(pred_x_val.reshape(-1,1), training = False)
		pred = tf.reduce_mean(tf.reshape(pred, (-1, n_draws)) ,axis = 1)
		loss = lossFUN(y_val, pred)
		if step % display_step == 0:
			print("step: %i, loss: %f, mean loss: %f, best_loss: %f" % (step, loss, mean_loss, best_loss))
		if k >= early_stopping:
			mean_loss = np.mean(losses)
			if mean_loss < best_loss:
				best_loss = mean_loss
				losses = np.zeros(early_stopping)
				k=0
			else:
				break
		total_losses[step-1] = loss
		losses[k] = loss
		k+=1

	if saving_path:
		pd.DataFrame({"losses":total_losses}).to_csv(saving_path, index = False)

	