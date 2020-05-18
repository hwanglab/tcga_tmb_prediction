import tensorflow as tf
import numpy as np

eps = 1e-10
log = lambda x: tf.log(x+eps)
relu = tf.nn.relu
matmul = tf.matmul
multiply = tf.multiply
sigmoid = tf.nn.sigmoid
softmax = tf.nn.softmax
flatten = tf.layers.flatten
dropout = tf.layers.dropout

conv2d = tf.nn.conv2d
bias_add = tf.nn.bias_add
max_pool = tf.nn.max_pool

def sigmoid_cent(logits, labels, avg=True):
    cent = tf.nn.sigmoid_cross_entropy_with_logits(labels=labels, logits=logits)
    if avg:
        cent = tf.reduce_mean(cent)
    return cent

def softmax_cent(logits, labels, avg=True):
    cent = tf.nn.softmax_cross_entropy_with_logits_v2(labels=labels, logits=logits)
    if avg:
        cent = tf.reduce_mean(cent)
    return cent

def get_staircase_lr(global_step, bdrs, vals):
    lr = tf.train.piecewise_constant(tf.to_int32(global_step), bdrs, vals)
    return lr
