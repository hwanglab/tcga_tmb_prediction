import os
import tensorflow as tf
import numpy as np
import time
from models.layers import relu, matmul, log, sigmoid, softmax
from models.layers import sigmoid_cent, softmax_cent, dropout
from utils.utils import auroc

class MLP(object):
    def __init__(self, dataset, dims, logger=None, measure='auc'):
        self.dataset = dataset
        self.measure = measure
        self.logger  = logger
        self.dims    = dims

    def build(self, reuse=False):
        dataset = self.dataset
        logger = self.logger
        dims = self.dims

        self.xs = tf.placeholder(tf.float32, [None, 2048])
        self.ys = tf.placeholder(tf.float32, [None, dataset.n_classes])

        logger.info('CREATE BASE PARAMETERS')
        with tf.variable_scope('base', reuse=reuse,
                initializer=tf.variance_scaling_initializer()):
            for i in range(3):
                tf.get_variable('layer%d/kernel'%i,
                                    [dims[i], dims[i+1]], trainable=True)
                tf.get_variable('layer%d/biases'%i,
                                    [dims[i+1]], trainable=True)

            tf.get_variable('softmax/kernel',
                    [dims[-1], dataset.n_classes], trainable=True)
            tf.get_variable('softmax/biases',
                    [dataset.n_classes], trainable=True)

    def classify(self):
        xs = self.xs
        ys = self.ys
        with tf.variable_scope('base', reuse=True):
            for i in range(3):
                w = tf.get_variable('layer%d/kernel'%i)
                b = tf.get_variable('layer%d/biases'%i)
                xs = relu(matmul(xs, w) + b)
                xs = dropout(xs, 0.5)

            w = tf.get_variable('softmax/kernel')
            b = tf.get_variable('softmax/biases')
            logits = matmul(xs, w) + b

        if ys.shape[1] == 1:
            cent = sigmoid_cent(logits, ys)
            yhats = sigmoid(logits)
        else:
            cent = softmax_cent(logits, ys)
            yhats = softmax(logits)

        self.logits = logits
        self.yhats  = yhats
        self.cent   = cent

    def build_loss(self, wd_rate=1e-4):
        params = tf.trainable_variables()
        weight_decay = tf.add_n([tf.nn.l2_loss(var) for var in params])
        loss = self.cent + wd_rate*weight_decay

        self.loss = loss
