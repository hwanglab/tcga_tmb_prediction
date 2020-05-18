import os
import glob
import numpy as np
import argparse
import logging
import tensorflow as tf

from utils.utils import cv_generate_index, load_tmb_files_v2, centering, auroc
from utils.paths import TMB_DATA_PATH, TMB_INFO_FNAME, CV_PATH_DIR
from models.mlp import MLP
from models.layers import get_staircase_lr
from datasets.dataset import DataSet

np.random.seed(6150)

parser = argparse.ArgumentParser()
parser.add_argument('--results_path', type=str,
                    default='../results/mlp')
parser.add_argument('--n_folds', type=int, default=5)
parser.add_argument('--n_exps', type=int, default=5)
parser.add_argument('--n_epochs', type=int, default=50)
parser.add_argument('--batch_size', type=int, default=100)
parser.add_argument('--init_lr', type=float, default=1e-2)
parser.add_argument('--mode', type=str, default='train')
parser.add_argument('--gpu_id', type=str, default='7')

args = parser.parse_args()

os.environ['CUDA_VISIBLE_DEVICES'] = args.gpu_id

if not os.path.exists(args.results_path):
    os.makedirs(args.results_path)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('mlp')
logger.addHandler(logging.FileHandler(
    os.path.join(args.results_path, '%s.log'%args.mode), mode='w')
)

n_exps = args.n_exps
n_folds = args.n_folds

xs, ys = load_tmb_files_v2(TMB_DATA_PATH, TMB_INFO_FNAME)
xs = [np.mean(xs_, axis=0, keepdims=True) for xs_ in xs]
xs = np.concatenate(xs, axis=0)
n_bags = xs.shape[0]

ys = (ys+1)/2
lr_rate = args.init_lr

def train():
    for rep in range(n_exps):
        idx_list = cv_generate_index(rep, n_folds, ys, CV_PATH_DIR)
        aucs_trn = []
        aucs_tst = []
        for k in range(n_folds):
            logger.info("Ex%d - %dth fold" % (rep, k))

            tst_idx = [int(idx) for idx in idx_list[k]]
            trn_idx = [idx for idx in range(n_bags) if idx not in tst_idx]

            for idx in tst_idx:
                assert idx not in trn_idx

            xs_trn, xs_tst, _ = centering(xs[trn_idx], xs[tst_idx])
            ys_trn = ys[trn_idx]
            ys_tst = ys[tst_idx]

            xs_tst = xs_tst[ys_tst != 0.5]
            ys_tst = ys_tst[ys_tst != 0.5]

            xs_trn = xs_trn[ys_trn != 0.5]
            ys_trn = ys_trn[ys_trn != 0.5]

            ys_trn = np.expand_dims(ys_trn, axis=1)
            ys_tst = np.expand_dims(ys_tst, axis=1)

            dataset = DataSet('dataset', TMB_DATA_PATH, 1)
            dataset.load(xs_trn, ys_trn, which='train')
            dataset.load(xs_tst, ys_tst, which='test')

            n_train_batches = int(xs_trn.shape[0] / args.batch_size)

            mlp = MLP(dataset, [2048, 500, 200, 100], logger)
            mlp.build()
            mlp.classify()
            mlp.build_loss()

            global_step = tf.train.get_or_create_global_step()
            bdrs = [int(n_train_batches*args.n_epochs*r) for r in [.5, .8]]
            vals = [lr_rate, lr_rate*1e-1, lr_rate*1e-2]
            lr = get_staircase_lr(global_step, bdrs, vals)

            update_ops = tf.get_collection(tf.GraphKeys.UPDATE_OPS)
            with tf.control_dependencies(update_ops):
                train_op = tf.train.AdamOptimizer(lr).minimize(mlp.loss,
                                global_step=global_step)

            saver = tf.train.Saver()
            with tf.Session() as sess:
                sess.run(tf.global_variables_initializer())
                for epoch in range(0, args.n_epochs):
                    cents = []
                    yhats = []
                    ys_trn = []

                    for batch_xs, batch_ys in dataset.get_next_batch(args.batch_size):
                        cent_, yhat_, _ = sess.run([mlp.cent, mlp.yhats, train_op],
                                feed_dict={mlp.xs: batch_xs, mlp.ys: batch_ys})
                        cents.append(cent_)
                        yhats.append(yhat_)
                        ys_trn.append(batch_ys)

                    yhats = np.squeeze(np.concatenate(yhats, axis=0))
                    ys_trn = np.squeeze(np.concatenate(ys_trn, axis=0))
                    auc_trn = auroc(ys_trn, yhats, pos_label=1)
                    cent_trn = np.mean(cents)

                    yhats = sess.run([mlp.yhats], feed_dict={mlp.xs: xs_tst})
                    yhats = np.squeeze(yhats)
                    auc_tst = auroc(ys_tst, yhats, pos_label=1)

                saver.save(sess, os.path.join(args.results_path,
                            '%d_%d.ckpt'%(rep, k)))

            aucs_trn.append(auc_trn)
            aucs_tst.append(auc_tst)

            logger.info("Ex%d::%d] MLP: TRAIN AUC = %.5f, TEST AUC = %.5f" %
                    (rep, k, auc_trn, auc_tst))

            tf.reset_default_graph()

        logger.info("Ex%d, TRAIN_AUC (MACRO) = %.5f, TEST AUC (MACRO) = %.5f" %
                    (rep, np.mean(aucs_trn), np.mean(aucs_tst)))

def test():
    aucs_trn = []
    aucs_tst = []

    for rep in range(n_exps):
        idx_list = cv_generate_index(rep, n_folds, ys, CV_PATH_DIR)
        for k in range(n_folds):
            logger.info("Ex%d - %dth fold" % (rep, k))

            tst_idx = [int(idx) for idx in idx_list[k]]
            trn_idx = [idx for idx in range(n_bags) if idx not in tst_idx]

            for idx in tst_idx:
                assert idx not in trn_idx

            xs_trn, xs_tst, _ = centering(xs[trn_idx], xs[tst_idx])
            ys_trn = ys[trn_idx]
            ys_tst = ys[tst_idx]

            xs_tst = xs_tst[ys_tst != 0.5]
            ys_tst = ys_tst[ys_tst != 0.5]

            xs_trn = xs_trn[ys_trn != 0.5]
            ys_trn = ys_trn[ys_trn != 0.5]

            ys_trn = np.expand_dims(ys_trn, axis=1)
            ys_tst = np.expand_dims(ys_tst, axis=1)

            dataset = DataSet('dataset', TMB_DATA_PATH, 1)
            dataset.load(xs_trn, ys_trn, which='train')
            dataset.load(xs_tst, ys_tst, which='test')

            n_train_batches = int(xs_trn.shape[0] / args.batch_size)

            mlp = MLP(dataset, [2048, 500, 200, 100], logger)
            mlp.build()
            mlp.classify()
            mlp.build_loss()

            saver = tf.train.Saver()
            with tf.Session() as sess:
                sess.run(tf.global_variables_initializer())
                saver.restore(sess, os.path.join(args.results_path,
                                '%d_%d.ckpt'%(rep, k)))
                for epoch in range(0, args.n_epochs):
                    cents = []
                    yhats = []
                    ys_trn = []

                    for batch_xs, batch_ys in dataset.get_next_batch(args.batch_size):
                        cent_, yhat_ = sess.run([mlp.cent, mlp.yhats],
                                feed_dict={mlp.xs: batch_xs, mlp.ys: batch_ys})
                        cents.append(cent_)
                        yhats.append(yhat_)
                        ys_trn.append(batch_ys)

                    yhats = np.squeeze(np.concatenate(yhats, axis=0))
                    ys_trn = np.squeeze(np.concatenate(ys_trn, axis=0))
                    auc_trn = auroc(ys_trn, yhats, pos_label=1)
                    cent_trn = np.mean(cents)

                    yhats = sess.run([mlp.yhats], feed_dict={mlp.xs: xs_tst})
                    yhats = np.squeeze(yhats)
                    auc_tst = auroc(ys_tst, yhats, pos_label=1)

            aucs_trn.append(auc_trn)
            aucs_tst.append(auc_tst)

            logger.info("Ex%d::%d] MLP: TRAIN AUC = %.5f, TEST AUC = %.5f" %
                    (rep, k, auc_trn, auc_tst))

            tf.reset_default_graph()

    logger.info("OVERALL, TRAIN_AUC (MACRO) = %.5f, TEST AUC (MACRO) = %.5f" %
                    (np.mean(aucs_trn), np.mean(aucs_tst)))

if __name__ == '__main__':
    if args.mode == 'train':
        train()
    elif args.mode == 'test':
        test()
    else:
        raise NotImplementedError()
