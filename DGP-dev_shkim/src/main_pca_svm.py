import os
import glob
import numpy as np
import argparse
import logging

from sklearn.decomposition import PCA
from utils.utils import cv_generate_index, load_tmb_files_v2, centering, auroc
from utils.paths import TMB_DATA_PATH, TMB_INFO_FNAME, CV_PATH_DIR
from models.svm import SVM

np.random.seed(6150)

parser = argparse.ArgumentParser()
parser.add_argument('--results_path', type=str,
                    default='../results/pca_svm')
parser.add_argument('--n_folds', type=int, default=5)
parser.add_argument('--n_exps', type=str, default=5)
parser.add_argument('--mode', type=str, default='train')

args = parser.parse_args()

if not os.path.exists(args.results_path):
    os.makedirs(args.results_path)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('pca_svm_tmb')
logger.addHandler(logging.FileHandler(
    os.path.join(args.results_path, '%s.log'%args.mode), mode='w')
)

n_exps = args.n_exps
n_folds = args.n_folds

xs, ys = load_tmb_files_v2(TMB_DATA_PATH, TMB_INFO_FNAME)
xs = [np.mean(xs_, axis=0, keepdims=True) for xs_ in xs]
xs = np.concatenate(xs, axis=0)
n_bags = xs.shape[0]

pca = PCA(n_components=100)
svm = SVM('rbf', np.array([-1, 1]), max_iter=2000, verbose=False)

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

            xs_tst = xs_tst[ys_tst != 0]
            ys_tst = ys_tst[ys_tst != 0]
            xs_trn = xs_trn[ys_trn != 0]
            ys_trn = ys_trn[ys_trn != 0]

            pca.fit(xs_trn)
            xs_trn_pca = pca.transform(xs_trn)
            xs_tst_pca = pca.transform(xs_tst)

            ys_pred_trn = svm.train(xs_trn_pca, ys_trn)
            ys_pred_tst = svm.test(xs_tst_pca)

            auc_tst = auroc(ys_tst, ys_pred_tst, pos_label=1)
            auc_trn = auroc(ys_trn, ys_pred_trn, pos_label=1)

            aucs_tst.append(auc_tst)
            aucs_trn.append(auc_trn)

            svm.save(os.path.join(args.results_path,
                            'model_%d_%d.bin'%(rep, k)))

            logger.info("Ex%d::%d] SVM: TRAIN AUC = %.5f, TEST AUC = %.5f" %
                    (rep, k, auc_trn, auc_tst))
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

            xs_tst = xs_tst[ys_tst != 0]
            ys_tst = ys_tst[ys_tst != 0]
            xs_trn = xs_trn[ys_trn != 0]
            ys_trn = ys_trn[ys_trn != 0]

            pca.fit(xs_trn)
            xs_trn_pca = pca.transform(xs_trn)
            xs_tst_pca = pca.transform(xs_tst)

            svm.restore(os.path.join(args.results_path,
                                      'model_%d_%d.bin'%(rep, k)))
            ys_pred_trn = svm.test(xs_trn_pca)
            ys_pred_tst = svm.test(xs_tst_pca)

            auc_tst = auroc(ys_tst, ys_pred_tst, pos_label=1)
            auc_trn = auroc(ys_trn, ys_pred_trn, pos_label=1)

            aucs_tst.append(auc_tst)
            aucs_trn.append(auc_trn)

            logger.info("Ex%d::%d] SVM: TRAIN AUC = %.5f, TEST AUC = %.5f" %
                            (rep, k, auc_trn, auc_tst))
    logger.info("OVERALL, TRAIN_AUC (MACRO) = %.5f, TEST AUC (MACRO) = %.5f" %
                    (np.mean(aucs_trn), np.mean(aucs_tst)))

if __name__ == '__main__':
    if args.mode == 'train':
        train()
    elif args.mode == 'test':
        test()
    else:
        raise NotImplementedError()
