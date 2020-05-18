import pickle
import numpy as np
from sklearn import svm

class SVM(object):
    def __init__(self, kernel_str='linear', fixed_labels=None,
                    max_iter=2000, verbose=False):
        if kernel_str == 'linear':
            self.classifer = svm.SVC(kernel='linear', probability=True,
                    C=0.01, max_iter=max_iter, class_weight='balanced',
                    verbose=verbose)
        else:
            self.classifer = svm.SVC(kernel='rbf', gamma='scale',
                    probability=True, class_weight='balanced',
                    max_iter=max_iter, verbose=verbose)

        self.fixed_labels = fixed_labels

    def train(self, X, y):
        self.classifer.fit(X, y)

        trnMat = self.classifer.predict_log_proba(X)

        if self.fixed_labels is not None:
            idx = np.where(self.classifer.classes_ == self.fixed_labels[1])

            trnVals = trnMat[:, idx[0]]

        return np.exp(trnVals)

    def test(self, X):
        tstMat = self.classifer.predict_log_proba(X)

        if self.fixed_labels is not None:
            idx = np.where(self.classifer.classes_ == self.fixed_labels[1])

            tstVals = tstMat[:, idx[0]]

        return np.exp(tstVals)

    def match_index(self):
        index = []

        for idx in self.fixed_labels:
            index.append(np.where(self.classifer.classes_ == idx))

        return np.hstack(index)

    def save(self, save_path):
        with open(save_path, 'wb') as f:
            pickle.dump(self.classifer, f, pickle.HIGHEST_PROTOCOL)

    def restore(self, save_path):
        with open(save_path, 'rb') as f:
            self.classifer = pickle.load(f)
