import tensorflow as tf
import os
import numpy as np

class DataSet(object):
    def __init__(self, name, data_dir, n_classes):
        self.name = name
        self.data_dir = data_dir
        self.n_classes = n_classes

    def data_selection(self, which='train'):
        if which == 'train':
            data, labels = self.train_data, self.train_labels
        elif which == 'test':
            data, labels = self.test_data, self.test_labels
        elif which == 'val':
            data, labels = self.val_data, self.val_labels
        else:
            raise NotImplementedError()
        return data, labels

    def load(self, xs, ys, which='train'):
        if which == 'train':
            self.train_data    =  xs
            self.train_labels  =  ys
        elif which == 'test':
            self.test_data     =  xs
            self.test_labels   =  ys
        elif which == 'val':
            self.val_data      =  xs
            self.val_labels    =  ys
        else:
            raise NotImplementedError()

    def get_next_batch(self, batch_size, shuffle=True, which='train'):
        data, labels = self.data_selection(which)

        batch_idx = 0
        idxs = np.arange(0, len(data))

        if shuffle:
            np.random.shuffle(idxs)

        for batch_idx in range(0, len(data), batch_size):
            cur_idxs = idxs[batch_idx:batch_idx + batch_size]
            data_batch = data[cur_idxs]
            labels_batch = labels[cur_idxs]

            yield data_batch, labels_batch

    def print_stats(self, logger, is_val):
        if is_val:
            train, val, test = self.train_data, self.val_data, self.test_data
            logger.info('[*] #inst (train, val, test): (%d, %d, %d)'
                            %(train.shape[0], val.shape[0], test.shape[0]))
        else:
            train, test = self.train_data, self.test_data
            logger.info('[*] #inst (train, test): (%d, %d)'
                            %(train.shape[0], test.shape[0]))
        logger.info('[*] Printing simple stats is done...')
