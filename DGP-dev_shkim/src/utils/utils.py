import os
import glob
import numpy as np
from scipy.io import loadmat
from sklearn.metrics import roc_curve, auc

def auroc(ytrue, pred, pos_label=1):
    fpr, tpr, th = roc_curve(ytrue, pred, pos_label=pos_label)
    _auc = auc(fpr, tpr)
    return _auc

def centering(xmat_trn, xmat_tst):
    x_mean = np.mean(xmat_trn, axis=0)
    xmat_trn -= x_mean
    xmat_tst -= x_mean

    return xmat_trn, xmat_tst, x_mean

def cv_generate_index(num_Exps, num_Kfold, labels, path_dir):
    check_flag = os.path.exists(path_dir + 'Ex_' + str(num_Exps))

    if check_flag == True:
        idx_list = []

        for f in glob.glob1(path_dir + 'Ex_' + str(num_Exps) + "/", "*.txt"):
            fold_idx = []
            with open(path_dir + 'Ex_' + str(num_Exps) + '/' + f) as ftxt:
                for l in ftxt:
                    fold_idx.append(l.strip())

            idx_list.append(np.hstack(fold_idx))
    else:
        os.mkdir(path_dir + 'Ex_' + str(num_Exps) + '/')

        num_total = np.size(labels)
        num_infold = np.int_(np.ceil(num_total / num_Kfold))

        pos_idx = np.reshape(np.asarray(np.where(np.equal(labels, 1))),  [-1])
        neg_idx = np.reshape(np.asarray(np.where(np.equal(labels, -1))), [-1])

        n_pos_tot = np.size(pos_idx)
        n_pos_ink = np.int_(np.ceil((n_pos_tot/num_total) * num_infold))

        n_neg_tot = np.size(neg_idx)
        n_neg_ink = num_infold - n_pos_ink

        pos_idx_rnd = pos_idx[np.random.choice(n_pos_tot, n_pos_tot, replace=False)]
        neg_idx_rnd = neg_idx[np.random.choice(n_neg_tot, n_neg_tot, replace=False)]

        idx_list = []
        for k in range(num_Kfold):
            if k == num_Kfold-1:
                pos_endidx = n_pos_tot
                neg_endidx = n_neg_tot
            else:
                pos_endidx = (k + 1) * n_pos_ink
                neg_endidx = (k + 1) * n_neg_ink

            idx1 = np.arange(k*n_pos_ink, pos_endidx)
            idx2 = np.arange(k*n_neg_ink, neg_endidx)

            tst_idx = np.sort(np.concatenate((pos_idx_rnd[idx1], neg_idx_rnd[idx2])))

            np.savetxt(path_dir + 'Ex_' + str(num_Exps) + '/' + str(k) + '.txt',
                       tst_idx, delimiter='\t', fmt='%d')

            idx_list.append(tst_idx)

    return idx_list

def load_tmb_files_v2(data_path, info_filename):
    cnt = 0

    file_info = []
    X_list = []
    Idx_list = []
    for f in glob.glob1(data_path, "*.mat"):
        data = loadmat(data_path + f)
        Xcur = data['image_features']

        X_list.append(Xcur)
        Idx_list.append(cnt*np.ones((np.shape(Xcur)[0], 1)))

        file_info.append([[f], [np.shape(Xcur)[0]]])
        cnt += 1

    X = np.vstack(X_list)
    Idx = np.vstack(Idx_list)
    Idx.astype(int)

    table_info = []
    with open(info_filename) as f:
        for l in f:
            table_info.append(l.strip().split("\t"))

    table_ids = [x[0] for x in table_info]

    # threshold
    SNV_num = np.vstack([float(x[1]) for x in table_info[1:]])
    InDel = np.vstack([float(x[2]) for x in table_info[1:]])

    thres_high = np.percentile(SNV_num + InDel, 75, axis=0)
    thres_low = np.percentile(SNV_num + InDel, 25, axis=0)

    cnt = 0
    summary_table = []
    y_labels = []
    for x in file_info:
        str_id = x[0][0][0:12]

        try:
            index = table_ids.index(str_id)

            y1 = float(table_info[index][1]) + float(table_info[index][2])

            if y1 >= thres_high:
                y2 = 1
            elif y1 < thres_low:
                y2 = -1
            else:
                y2 = 0

            y_labels.append(y2)

            cur_table = np.array([cnt, x[0][0], str_id, y1, y2,
                                  table_info[index][0], table_info[index][1], table_info[index][2]])

        except ValueError:
            cur_table = [cnt, x[0][0], str_id, np.nan, np.nan, '', '', '']

        summary_table.append(cur_table)
        cnt += 1

    y_labels = np.asarray(y_labels)

    assert len(X_list) == y_labels.shape[0]

    return X_list, y_labels

def load_tmb_files(data_path, info_filename):
    cnt = 0

    file_info = []
    X_list = []
    Idx_list = []
    for f in glob.glob1(data_path, "*.mat"):
        data = loadmat(data_path + f)
        Xcur = data['image_features']

        X_list.append(Xcur)
        Idx_list.append(cnt*np.ones((np.shape(Xcur)[0], 1)))

        file_info.append([[f], [np.shape(Xcur)[0]]])
        cnt += 1

    X = np.vstack(X_list)
    Idx = np.vstack(Idx_list)
    Idx.astype(int)

    table_info = []
    with open(info_filename) as f:
        for l in f:
            table_info.append(l.strip().split("\t"))

    table_ids = [x[0] for x in table_info]

    # threshold
    SNV_num = np.vstack([float(x[1]) for x in table_info[1:]])
    InDel = np.vstack([float(x[2]) for x in table_info[1:]])

    thres = np.percentile(SNV_num + InDel, 75, axis=0)

    cnt = 0
    summary_table = []
    for x in file_info:
        str_id = x[0][0][0:12]

        try:
            index = table_ids.index(str_id)

            y1 = float(table_info[index][1]) + float(table_info[index][2])
            y2 = 1 if y1 >= thres else -1

            cur_table = np.array([cnt, x[0][0], str_id, y1, y2,
                                  table_info[index][0], table_info[index][1], table_info[index][2]])

        except ValueError:
            cur_table = [cnt, x[0][0], str_id, np.nan, np.nan, '', '', '']

        summary_table.append(cur_table)
        cnt += 1

    Y = np.zeros((np.size(Idx),2), dtype=int)
    Y[:, 1] = np.reshape(Idx, [-1])
    for x in summary_table:
        cur_idx = x[0]
        chk_idx = Idx == int(cur_idx)
        Y[np.reshape(chk_idx, [-1]), 0] = int(x[4])

    return X, Y, summary_table
