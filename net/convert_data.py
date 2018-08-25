# -*- coding: utf-8 -*-
"""
__title__ = ''
__author__ = 'duolaoa'
__date_ = '2018.08.02'
# code is far away from bugs with the god animal protecting
    I love animals. They taste delicious.
             ┏┓   ┏┓
            ┏┛┻━━━┛┻┓
            ┃       ┃
            ┃┳┛  ┗┳ ┃
            ┃┻    ┻ ┃
            ┗━┓   ┏━┛
              ┃   ┗━━━┓
              ┃神兽保佑 ┣┓
              ┃永无BUG ┏┛
              ┗┓┓┏━┳┓┏┛
               ┃┫┫ ┃┫┫
               ┗┻┛ ┗┻┛
"""

import os
import math
import xlrd
import numpy as np

data_dir = os.path.sep.join(["..", "data"])
data_dir = os.path.abspath(data_dir)


def read_file(file_path):
    print("Reading file from {}...".format(file_path))
    name = os.path.basename(file_path).split(".")[0]
    ExcelFile = xlrd.open_workbook(file_path)
    sheet = ExcelFile.sheet_by_index(0)
    col_names = sheet.row_values(0)
    pf = sheet.col_values(col_names.index("Product1_ECFP4"))[1:]
    rf = sheet.col_values(col_names.index("Reaction_site_ECFP4"))[1:]
    label = [1 for _ in range(len(pf))]
    pf = [[int(c) for c in fingerprint] for fingerprint in pf]
    rf = [[int(c) for c in fingerprint] for fingerprint in rf]
    d = dict()
    d['name'] = name
    d['product_fingerprint'] = pf
    d['reaction_fingerprint'] = rf
    d['label'] = label
    print("Completed.")
    return d


def read_all_data(dir_path=data_dir):
    files = [path if path.split(".")[-1] in ("xls", "xlsx") else None for path in os.listdir(dir_path)]
    file_paths = [os.path.join(dir_path, f) for f in files]
    data_list = []
    for file_path in file_paths:
        data_list.append(read_file(file_path))
    return data_list


def get_train_and_test(data, train_ratio=0.7):
    train, test = [], []
    for l in data:
        pos = math.floor(len(l["product_fingerprint"]) * train_ratio)
        train.append(dict(name=l["name"], product_fingerprint=l["product_fingerprint"][:pos],
                          reaction_fingerprint=l["reaction_fingerprint"][:pos],
                          label=l["label"][:pos]))
        test.append(dict(name=l["name"], product_fingerprint=l["product_fingerprint"][pos:],
                         reaction_fingerprint=l["reaction_fingerprint"][pos:],
                         label=l["label"][pos:]))
    return train, test


def get_negative_data(train, test, shift_size=1, change_property="product_fingerprint"):
    fake_train, fake_test = [], []
    if change_property == "product_fingerprint":
        for t1 in train:
            fake_train.append(dict(name=t1["name"],
                                   product_fingerprint=t1["product_fingerprint"][shift_size:] + t1[
                                                                                                    "product_fingerprint"][
                                                                                                :shift_size],
                                   reaction_fingerprint=t1["reaction_fingerprint"],
                                   label=[0 if _ == 1 else 1 for _ in t1["label"]]))
        for t2 in test:
            fake_test.append(dict(name=t2["name"],
                                  product_fingerprint=t2["product_fingerprint"][shift_size:] + t2[
                                                                                                   "product_fingerprint"][
                                                                                               :shift_size],
                                  reaction_fingerprint=t2["reaction_fingerprint"],
                                  label=[0 if _ == 1 else 1 for _ in t2["label"]]))
        return fake_train, fake_test
    elif change_property == "reaction_fingerprint":
        for t1 in train:
            fake_train.append(dict(name=t1["name"],
                                   product_fingerprint=t1["product_fingerprint"],
                                   reaction_fingerprint=t1["reaction_fingerprint"][shift_size:] + t1[
                                                                                                      "reaction_fingerprint"][
                                                                                                  :shift_size],
                                   label=[0 if _ == 1 else 1 for _ in t1["label"]]))
        for t2 in test:
            fake_test.append(dict(name=t2["name"],
                                  product_fingerprint=t2["product_fingerprint"],
                                  reaction_fingerprint=t2["reaction_fingerprint"][shift_size:] + t2[
                                                                                                     "reaction_fingerprint"][
                                                                                                 :shift_size],
                                  label=[0 if _ == 1 else 1 for _ in t2["label"]]))
        return fake_train, fake_test
        # elif change_property == "label":
        #     for t1 in train:
        #         fake_train.append(dict(name=t1["name"],
        #                                product_fingerprint=t1["product_fingerprint"],
        #                                reaction_fingerprint=t1["reaction_fingerprint"],
        #                                label=t1["label"][shift_size:] + t1["label"][:shift_size]))
        #     for t2 in test:
        #         fake_test.append(dict(name=t2["name"],
        #                               product_fingerprint=t2["product_fingerprint"],
        #                               reaction_fingerprint=t2["reaction_fingerprint"],
        #                               label=t2["label"][shift_size:] + t2["label"][:shift_size]))
        # return fake_train, fake_test
    raise (Exception("Parameter error"))


def get_shuffled_data(data=None):
    if data is None:
        data = read_all_data()
    train, test = get_train_and_test(data)
    neg_train, neg_test = get_negative_data(train, test)
    all_train, all_test = [], []
    for sheet in train:
        for couple in zip(sheet["product_fingerprint"], sheet["reaction_fingerprint"], sheet["label"]):
            all_train.append(couple)
    for sheet in neg_train:
        for couple in zip(sheet["product_fingerprint"], sheet["reaction_fingerprint"], sheet["label"]):
            all_train.append(couple)
    for sheet in test:
        for couple in zip(sheet["product_fingerprint"], sheet["reaction_fingerprint"], sheet["label"]):
            all_test.append(couple)
    for sheet in neg_test:
        for couple in zip(sheet["product_fingerprint"], sheet["reaction_fingerprint"], sheet["label"]):
            all_test.append(couple)
    shuffled_train = np.random.permutation(all_train)
    shuffled_test = np.random.permutation(all_test)
    final_train = dict(product_fingerprint=[], reaction_fingerprint=[], label=[])
    final_test = dict(product_fingerprint=[], reaction_fingerprint=[], label=[])
    for triple in shuffled_train:
        final_train["product_fingerprint"].append(triple[0])
        final_train["reaction_fingerprint"].append(triple[1])
        final_train["label"].append(triple[2])
    for triple in shuffled_test:
        final_test["product_fingerprint"].append(triple[0])
        final_test["reaction_fingerprint"].append(triple[1])
        final_test["label"].append(triple[2])
    return final_train, final_test


def get_shuffled_data2(data=None):
    if data is None:
        data = read_all_data()
    train, test = get_train_and_test(data)
    neg_train, neg_test = get_negative_data(train, test)
    all_train, all_test = [], []
    for sheet in train:
        for couple in zip(sheet["product_fingerprint"], sheet["reaction_fingerprint"], sheet["label"]):
            all_train.append(couple)
    for sheet in neg_train:
        for couple in zip(sheet["product_fingerprint"], sheet["reaction_fingerprint"], sheet["label"]):
            all_train.append(couple)
    for sheet in test:
        for couple in zip(sheet["product_fingerprint"], sheet["reaction_fingerprint"], sheet["label"]):
            all_test.append(couple)
    for sheet in neg_test:
        for couple in zip(sheet["product_fingerprint"], sheet["reaction_fingerprint"], sheet["label"]):
            all_test.append(couple)
    shuffled_train = np.random.permutation(all_train)
    shuffled_test = np.random.permutation(all_test)
    final_train = dict(product_fingerprint=[], reaction_fingerprint=[], spliced_fingerprint=[], label=[])
    final_test = dict(product_fingerprint=[], reaction_fingerprint=[], spliced_fingerprint=[], label=[])
    for triple in shuffled_train:
        final_train["product_fingerprint"].append(triple[0])
        final_train["reaction_fingerprint"].append(triple[1])
        final_train["spliced_fingerprint"].append(triple[0] + triple[1])
        final_train["label"].append(triple[2])
    for triple in shuffled_test:
        final_test["product_fingerprint"].append(triple[0])
        final_test["reaction_fingerprint"].append(triple[1])
        final_test["spliced_fingerprint"].append(triple[0] + triple[1])
        final_test["label"].append(triple[2])

    return final_train, final_test


if __name__ == "__main__":
    # data_dir = os.path.sep.join(["..", "data"])
    # data_dir = os.path.abspath(data_dir)
    # l = read_all_data()
    data = [read_file(os.path.join("..", "data", "Suzuki_2.xlsx"))]
    # train, test = get_train_and_test(data)
    # train = [dict(name="name",
    #               product_fingerprint=["00", "01", "10"],
    #               reaction_fingerprint=["01", "10", "11"],
    #               label=[1, 1, 0]),
    #          dict(name="name",
    #               product_fingerprint=["00", "01", "10"],
    #               reaction_fingerprint=["01", "10", "11"],
    #               label=[1, 1, 0])
    #          ]
    # test = [dict(name="name",
    #              product_fingerprint=["00", "01", "10"],
    #              reaction_fingerprint=["01", "10", "11"],
    #              label=[1, 1, 0]),
    #         dict(name="name",
    #              product_fingerprint=["00", "01", "10"],
    #              reaction_fingerprint=["01", "10", "11"],
    #              label=[1, 1, 0])
    #         ]
    t1, t2 = get_shuffled_data(data)
    print(t1, t2)
    pass
