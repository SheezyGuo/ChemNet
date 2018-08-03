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

data_dir = os.path.sep.join(["..", "data"])
data_dir = os.path.abspath(data_dir)


def read_file(file_path):
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
        pos = math.floor(len(l["product_fingerprint"]) * 0.7)
        train.append(dict(name=l["name"], product_fingerprint=l["product_fingerprint"][:pos],
                          reaction_fingerprint=l["reaction_fingerprint"][:pos],
                          label=l["label"][:pos]))
        test.append(dict(name=l["name"], product_fingerprint=l["product_fingerprint"][pos:],
                         reaction_fingerprint=l["reaction_fingerprint"][pos:],
                         label=l["label"][pos:]))
    return train, test


if __name__ == "__main__":
    data_dir = os.path.sep.join(["..", "data"])
    data_dir = os.path.abspath(data_dir)
    l = read_all_data()
    pass
