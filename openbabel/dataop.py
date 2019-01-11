# -*- coding: utf-8 -*-
"""
__title__ = 'dataop'
__package = ''
__project__ ='ChemNet"
__author__ = 'gsh'
__date_ = '2018.12.04'
__IDE__ = 'PyCharm'
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

import xlrd

data_dir = os.path.join("..", "obdata")


def read_file(file_path):
    print("Reading file from {}...".format(file_path))
    name = os.path.basename(file_path).split(".")[0]
    ExcelFile = xlrd.open_workbook(file_path)
    sheet = ExcelFile.sheet_by_index(0)
    col_names = sheet.row_values(0)
    osmi = sheet.col_values(col_names.index("Oraginal_Smiles"))[1:]
    smi = sheet.col_values(col_names.index("Smiles"))[1:]
    l = list()
    for o, s in zip(osmi, smi):
        if o and s and len(o) > 0 and len(s) > 0:
            l.append((o, s))
    return l


def read_all_data(dir_path=data_dir, num=None):
    files = [path if path.split(".")[-1] in ("xls", "xlsx") else None for path in os.listdir(dir_path)]
    file_paths = [os.path.join(dir_path, f) for f in files]
    data_list = []
    if num and num > 0:
        count = 0
        for file_path in file_paths:
            data_list.append(read_file(file_path))
            count += 1
            if count >= num:
                break
    else:
        for file_path in file_paths:
            data_list.append(read_file(file_path))
    return data_list


if __name__ == "__main__":
    data = read_all_data(num=1)
    pass
