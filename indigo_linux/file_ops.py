# -*- coding: utf-8 -*-
"""
__title__ = 'file_ops'
__package = ''
__project__ ='ChemNet"
__author__ = 'gsh'
__date_ = '2018.10.18'
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

data_dir = os.path.sep.join(["..", "data"])
data_dir = os.path.abspath(data_dir)

automap_dir = os.path.join(data_dir, "..", "automap")
automap_dir = os.path.abspath(automap_dir)

error_dir = os.path.join(data_dir, "..", "automap_error")
error_dir = os.path.abspath(error_dir)


def get_files_in(dir_path, suffix):
    if not os.path.isdir(dir_path):
        os.mkdir(dir_path)
    files = []
    for path in os.listdir(dir_path):
        if path.split(".")[-1] in suffix:
            files.append(path)
        else:
            print(path)
    return files


def get_unhandle_files():
    return get_files_in(data_dir, ("xls", "xlsx"))


def get_handled_files():
    return get_files_in(automap_dir, ("csv",))


def read_file(file_path):
    print("Reading file from {}...".format(file_path))
    name = os.path.basename(file_path).split(".")[0]
    ExcelFile = xlrd.open_workbook(file_path)
    sheet = ExcelFile.sheet_by_index(0)
    col_names = sheet.row_values(0)
    assert "Reaction_Smiles" in col_names
    reactions = sheet.col_values(col_names.index("Reaction_Smiles"))[1:]
    print("Completed.")
    return reactions


def read_all_data(dir_path, num=None):
    files = get_files_in(data_dir)
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


def read_certain_file(dir_path, id):
    files = get_files_in(data_dir)
    file_paths = [os.path.join(dir_path, f) for f in files]
    assert id < len(file_paths)
    file_path = file_paths[id]
    data = read_file(file_path)
    return data, file_path


def get_unhandled_file():
    handled_files = get_files_in(automap_dir, ["csv"])
    handled_files_name = [file[:file.rindex(".")] for file in handled_files]
    error_files = get_files_in(error_dir, ["csv"])
    error_files_name = [file[:file.rindex(".")].lstrip("error_") for file in error_files]
    unhandled_files = get_files_in(data_dir, ("xls", "xlsx"))
    handled_files_name.extend(error_files_name)
    print("All files num:{:4d}".format(len(unhandled_files)))
    to_remove_list = []
    for unhandled_file in unhandled_files:
        unhandled_file_name = unhandled_file[:unhandled_file.rindex(".")]
        if unhandled_file_name in handled_files_name:
            to_remove_list.append(unhandled_file)
    for file in to_remove_list:
        unhandled_files.remove(file)
    print("Unhandled files num:{:4d}".format(len(unhandled_files)))
    print("Handled files num:{:4d}".format(len(handled_files)))
    return unhandled_files


def files_to_basename(files):
    return [file[:file.rindex(".")] for file in files]


def basename_to_abspath(names, dir_path):
    return [os.path.join(dir_path, name) for name in names]


def add_suffix(file_path):
    if os.path.isfile(file_path + ".xls"):
        file_path += ".xls"
    elif os.path.isfile(file_path + ".xlsx"):
        file_path += ".xlsx"
    else:
        raise FileNotFoundError("No such file in " + file_path)
    return file_path


class Recorder(object):
    def __init__(self):
        self.corruptlist = "corruptlist.txt"
        self.donelist = "donelist.txt"

    def _get_list(self, path):
        file = open(path)
        l = []
        for line in file:
            l.append(line)
        file.close()
        return l

    def get_corrupt(self):
        return self._get_list(self.corruptlist)

    def get_done(self):
        return self._get_list(self.donelist)

    def append_corrupt(self, l):
        file = open(self.corruptlist, mode="a+")
        for line in l:
            file.write(line + "\n")
        file.close()

    def write_corrupt(self, l):
        file = open(self.corruptlist, mode="w")
        for line in l:
            if line.endswith("\n"):
                file.write(line)
            else:
                file.write(line + "\n")
        file.close()

    def write_done(self, l):
        file = open(self.donelist, mode="w")
        for line in l:
            if line.endswith("\n"):
                file.write(line)
            else:
                file.write(line + "\n")
        file.close()

    def maintain_done(self):
        handled = get_handled_files()
        handled = files_to_basename(handled)
        self.write_done(handled)

    def get_no_corrupt(self, l):
        corrupt = self.get_corrupt()
        corrupt = [c.strip("\n") for c in corrupt]
        for c in corrupt:
            if c in l:
                l.remove(c)
        return l

    def maintain_corrupt(self):
        corruption = self.get_corrupt()
        done = self.get_done()
        to_remove = []
        for c in corruption:
            if c in done:
                to_remove.append(c)
        for r in to_remove:
            corruption.remove(r)
        corruption = list(set(corruption))
        self.write_corrupt(corruption)


if __name__ == "__main__":
    r = Recorder()
    r.maintain_done()
    r.maintain_corrupt()
