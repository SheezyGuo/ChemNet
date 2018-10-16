import threading
from queue import Queue

import pandas as pd
import xlrd

from convert_data import data_dir
from indigo import *

count_num = 0
total = 0


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


def read_certain_file(dir_path, id):
    files = []
    for path in os.listdir(dir_path):
        if path.split(".")[-1] in ("xls", "xlsx"):
            files.append(path)
    file_paths = [os.path.join(dir_path, f) for f in files]
    assert id < len(file_paths)
    file_path = file_paths[id]
    data = read_file(file_path)
    return data, file_path


class MyThread(threading.Thread):
    def __init__(self, queue):
        threading.Thread.__init__(self)
        self.queue = queue

    def run(self):
        while not self.queue.empty():
            file_path = self.queue.get()
            data = read_file(file_path)
            file_name = file_path[:file_path.rindex(".")]
            self.setName("{:04d}--->{:s}".format(self.num, file_name))
            print("Starting " + self.name)
            indigo = Indigo()
            raw_reactions = []
            transformed_reactions = []
            automap_list = []
            error_reactions = []
            error_info = []
            reaction_list = data
            for index, reaction in enumerate(reaction_list):
                print(self.name, index)
                try:
                    rxn = indigo.loadReaction(reaction)
                    transformed = rxn.smiles()
                    rxn.automap("discard")
                    automap = rxn.smiles()
                    raw_reactions.append(reaction)
                    transformed_reactions.append(transformed)
                    automap_list.append(automap)
                except IndigoException as e:
                    error_reactions.append(reaction)
                    error_info.append(str(e))
                    print(self.name, reaction, e)
                except Exception as e:
                    error_reactions.append(reaction)
                    error_info.append(str(e))
                    print(self.name, reaction, e)
            dataframe = pd.DataFrame({
                "raw_ractions": raw_reactions, "transformed_reactions": transformed_reactions, "automap": automap_list})
            errorframe = pd.DataFrame({
                "error_ractions": error_reactions, "error_info": error_info})
            if not os.path.isdir(os.path.join("..", "automap")):
                os.mkdir(os.path.join("..", "automap"))
            if not os.path.isdir(os.path.join("..", "automap_error")):
                os.mkdir(os.path.join("..", "automap_error"))
            dataframe.to_csv(
                "{}.csv".format(os.path.join("..", "automap", file_name), sep=','))
            if len(error_reactions) > 0 or len(error_info) > 0:
                errorframe.to_csv(
                    "{}.csv".format(os.path.join("..", "automap_error", "error_" + file_name)), sep=',')
            count_num += 1
            print("{:04d}/{:04d}   Exiting {:s}".format(count_num, total, self.name))
            self.queue.task_done()


def main():
    files = []
    for path in os.listdir(data_dir):
        if path.split(".")[-1] in ("xls", "xlsx"):
            files.append(path)
    task_queue = Queue()
    file_num = len(files)
    for file in files:
        task_queue.put(file)
    total = file_num
    thread_list = []
    for i in range(file_num):
        t = MyThread(task_queue)
        thread_list.append(t)
    for t in thread_list:
        t.start()
        t.join()
    task_queue.join()


if __name__ == "__main__":
    main()
