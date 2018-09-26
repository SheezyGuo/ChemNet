import threading

import pandas as pd
import xlrd

from indigo import *


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
    files = [path if path.split(".")[-1] in ("xls", "xlsx") else None for path in os.listdir(dir_path)]
    file_paths = [os.path.join(dir_path, f) for f in files]
    assert id < len(file_paths)
    file_path = file_paths[id]
    data = read_file(file_path)
    return data, file_path


class myThread(threading.Thread):
    def __init__(self, num):
        threading.Thread.__init__(self)
        self.num = num

    def run(self):
        print("Starting " + self.name)
        data_dir = os.path.join("..", "data")
        data, file_path = read_certain_file(data_dir, id=self.num)
        indigo = Indigo()
        raw_reactions = []
        transformed_reactions = []
        automap_list = []
        reaction_list = data
        count = 0
        for reaction in reaction_list:
            count += 1
            if count % 100 == 0:
                print(self.name, count)
            try:
                rxn = indigo.loadReaction(reaction)
                transformed = rxn.smiles()
                rxn.automap("discard")
                automap = rxn.smiles()
                raw_reactions.append(reaction)
                transformed_reactions.append(transformed)
                automap_list.append(automap)
            except IndigoException as e:
                print(self.name, reaction, e)
            except Exception as e:
                print(e)
                continue
        dataframe = pd.DataFrame({
            "raw_ractions": raw_reactions, "transformed_reactions": transformed_reactions, "automap": automap_list})
        if not os.path.isdir(os.path.join("..", "automap")):
            os.mkdir(os.path.join("..", "automap"))
        dataframe.to_csv(
            "{}.csv".format(os.path.join("..", "automap", os.path.basename(file_path).split(".")[0])), sep=',')
        print("Exiting " + self.name)


if __name__ == "__main__":
    file_num = 18
    file_list = [2, 7]
    thread_list = []
    for i in file_list:
        t = myThread(i)
        thread_list.append(t)
    for t in thread_list:
        t.start()
