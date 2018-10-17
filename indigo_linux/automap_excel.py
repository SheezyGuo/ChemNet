import ctypes
import inspect
import math
import threading
from functools import wraps
from queue import Queue
from threading import Thread

import pandas as pd
import xlrd

from convert_data import data_dir
from indigo import *


def get_excel_files_in(dir_path):
    if not os.path.isdir(dir_path):
        os.mkdir(dir_path)
    files = []
    for path in os.listdir(dir_path):
        if path.split(".")[-1] in ("xls", "xlsx"):
            files.append(path)
        else:
            print(path)
    return files


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
    files = get_excel_files_in(data_dir)
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
    files = get_excel_files_in(data_dir)
    file_paths = [os.path.join(dir_path, f) for f in files]
    assert id < len(file_paths)
    file_path = file_paths[id]
    data = read_file(file_path)
    return data, file_path


indigo = Indigo()


class MyThread(threading.Thread):

    def __init__(self, queue, count_list):
        threading.Thread.__init__(self)
        self.queue = queue
        self.count_list = count_list

    def run(self):
        global indigo
        while not self.queue.empty():
            file_path = self.queue.get()
            data = read_file(file_path)
            self.queue.task_done()
            file_name = os.path.basename(file_path)
            file_name = file_name[:file_name.rindex(".")]
            self.setName("{:s}".format(file_name))
            print("Starting " + self.name)

            raw_reactions = []
            transformed_reactions = []
            automap_list = []
            error_reactions = []
            error_info = []
            reaction_list = data

            for index, reaction in enumerate(reaction_list):
                if reaction.startswith(">>") or reaction.count(".") > 1:
                    continue
                print(self.name, index)
                try:
                    rxn = indigo.loadReaction(reaction)
                    transformed = rxn.smiles()

                    @time_limited_func
                    def rxn_auto_map(_rxn):
                        _rxn.automap("discard")

                    result, exception = rxn_auto_map(rxn)
                    if exception is not None:
                        raise Exception(str(exception))
                    if not result:
                        raise TimeoutError("Automap timeout.")
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
            if len(raw_reactions) > 0 or len(transformed_reactions) > 0 or len(automap_list) > 0:
                dataframe.to_csv(
                    "{}.csv".format(os.path.join("..", "automap", file_name), sep=','))
            if len(error_reactions) > 0 or len(error_info) > 0:
                errorframe.to_csv(
                    "{}.csv".format(os.path.join("..", "automap_error", "error_" + file_name)), sep=',')
            self.count_list[0] = self.count_list[0] + 1
            print("{:04d}/{:04d}   Exiting {:s}".format(self.count_list[0], self.count_list[1], self.name))


def main(slice=None):
    files = get_excel_files_in(data_dir)
    file_num = len(files)
    thread_num = 32
    if slice == "all":
        task_queue = Queue()
        for file in files:
            task_queue.put(file)
        thread_list = []
        count_list = [0, file_num]
        for i in range(thread_num):
            t = MyThread(task_queue, count_list)
            thread_list.append(t)
        for t in thread_list:
            t.start()
        for t in thread_list:
            t.join()
        task_queue.join()
    elif slice is None or slice == 0:
        for i in range(int(math.ceil(file_num / thread_num))):
            print("Round", i + 1, "#" * 60)
            task_queue = Queue()
            sub_files = files[i * thread_num: file_num if (i + 1) * thread_num > file_num else (i + 1) * thread_num]
            for file in sub_files:
                task_queue.put(file)
            thread_list = []
            count_list = [0, file_num]
            for i in range(thread_num):
                t = MyThread(task_queue, count_list)
                thread_list.append(t)
            for t in thread_list:
                t.start()
            for t in thread_list:
                t.join()
            task_queue.join()
    elif slice > 0:
        task_queue = Queue()
        sub_files = files[
                    slice * thread_num: file_num if (slice + 1) * thread_num > file_num else (slice + 1) * thread_num]
        for file in sub_files:
            task_queue.put(file)
        thread_list = []
        count_list = [0, file_num]
        for i in range(thread_num):
            t = MyThread(task_queue, count_list)
            thread_list.append(t)
        for t in thread_list:
            t.start()
        for t in thread_list:
            t.join()
        task_queue.join()
        print(slice)


def time_limited_func(func, *args, **kargs):
    interval = 100

    def _async_raise(tid, exctype):
        """raises the exception, performs cleanup if needed"""
        tid = ctypes.c_long(tid)
        if not inspect.isclass(exctype):
            exctype = type(exctype)
        res = ctypes.pythonapi.PyThreadState_SetAsyncExc(tid, ctypes.py_object(exctype))
        if res == 0:
            raise ValueError("invalid thread id")
        elif res != 1:
            # """if it returns a number greater than one, you're in trouble,
            # and you should call it again with exc=NULL to revert the effect"""
            ctypes.pythonapi.PyThreadState_SetAsyncExc(tid, None)
            raise SystemError("PyThreadState_SetAsyncExc failed")

    def stop_thread(thread):
        _async_raise(thread.ident, SystemExit)

    @wraps(func)
    def wrapper(*args, **kargs):
        """Do func in interval"""

        class time_limited_class(Thread):
            def __init__(self):
                Thread.__init__(self)
                self.is_finished = False
                self.exception = None

            def run(self):
                timer = threading.Timer(interval, stop_thread, [self])
                timer.start()
                try:
                    func(*args, **kargs)
                except IndigoException as e:
                    self.exception = e
                except Exception as e:
                    self.exception = e
                timer.cancel()
                self.is_finished = True

        t = time_limited_class()
        t.start()
        t.join()
        return t.is_finished, t.exception

    return wrapper


@time_limited_func
def f(a):
    from time import sleep
    sleep(a)


def get_unhandled_file():
    handled_dir = os.path.abspath(os.path.join(data_dir, "..", "automap"))
    handled_files = get_excel_files_in(handled_dir)
    handled_files_name = [file[:file.rindex(".")] for file in handled_files]
    unhandled_files = get_excel_files_in(data_dir)
    print("All files num:{:4d}".format(len(unhandled_files)))
    for unhandled_file in unhandled_files:
        unhandled_file_name = unhandled_file[:unhandled_file.rindex(".")]
        if unhandled_file_name in handled_files_name:
            unhandled_files.remove(unhandled_file)
    print("Unhandled files num:{:4d}".format(len(unhandled_files)))
    print("Handled files num:{:4d}".format(len(handled_files)))
    return unhandled_files


if __name__ == "__main__":
    try:
        result = f(1)
        print(result)
    except Exception as e:
        print(e)
    main(3)
