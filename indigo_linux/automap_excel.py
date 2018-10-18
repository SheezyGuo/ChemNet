import math
import threading
from queue import Queue

import pandas as pd

from file_ops import get_unhandled_file, Recorder, files_to_basename, add_suffix
from file_ops import read_file, data_dir
from indigo import *
from time_limited_wrapper import time_limited_func

indigo = Indigo()


class MyThread(threading.Thread):

    def __init__(self, queue, count_list):
        threading.Thread.__init__(self)
        self.queue = queue
        self.count_list = count_list

    def run(self):
        global indigo
        while not self.queue.empty():
            file_name = self.queue.get()
            file_path = os.path.join(data_dir, file_name)
            file_path = add_suffix(file_path)
            data = read_file(file_path)
            self.queue.task_done()
            self.setName("{:s}".format(file_name))
            print("Starting " + self.name)

            raw_reactions = []
            transformed_reactions = []
            automap_list = []
            error_reactions = []
            error_info = []
            reaction_list = data

            for index, reaction in enumerate(reaction_list):
                if reaction.startswith(">>") or reaction.endswith(">>") or reaction.count(".") > 1:
                    error_reactions.append(reaction)
                    error_info.append("Ignored this reaction")
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


def main(mode=None):
    files = get_unhandled_file()
    files = files_to_basename(files)
    recoder = Recorder()
    recoder.maintain_corrupt()
    files = recoder.get_no_corrupt(files)
    file_num = len(files)
    thread_num = 10
    if mode == "all":
        recoder.append_corrupt(files)
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
    elif mode is None:
        for i in range(int(math.ceil(file_num / thread_num))):
            print("Round", i + 1, "#" * 60)
            task_queue = Queue()
            sub_files = files[i * thread_num: file_num if (i + 1) * thread_num > file_num else (i + 1) * thread_num]
            recoder.append_corrupt(sub_files)
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
    elif mode >= 0:
        task_queue = Queue()
        sub_files = files[
                    mode * thread_num: file_num if (mode + 1) * thread_num > file_num else (mode + 1) * thread_num]
        recoder.append_corrupt(sub_files)
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
        print(mode)


def handle_all():
    while len(get_unhandled_file()) > 0:
        t = threading.Thread(target=lambda: main(None))
        t.start()
        t.join()


if __name__ == "__main__":
    handle_all()
