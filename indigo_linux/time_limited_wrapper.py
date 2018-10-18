# -*- coding: utf-8 -*-
"""
__title__ = 'time_limited_wrapper'
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
import ctypes
import inspect
import threading
from functools import wraps
from threading import Thread

from indigo import IndigoException


def time_limited_func(func, *args, **kargs):
    interval = 10

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


if __name__ == "__main__":
    f(3)
    f(1)
