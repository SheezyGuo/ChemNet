# -*- coding: utf-8 -*-
"""
__title__ = 'rdf_test'
__package = ''
__project__ ='ChemNet"
__author__ = 'gsh'
__date_ = '2018.10.23'
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

from indigo import Indigo

indigo = Indigo()

for item in indigo.iterateRDFile(os.path.join("..", "rdf", "amadori.rdf")):
    for prop in item.iterateProperties():
        try:
            print(prop.name(), ":", prop.rawData())
        except Exception as e:
            print(e)
