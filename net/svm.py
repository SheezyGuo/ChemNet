# -*- coding: utf-8 -*-
"""
__title__ = ''
__author__ = 'duolaoa'
__date_ = '2018.08.26'
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

import numpy as np
from sklearn.svm import SVC
from convert_data import *

train, test = get_shuffled_data2()
clf = SVC()
clf.fit(train["spliced_fingerprint"], train["label"])
result = clf.predict(test["spliced_fingerprint"])
equals = np.equal(result, test["label"])
mean = np.mean(equals)
print(mean)
