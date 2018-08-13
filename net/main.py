# -*- coding: utf-8 -*-
"""
__title__ = ''
__author__ = 'duolaoa'
__date_ = '2018.08.03'
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
from convert_data import *
from net import SimilarityNet


def main():
    net = SimilarityNet(True)
    data = read_all_data()
    final_train, final_test = get_shuffled_data(data)
    net.train(product_fingerprint=final_train["product_fingerprint"],
              reaction_fingerprint=final_train["reaction_fingerprint"], label=final_train["label"])
    net.predict(product_fingerprint=final_test["product_fingerprint"],
                reaction_fingerprint=final_test["reaction_fingerprint"], label=final_test["label"])


if __name__ == "__main__":
    main()
