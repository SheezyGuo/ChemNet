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
    # pf = np.random.uniform(-1, 1, [1280, 1024]).astype(np.float32)
    # rf = np.random.uniform(-1, 1, [1280, 1024]).astype(np.float32)
    # label = np.ones([1280, 1], dtype=np.float32)
    data = read_all_data()
    train, test = get_train_and_test(data)
    fake_train, fake_test = get_negative_data(train, test)
    # for sheet in train:
    #     net.train(product_fingerprint=sheet["product_fingerprint"],
    #               reaction_fingerprint=sheet["reaction_fingerprint"], label=sheet["label"])
    # for sheet in fake_train:
    #     net.train(product_fingerprint=sheet["product_fingerprint"],
    #               reaction_fingerprint=sheet["reaction_fingerprint"], label=sheet["label"])
    # for sheet in test:
    #     net.predict(product_fingerprint=sheet["product_fingerprint"],
    #                 reaction_fingerprint=sheet["reaction_fingerprint"], label=sheet["label"])
    # for sheet in fake_test:
    #     net.predict(product_fingerprint=sheet["product_fingerprint"],
    #                 reaction_fingerprint=sheet["reaction_fingerprint"], label=sheet["label"])
    for sheet in train:
        net.train(product_fingerprint=sheet["product_fingerprint"],
                  reaction_fingerprint=sheet["reaction_fingerprint"], label=sheet["label"])

    for sheet in test:
        net.predict(product_fingerprint=sheet["product_fingerprint"],
                    reaction_fingerprint=sheet["reaction_fingerprint"], label=sheet["label"])
    for sheet in fake_train:
        net.train(product_fingerprint=sheet["product_fingerprint"],
                  reaction_fingerprint=sheet["reaction_fingerprint"], label=sheet["label"])
    for sheet in fake_test:
        net.predict(product_fingerprint=sheet["product_fingerprint"],
                    reaction_fingerprint=sheet["reaction_fingerprint"], label=sheet["label"])


if __name__ == "__main__":
    main()
