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
    # data = read_all_data()
    data = [read_file("../data/test_data_2000.xlsx")]
    # final_train, final_test = get_shuffled_data(data)
    final_train, final_test = get_train_and_test(data)
    fake_train, fake_test = get_negative_data(final_train, final_test)
    final_train = final_train[0]
    final_test = final_test[0]
    fake_train = fake_train[0]
    fake_test = fake_test[0]
    net.train(product_fingerprint=final_train["product_fingerprint"],
              reaction_fingerprint=final_train["reaction_fingerprint"], label=final_train["label"])
    net.predict(product_fingerprint=final_test["product_fingerprint"],
                reaction_fingerprint=final_test["reaction_fingerprint"], label=final_test["label"])
    print("*****")
    net.predict(product_fingerprint=fake_test["product_fingerprint"],
                reaction_fingerprint=fake_test["reaction_fingerprint"], label=fake_test["label"])
    net.train(product_fingerprint=fake_train["product_fingerprint"],
              reaction_fingerprint=fake_train["reaction_fingerprint"], label=fake_train["label"])
    net.predict(product_fingerprint=fake_test["product_fingerprint"],
                reaction_fingerprint=fake_test["reaction_fingerprint"], label=fake_test["label"])
    print("*****")
    net.predict(product_fingerprint=final_test["product_fingerprint"],
                reaction_fingerprint=final_test["reaction_fingerprint"], label=final_test["label"])
    s_train, s_test = get_shuffled_data(data)
    net.train(product_fingerprint=s_train["product_fingerprint"],
              reaction_fingerprint=s_train["reaction_fingerprint"], label=s_train["label"])
    net.train(product_fingerprint=s_train["product_fingerprint"],
              reaction_fingerprint=s_train["reaction_fingerprint"], label=s_train["label"])
    net.predict(product_fingerprint=s_test["product_fingerprint"],
                reaction_fingerprint=s_test["reaction_fingerprint"], label=s_test["label"])

def main2():
    net = SimilarityNet(True)
    final_train, final_test = get_shuffled_data()
    net.train(product_fingerprint=final_train["product_fingerprint"],
              reaction_fingerprint=final_train["reaction_fingerprint"], label=final_train["label"])
    net.predict(product_fingerprint=final_test["product_fingerprint"],
                reaction_fingerprint=final_test["reaction_fingerprint"], label=final_test["label"])


if __name__ == "__main__":
    main()
