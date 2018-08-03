#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
__title__ = ''
__author__ = 'duolaoa'
__mtime__ = '2018.07.31'
# code is far away from bugs with the god animal protecting
    I love animals. They taste delicious.
             ┏┓   ┏┓
            ┏┛┻━━━┛┻┓
            ┃☃      ┃
            ┃┳┛  ┗┳ ┃
            ┃┻      ┃
            ┗━┓   ┏━┛
              ┃   ┗━━━┓
              ┃神兽保佑 ┣┓
              ┃永无BUG ┏┛
              ┗┓┓┏━┳┓┏┛
               ┃┫┫ ┃┫┫
               ┗┻┛ ┗┻┛
"""

import tensorflow as tf


def similarity_layer(self, product_fingerprint, reaction_fingerprint, train_flag):
    #############################left_layer###############################

    self.dense_left = tf.layers.dense(inputs=product_fingerprint, units=1024, activation=tf.nn.elu)
    if not train_flag:
        self.dropout = tf.layers.dropout(self.dense_left, rate=0.3)
        self.dense_left = self.dropout
        pre = self.dense_left
    else:
        pre = product_fingerprint
    self.gate_1 = tf.layers.dense(inputs=self.dense_left, units=1024, activation=tf.nn.sigmoid)
    self.highway_1 = self.gate_1 * self.dense_left + (1 - self.gate_1) * pre
    self.gate_2 = tf.layers.dense(inputs=self.highway_1, units=1024, activation=tf.nn.sigmoid)
    self.highway_2 = self.gate_2 * self.highway_1 + (1 - self.gate_2) * self.dense_left
    self.gate_3 = tf.layers.dense(inputs=self.highway_2, units=1024, activation=tf.nn.sigmoid)
    self.highway_3 = self.gate_3 * self.highway_2 + (1 - self.gate_3) * self.highway_1
    self.gate_4 = tf.layers.dense(inputs=self.highway_3, units=1024, activation=tf.nn.sigmoid)
    self.highway_4 = self.gate_4 * self.highway_3 + (1 - self.gate_4) * self.highway_2
    self.gate_5 = tf.layers.dense(inputs=self.highway_4, units=1024, activation=tf.nn.sigmoid)
    self.highway_5 = self.gate_5 * self.highway_4 + (1 - self.gate_5) * self.highway_3

    #############################left_layer###############################

    #############################right_layer###############################

    self.dense_right = tf.layers.dense(inputs=reaction_fingerprint, units=1024, activation=tf.nn.elu)

    #############################right_layer###############################

    model1 = tf.sqrt(tf.reduce_sum(tf.square(self.highway_5)))
    model2 = tf.sqrt(tf.reduce_sum(tf.square(self.dense_right)))
    multi_sum = tf.reduce_sum(tf.multiply(self.highway_5, self.dense_right))
    cosine = multi_sum / (model1 * model2)

    self.value = tf.nn.sigmoid(cosine)

    return self.value


if __name__ == "__main__":
    sess = tf.Session()
    batch_size = 128
    pf = tf.placeholder(dtype=tf.float32, shape=[batch_size, 1024], name="product_fingerprint")
    rf = tf.placeholder(dtype=tf.float32, shape=[batch_size, 1024], name="reaction_fingerprint")
    value = similarity_layer(pf, rf, True)
    initializer = tf.global_variables_initializer()
    sess.run(initializer)
    print(sess.run(value, feed_dict={""}))
