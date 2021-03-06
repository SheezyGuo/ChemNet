# -*- coding: utf-8 -*-
"""
__title__ = ''
__author__ = 'duolaoa'
__date_ = '2018.08.16'
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

import tensorflow as tf
import numpy as np


def test_loss():
    with tf.Session() as sess:
        sim = tf.constant([0.1, 0.4, 0.5, 0.8], shape=[4, 1], dtype=tf.float32)
        label = tf.constant([0, 0, 1, 1], shape=[4, 1], dtype=tf.float32)
        similarity_loss = tf.reduce_mean(tf.square(
            tf.maximum(sim - 0.3, np.zeros_like(sim)) * (1 - label) + tf.minimum(sim - 0.7, np.zeros_like(sim)) * label
        ))
        sess.run(tf.global_variables_initializer())
        print(sess.run(similarity_loss))


if __name__ == "__main__":
    test_loss()
