# -*- coding: utf-8 -*-
"""
__title__ = 'Siamese'
__package = ''
__project__ ='ChemNet"
__author__ = 'gsh'
__date_ = '2018.09.03'
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

import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
import time

from convert_data import get_shuffled_data2
import os

flags = tf.app.flags
flags.DEFINE_float("learning_rate", 1e-4, "Learning rate")
flags.DEFINE_string("checkpoint_dir", os.path.join("..", "checkpoint"), "dir path to checkpoint")
flags.DEFINE_float("stddev", 0.01, "Stander deviation")
flags.DEFINE_integer("seq_len", 2048, "Sequence length")
flags.DEFINE_integer("batch_size", 128, "Batch size")
flags.DEFINE_float("Q", 0.5, "Value of Q")
flags.DEFINE_float("train_drop_rate", 0.3, "Training drop rate")
flags.DEFINE_float("predict_drop_rate", 0.0, "Predict drop rate")
FLAGS = flags.FLAGS


class SiameseNet(object):
    def __init__(self):
        with tf.name_scope('parameter') as scope:
            self.drop_rate = tf.placeholder(dtype=tf.float32, shape=[1, 1], name="drop_rate")
            self.global_step = tf.Variable(initial_value=0, dtype=tf.int32, name="global_step", trainable=False)
        with tf.name_scope("data") as scope:
            self.data1 = tf.placeholder(dtype=tf.float32, shape=[None, FLAGS.seq_len], name="data1")
            self.data2 = tf.placeholder(dtype=tf.float32, shape=[None, FLAGS.seq_len], name="data2")
            # self.data3 = tf.placeholder(dtype=tf.float32, shape=[None, FLAGS.seq_len], name="data3")
            self.label = tf.placeholder(dtype=tf.float32, shape=[None, 1], name="label")
        config = tf.ConfigProto()
        config.gpu_options.allow_growth = True
        self.sess = tf.Session(config=config)
        self.loss = self.__loss()
        self.optimizer = tf.train.AdamOptimizer(learning_rate=FLAGS.learning_rate)
        self.train_op = self.optimizer.minimize(self.loss, global_step=self.global_step)
        try:
            self.sess.run(tf.global_variables_initializer())
        except Exception:
            self.sess.run(tf.initialize_all_variables())
        # print(tf.trainable_variables())

    def __core(self, inputs):
        with tf.name_scope('fc1') as scope:
            fc1 = tf.layers.dense(inputs=inputs, units=1024, activation=tf.nn.elu,
                                  kernel_initializer=tf.truncated_normal_initializer(stddev=FLAGS.stddev), name="fc1")
            dropout1 = tf.layers.dropout(inputs=fc1, rate=self.drop_rate, name="dropout1")
        with tf.name_scope('fc2') as scope:
            fc2 = tf.layers.dense(inputs=dropout1, units=512, activation=tf.nn.elu,
                                  kernel_initializer=tf.truncated_normal_initializer(stddev=FLAGS.stddev), name="fc2")
        return fc2

    def __model(self, input1, input2):
        with tf.variable_scope('out') as scope:
            out1 = self.__core(input1)
            scope.reuse_variables()
            out2 = self.__core(input2)
        with tf.name_scope('distance') as scope:
            distance = tf.sqrt(tf.reduce_sum(tf.square(out1 - out2), axis=1))
        self.distance = distance
        return distance

    def __loss(self):
        with tf.name_scope("error") as scope:
            e_w = self.__model(self.data1, self.data2)
            # e_w_g = self.__model(data1, data2, drop_rate)
            # e_w_i = self.__model(data1, data3, drop_rate)
        with tf.name_scope("loss") as scope:
            Q = tf.constant(FLAGS.Q, dtype=tf.float32, name="Q")
            increase = tf.multiply(tf.multiply(self.label, 2 / Q), tf.square(e_w), name="increase")
            decrease = tf.multiply(tf.multiply(1 - self.label, 2 * Q), tf.exp(-2.77 / Q * e_w), name="decrease")
            loss = tf.reduce_mean(increase + decrease, name="loss")
        return loss

    def __accuracy(self, input1, input2, label):
        # todo add accuracy
        pass

    def train(self, input1, input2, label):
        FLAGS.batch_size = 1
        count = int(np.ceil(len(input1) / FLAGS.batch_size))
        for i in range(count):
            batch_pf = input1[i * FLAGS.batch_size:(i + 1) * FLAGS.batch_size]
            batch_rf = input2[i * FLAGS.batch_size:(i + 1) * FLAGS.batch_size]
            batch_label = label[i * FLAGS.batch_size:(i + 1) * FLAGS.batch_size]
            _, loss, distance, label_ = self.sess.run([self.train_op, self.loss, self.distance, self.label],
                                                      feed_dict={self.data1: batch_pf, self.data2: batch_rf,
                                                                 self.label: np.reshape(batch_label, [-1, 1]),
                                                                 self.drop_rate: np.reshape(FLAGS.train_drop_rate,
                                                                                            [1, 1])})
            print("Loss:{:8.4f},Distance:{:8.4f},Label:{:3.1f}".format(loss, distance[0], label_[0][0]))

    def get_distance(self, input1, input2, label):
        FLAGS.batch_size = 1
        count = int(np.ceil(len(input1) / FLAGS.batch_size))
        trues = []
        falses = []
        for i in range(count):
            batch_pf = input1[i * FLAGS.batch_size:(i + 1) * FLAGS.batch_size]
            batch_rf = input2[i * FLAGS.batch_size:(i + 1) * FLAGS.batch_size]
            batch_label = label[i * FLAGS.batch_size:(i + 1) * FLAGS.batch_size]
            loss, distance, label_ = self.sess.run([self.loss, self.distance, self.label],
                                                   feed_dict={self.data1: batch_pf, self.data2: batch_rf,
                                                              self.label: np.reshape(batch_label, [-1, 1]),
                                                              self.drop_rate: np.reshape(FLAGS.predict_drop_rate,
                                                                                         [1, 1])})
            distance = np.reshape(distance, [-1])
            label_ = np.reshape(label_, [-1])
            assert len(distance) == len(label_)
            for i in range(len(label_)):
                if label_[i] == 1.0:
                    trues.append(distance[i])
                elif label_[i] == 0.0:
                    falses.append(distance[i])
                else:
                    raise Exception("Wrong label")

        figure = plt.figure()
        subplot = figure.add_subplot(111)
        subplot.scatter(x=trues, y=[1] * len(trues), marker="x", color='blue', alpha=0.7, label='true')
        subplot.scatter(x=falses, y=[0] * len(falses), marker="o", color='red', alpha=0.7, label='false')
        mean_true = np.mean(trues)
        mean_false = np.mean(falses)
        subplot.scatter(x=mean_true, y=0.5, marker="x", color='green', alpha=0.7, label='mean_true')
        subplot.scatter(x=mean_false, y=0.5, marker="o", color='black', alpha=0.7,
                        label='mean_false')
        plt.xlabel("distance")
        plt.ylabel("None")
        plt.legend(loc="best")
        now = time.localtime(time.time())
        figure.savefig(
            os.path.join("..", "picture",
                         "{:04d}{:02d}{:02d}_{:02d}{:02d}{:02d}.png".format(now[0], now[1], now[2], now[3], now[4],
                                                                            now[5])))

        return trues, falses, mean_true, mean_false

    def predict(self, input1, input2, label):
        _1, _2, mean_true, mean_false = self.get_distance(input1, input2, label)
        FLAGS.batch_size = 1000
        count = int(np.ceil(len(input1) / FLAGS.batch_size))
        for i in range(count):
            batch_pf = input1[i * FLAGS.batch_size:(i + 1) * FLAGS.batch_size]
            batch_rf = input2[i * FLAGS.batch_size:(i + 1) * FLAGS.batch_size]
            batch_label = label[i * FLAGS.batch_size:(i + 1) * FLAGS.batch_size]
            loss, distance, label_ = self.sess.run([self.loss, self.distance, self.label],
                                                   feed_dict={self.data1: batch_pf, self.data2: batch_rf,
                                                              self.label: np.reshape(batch_label, [-1, 1]),
                                                              self.drop_rate: np.reshape(FLAGS.predict_drop_rate,
                                                                                         [1, 1])})
            f = lambda a, b, x: (b - a) * (a + b - 2 * x)
            batch_predict = np.int32(np.greater_equal(f(mean_true, mean_false, distance), 0))
            accuracy = np.mean(np.equal(batch_label, batch_predict))
            print("Accuracy:{:5.3f}".format(accuracy))


if __name__ == "__main__":
    net = SiameseNet()
    train, test = get_shuffled_data2(num=1)
    net.train(input1=train["product_fingerprint"],
              input2=train["reaction_fingerprint"], label=train["label"])
    net.predict(input1=test["product_fingerprint"],
                input2=test["reaction_fingerprint"], label=test["label"])
