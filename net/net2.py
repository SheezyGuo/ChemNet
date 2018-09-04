# -*- coding: utf-8 -*-
"""
__title__ = ''
__author__ = 'duolaoa'
__date_ = '2018.08.25'
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
import os
from convert_data import *

flags = tf.app.flags
flags.DEFINE_float("learning_rate", 0.1, "Learning rate of for adam [0.0002]")
flags.DEFINE_integer("batch_size", 1, "MiniBatch size set to be 128")
flags.DEFINE_string("checkpoint_dir", os.path.join("..", "checkpoint2"), "dir path to checkpoint")
flags.DEFINE_float("stddev", 0.01, "Stander deviation")
FLAGS = flags.FLAGS


class SimilarityNet(object):
    def __init__(self, train_flag=False, model_name=None):
        self.train_flag = train_flag
        self.sess = tf.Session()
        self.drop_prob = tf.placeholder(dtype=tf.float32, shape=[1, 1], name="drop_prob")
        self.__build_model()
        self.__loss()
        self.__accuracy()
        if train_flag:
            self.global_step = tf.Variable(0, name='global_step', trainable=False)
            self.optimizer = tf.train.GradientDescentOptimizer(learning_rate=FLAGS.learning_rate)
            self.train_op = self.optimizer.minimize(self.loss, global_step=self.global_step)
            self.saver = tf.train.Saver(max_to_keep=4)
        else:
            self.load(model_name)
        try:
            self.sess.run(tf.global_variables_initializer())
        except Exception:
            self.sess.run(tf.initialize_all_variables())

    def __build_model(self):
        self.spliced_fingerprint = tf.placeholder(dtype=tf.float32, shape=[None, 4096], name="data")
        self.label = tf.placeholder(dtype=tf.float32, shape=[None], name="label")
        self.dense1 = tf.layers.dense(inputs=self.spliced_fingerprint, units=4096, activation=tf.nn.relu,
                                      kernel_initializer=tf.truncated_normal_initializer(stddev=FLAGS.stddev))
        self.dense1 = tf.layers.batch_normalization(self.dense1, momentum=0.8)
        self.dense2 = tf.layers.dense(inputs=self.dense1, units=4096, activation=tf.nn.relu,
                                      kernel_initializer=tf.truncated_normal_initializer(stddev=FLAGS.stddev))
        self.dense2 = tf.layers.batch_normalization(self.dense2, momentum=0.8)
        self.dense3 = tf.layers.dense(inputs=self.dense2, units=1, activation=tf.nn.sigmoid,
                                      kernel_initializer=tf.truncated_normal_initializer(stddev=FLAGS.stddev))

    def __loss(self):
        # self.loss = tf.reduce_sum(tf.nn.sparse_softmax_cross_entropy_with_logits(labels=tf.cast(self.label, tf.int32),
        #                                                                            logits=self.dense3, name="loss"))
        self.loss = tf.reduce_mean(tf.square(self.dense3 - self.label))

    def __accuracy(self):
        # self.label_pred = tf.argmax(self.dense3, axis=1)
        self.label_pred = tf.round(self.dense3)
        self.equals = tf.equal(tf.reshape(tf.cast(self.label, tf.float32), [-1, 1]), self.label_pred)
        self.accuracy = tf.reduce_mean(tf.cast(self.equals, tf.float32))

    def train(self, spliced_fingerprint, label):
        FLAGS.batch_size = 1
        count = int(np.ceil(len(spliced_fingerprint) / FLAGS.batch_size))
        for i in range(count):
            batch_sf = spliced_fingerprint[i * FLAGS.batch_size:(i + 1) * FLAGS.batch_size]
            batch_label = label[i * FLAGS.batch_size:(i + 1) * FLAGS.batch_size]
            a = np.array(batch_sf)
            _, loss, accuracy, dense3, label_pred = self.sess.run(
                [self.train_op, self.loss, self.accuracy, self.dense3, self.label_pred],
                feed_dict={self.spliced_fingerprint: batch_sf,
                           self.label: batch_label})
            # print(
            #     "Loss:{:8.4f},Accuracy:{:8.4f},class1={:8.4f},class2={:8.4f}".format(loss, accuracy, dense3[0][0],
            #                                                                          dense3[0][1]))
            print(
                "Loss:{:8.4f},Accuracy:{:8.4f},dense={:8.4f}".format(loss, accuracy, dense3[0][0]))

            # print("y_:", label_pred)

    def predict(self, spliced_fingerprint, label):
        FLAGS.batch_size = 5000
        count = int(np.ceil(len(spliced_fingerprint) / FLAGS.batch_size))
        for i in range(count):
            batch_sf = spliced_fingerprint[i * FLAGS.batch_size:(i + 1) * FLAGS.batch_size]
            batch_label = label[i * FLAGS.batch_size:(i + 1) * FLAGS.batch_size]
            accuracy = self.sess.run(self.accuracy,
                                     feed_dict={self.spliced_fingerprint: batch_sf, self.label: batch_label})
            print("Accuracy:{:8.4f}".format(accuracy))

    def save(self, global_step):
        print("Saving model of step {}...".format(global_step))
        self.saver.save(self.sess, os.path.join(FLAGS.checkpoint_dir, "similarity_net_model"), global_step=global_step)
        print("Saved.".format(global_step))

    def load(self, file_name):
        print("Loading model from {}...".format(file_name))
        self.saver = tf.train.import_meta_graph(os.path.join(FLAGS.checkpoint_dir, file_name))
        self.saver.restore(self.sess, tf.train.latest_checkpoint(FLAGS.checkpoint_dir))
        print("Completed")


if __name__ == "__main__":
    train, test = get_shuffled_data2()
    net = SimilarityNet(True)
    net.train(spliced_fingerprint=train["spliced_fingerprint"], label=train["label"])
    net.predict(spliced_fingerprint=test["spliced_fingerprint"], label=test["label"])
