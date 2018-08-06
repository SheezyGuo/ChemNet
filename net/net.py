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
import numpy as np
import os

flags = tf.app.flags

# flags.DEFINE_integer("epoch", 25, "Epoch to train [25]")
flags.DEFINE_float("learning_rate", 0.0002, "Learning rate of for adam [0.0002]")
flags.DEFINE_integer("batch_size", 128, "MiniBatch size set to be 128")
flags.DEFINE_string("checkpoint_dir", os.path.join("..", "checkpoint"), "dir path to checkpoint")
FLAGS = flags.FLAGS


class SimilarityNet(object):
    def __init__(self, train_flag=False, model_name=None):
        self.train_flag = train_flag
        self.sess = tf.Session()
        self.keep_prob = tf.placeholder(dtype=tf.float32, shape=[1, 1], name="keep_prob")
        self.__build_model()
        self.__loss()
        self.__accuracy()
        if train_flag:
            self.global_step = tf.Variable(0, name='global_step', trainable=False)
            self.optimizer = tf.train.AdamOptimizer(learning_rate=FLAGS.learning_rate)
            self.trainop = self.optimizer.minimize(self.similarity_loss, global_step=self.global_step)
            self.saver = tf.train.Saver(max_to_keep=4)
        else:
            self.load(model_name)
        try:
            self.sess.run(tf.global_variables_initializer())
        except:
            self.sess.run(tf.initialize_all_variables())

    def __build_model(self):
        self.product_fingerprint = tf.placeholder(dtype=tf.float32, shape=[None, 2048], name="product_fingerprint")
        self.reaction_fingerprint = tf.placeholder(dtype=tf.float32, shape=[None, 2048], name="reaction_fingerprint")
        self.label = tf.placeholder(dtype=tf.float32, shape=[None, 1], name="label")
        self.dense_left = tf.layers.dense(inputs=self.product_fingerprint,
                                          units=1024,
                                          kernel_initializer=tf.truncated_normal_initializer(stddev=0.01),
                                          activation=tf.nn.elu)
        self.dropout = tf.layers.dropout(self.dense_left, rate=self.keep_prob)
        self.dense_left = self.dropout

        self.gate_1 = tf.layers.dense(inputs=self.dense_left,
                                      units=1024,
                                      kernel_initializer=tf.truncated_normal_initializer(stddev=0.01),
                                      activation=tf.nn.sigmoid)
        self.highway_y_1 = tf.layers.dense(inputs=self.dense_left,
                                           units=1024,
                                           kernel_initializer=tf.truncated_normal_initializer(stddev=0.01),
                                           activation=tf.nn.elu)
        self.highway_1 = self.gate_1 * self.highway_y_1 + (1 - self.gate_1) * self.dense_left

        self.gate_2 = tf.layers.dense(inputs=self.highway_1,
                                      units=1024,
                                      kernel_initializer=tf.truncated_normal_initializer(stddev=0.01),
                                      activation=tf.nn.sigmoid)
        self.highway_y_2 = tf.layers.dense(inputs=self.highway_1,
                                           units=1024,
                                           kernel_initializer=tf.truncated_normal_initializer(stddev=0.01),
                                           activation=tf.nn.elu)
        self.highway_2 = self.gate_2 * self.highway_y_2 + (1 - self.gate_2) * self.highway_1

        self.gate_3 = tf.layers.dense(inputs=self.highway_2,
                                      units=1024,
                                      kernel_initializer=tf.truncated_normal_initializer(stddev=0.01),
                                      activation=tf.nn.sigmoid)
        self.highway_y_3 = tf.layers.dense(inputs=self.highway_2,
                                           units=1024,
                                           kernel_initializer=tf.truncated_normal_initializer(stddev=0.01),
                                           activation=tf.nn.elu)
        self.highway_3 = self.gate_3 * self.highway_y_3 + (1 - self.gate_3) * self.highway_2

        self.gate_4 = tf.layers.dense(inputs=self.highway_3,
                                      units=1024,
                                      kernel_initializer=tf.truncated_normal_initializer(stddev=0.01),
                                      activation=tf.nn.sigmoid)
        self.highway_y_4 = tf.layers.dense(inputs=self.highway_3,
                                           units=1024,
                                           kernel_initializer=tf.truncated_normal_initializer(stddev=0.01),
                                           activation=tf.nn.elu)
        self.highway_4 = self.gate_4 * self.highway_y_4 + (1 - self.gate_4) * self.highway_3

        self.gate_5 = tf.layers.dense(inputs=self.highway_4,
                                      units=1024,
                                      kernel_initializer=tf.truncated_normal_initializer(stddev=0.01),
                                      activation=tf.nn.sigmoid)
        self.highway_y_5 = tf.layers.dense(inputs=self.highway_4,
                                           units=1024,
                                           kernel_initializer=tf.truncated_normal_initializer(stddev=0.01),
                                           activation=tf.nn.elu)
        self.highway_5 = self.gate_5 * self.highway_y_5 + (1 - self.gate_5) * self.highway_4

        #############################left_layer###############################

        #############################right_layer###############################

        self.dense_right = tf.layers.dense(inputs=self.reaction_fingerprint,
                                           units=1024,
                                           kernel_initializer=tf.truncated_normal_initializer(stddev=0.01),
                                           activation=tf.nn.elu)

        #############################right_layer###############################

        model1 = tf.sqrt(tf.reduce_sum(tf.square(self.highway_5), axis=1))
        model2 = tf.sqrt(tf.reduce_sum(tf.square(self.dense_right), axis=1))
        multi_sum = tf.reduce_sum(tf.multiply(self.highway_5, self.dense_right), axis=1)
        cosine = multi_sum / (model1 * model2)

        self.similarity = tf.nn.sigmoid(cosine)[:, np.newaxis]

        return self.similarity

    def __loss(self):
        # self.similarity_loss = tf.nn.sigmoid_cross_entropy_with_logits(labels=self.y, logits=self.similarity)
        self.similarity_loss = tf.reduce_sum(tf.square(self.label - self.similarity))

    def __accuracy(self):
        self.equals = tf.equal(self.label, tf.round(self.similarity))
        self.accuracy = tf.reduce_mean(tf.cast(self.equals, tf.float32))

    def train(self, product_fingerprint, reaction_fingerprint, label):
        self.train_flag = True
        epoch = int(np.ceil(len(product_fingerprint) / FLAGS.batch_size))
        for i in range(epoch):
            batch_pf = product_fingerprint[i * FLAGS.batch_size:(i + 1) * FLAGS.batch_size]
            batch_rf = reaction_fingerprint[i * FLAGS.batch_size:(i + 1) * FLAGS.batch_size]
            batch_label = label[i * FLAGS.batch_size:(i + 1) * FLAGS.batch_size]
            accuracy, _, loss = \
                self.sess.run([self.accuracy, self.trainop, self.similarity_loss],
                              feed_dict={self.product_fingerprint: batch_pf,
                                         self.reaction_fingerprint: batch_rf,
                                         self.label: batch_label,
                                         self.keep_prob: np.array(0.7).reshape([1, 1])})
            step = self.sess.run(self.global_step)
            if step % 10 == 0:
                self.save(global_step=step)
            print("#", step, accuracy, loss)

    def predict(self, product_fingerprint, reaction_fingerprint, label):
        self.train_flag = False
        epoch = int(np.ceil(len(product_fingerprint) / FLAGS.batch_size))
        for i in range(epoch):
            batch_pf = product_fingerprint[i * FLAGS.batch_size:(i + 1) * FLAGS.batch_size]
            batch_rf = reaction_fingerprint[i * FLAGS.batch_size:(i + 1) * FLAGS.batch_size]
            batch_label = label[i * FLAGS.batch_size:(i + 1) * FLAGS.batch_size]
            accuracy, loss, similarity = \
                self.sess.run([self.accuracy, self.similarity_loss, self.similarity],
                              feed_dict={self.product_fingerprint: batch_pf,
                                         self.reaction_fingerprint: batch_rf,
                                         self.label: batch_label,
                                         self.keep_prob: np.array(1.0).reshape([1, 1])})
            print("############", accuracy, loss, np.mean(similarity, axis=0)[0])

    def save(self, global_step):
        self.saver.save(self.sess, os.path.join(FLAGS.checkpoint_dir, "similarity_net_model"), global_step=global_step)

    def load(self, file_name):
        self.saver = tf.train.import_meta_graph(os.path.join(FLAGS.checkpoint_dir, file_name))
        self.saver.restore(self.sess, tf.train.latest_checkpoint(FLAGS.checkpoint_dir))


if __name__ == "__main__":
    net = SimilarityNet(True)
    pf = np.random.uniform(-1, 1, [1280, 2048]).astype(np.float32)
    rf = np.random.uniform(-1, 1, [1280, 2048]).astype(np.float32)
    label = np.ones([1280, 1], dtype=np.float32)

    net.train(pf, rf, label)
