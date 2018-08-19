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
flags.DEFINE_float("learning_rate", 1e-4, "Learning rate of for adam [0.0002]")
flags.DEFINE_integer("batch_size", 1, "MiniBatch size set to be 128")
flags.DEFINE_string("checkpoint_dir", os.path.join("..", "checkpoint"), "dir path to checkpoint")
flags.DEFINE_float("stddev", 1.0, "Stander deviation")
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
            self.optimizer = tf.train.AdamOptimizer(learning_rate=FLAGS.learning_rate)
            self.trainop = self.optimizer.minimize(self.similarity_loss, global_step=self.global_step)
            self.saver = tf.train.Saver(max_to_keep=4)
        else:
            self.load(model_name)
        try:
            self.sess.run(tf.global_variables_initializer())
        except Exception:
            self.sess.run(tf.initialize_all_variables())

    def __build_model(self):
        self.product_fingerprint = tf.placeholder(dtype=tf.float32, shape=[None, 2048], name="product_fingerprint")
        self.reaction_fingerprint = tf.placeholder(dtype=tf.float32, shape=[None, 2048], name="reaction_fingerprint")
        self.label = tf.placeholder(dtype=tf.float32, shape=[None, 1], name="label")
        self.dense_left1 = tf.layers.dense(inputs=self.product_fingerprint,
                                           units=1024,
                                           kernel_initializer=tf.truncated_normal_initializer(stddev=FLAGS.stddev),
                                           activation=tf.nn.sigmoid)
        self.dropout_left1 = tf.layers.dropout(self.dense_left1, rate=self.drop_prob)

        self.dense_left2 = tf.layers.dense(inputs=self.dropout_left1,
                                           units=1024,
                                           kernel_initializer=tf.truncated_normal_initializer(stddev=FLAGS.stddev),
                                           activation=tf.nn.sigmoid)
        self.dropout_left2 = tf.layers.dropout(self.dense_left2, rate=self.drop_prob)

        self.gate_1 = tf.layers.dense(inputs=self.dropout_left2,
                                      units=1024,
                                      kernel_initializer=tf.truncated_normal_initializer(stddev=FLAGS.stddev),
                                      activation=tf.nn.sigmoid)
        self.highway_y_1 = tf.layers.dense(inputs=self.dropout_left2,
                                           units=1024,
                                           kernel_initializer=tf.truncated_normal_initializer(stddev=FLAGS.stddev),
                                           activation=tf.nn.sigmoid)
        self.highway_1 = self.gate_1 * self.highway_y_1 + (1 - self.gate_1) * self.dropout_left2

        self.gate_2 = tf.layers.dense(inputs=self.highway_1,
                                      units=1024,
                                      kernel_initializer=tf.truncated_normal_initializer(stddev=FLAGS.stddev),
                                      activation=tf.nn.sigmoid)
        self.highway_y_2 = tf.layers.dense(inputs=self.highway_1,
                                           units=1024,
                                           kernel_initializer=tf.truncated_normal_initializer(stddev=FLAGS.stddev),
                                           activation=tf.nn.sigmoid)
        self.highway_2 = self.gate_2 * self.highway_y_2 + (1 - self.gate_2) * self.highway_1

        self.gate_3 = tf.layers.dense(inputs=self.highway_2,
                                      units=1024,
                                      kernel_initializer=tf.truncated_normal_initializer(stddev=FLAGS.stddev),
                                      activation=tf.nn.sigmoid)
        self.highway_y_3 = tf.layers.dense(inputs=self.highway_2,
                                           units=1024,
                                           kernel_initializer=tf.truncated_normal_initializer(stddev=FLAGS.stddev),
                                           activation=tf.nn.sigmoid)
        self.highway_3 = self.gate_3 * self.highway_y_3 + (1 - self.gate_3) * self.highway_2

        self.gate_4 = tf.layers.dense(inputs=self.highway_3,
                                      units=1024,
                                      kernel_initializer=tf.truncated_normal_initializer(stddev=FLAGS.stddev),
                                      activation=tf.nn.sigmoid)
        self.highway_y_4 = tf.layers.dense(inputs=self.highway_3,
                                           units=1024,
                                           kernel_initializer=tf.truncated_normal_initializer(stddev=FLAGS.stddev),
                                           activation=tf.nn.sigmoid)
        self.highway_4 = self.gate_4 * self.highway_y_4 + (1 - self.gate_4) * self.highway_3

        self.gate_5 = tf.layers.dense(inputs=self.highway_4,
                                      units=1024,
                                      kernel_initializer=tf.truncated_normal_initializer(stddev=FLAGS.stddev),
                                      activation=tf.nn.sigmoid)
        self.highway_y_5 = tf.layers.dense(inputs=self.highway_4,
                                           units=1024,
                                           kernel_initializer=tf.truncated_normal_initializer(stddev=FLAGS.stddev),
                                           activation=tf.nn.sigmoid)
        self.highway_5 = self.gate_5 * self.highway_y_5 + (1 - self.gate_5) * self.highway_4

        #############################left_layer###############################

        #############################right_layer###############################

        self.dense_right1 = tf.layers.dense(inputs=self.reaction_fingerprint,
                                            units=1024,
                                            kernel_initializer=tf.truncated_normal_initializer(stddev=FLAGS.stddev),
                                            activation=tf.nn.sigmoid)
        self.dropout_right1 = tf.layers.dropout(self.dense_right1, rate=self.drop_prob)

        self.dense_right2 = tf.layers.dense(inputs=self.dropout_right1,
                                            units=1024,
                                            kernel_initializer=tf.truncated_normal_initializer(stddev=FLAGS.stddev),
                                            activation=tf.nn.sigmoid)
        #############################right_layer###############################

        self.model1 = tf.sqrt(tf.reduce_sum(tf.square(self.highway_5), axis=1))
        self.model2 = tf.sqrt(tf.reduce_sum(tf.square(self.dense_right2), axis=1))
        self.multi_sum = tf.reduce_sum(tf.multiply(self.highway_5, self.dense_right2), axis=1)
        self.cosine = tf.divide(tf.abs(self.multi_sum), tf.multiply(self.model1, self.model2))
        # self.similarity = tf.nn.sigmoid(cosine)[:, np.newaxis]
        self.similarity = tf.reshape(self.cosine, shape=[-1, 1])

        return self.similarity

    def __loss(self):
        # self.similarity_loss = tf.reduce_mean(tf.square(
        #     tf.maximum(self.similarity - 0.3, np.zeros_like(self.similarity)) * (1 - self.label) +
        #     tf.minimum(self.similarity - 0.7, np.zeros_like(self.similarity)) * self.label
        # ))
        self.similarity_loss = tf.reduce_sum((1 - self.similarity) * self.label + self.similarity * (1 - self.label))
        # self.similarity_loss = tf.reduce_mean(tf.square(self.similarity - self.label))

    def __accuracy(self):
        # cosine = tf.map_fn(lambda x: 1.0 if tf.greater_equal(x, 0.7) else x, cosine)
        self.label_pred = tf.round(self.similarity)
        self.equals = tf.equal(self.label, self.label_pred)
        self.accuracy = tf.reduce_mean(tf.cast(self.equals, tf.float32))

    def train(self, product_fingerprint, reaction_fingerprint, label):
        self.train_flag = True
        FLAGS.batch_size = 1
        epoch = int(np.ceil(len(product_fingerprint) / FLAGS.batch_size))
        train_all_flag = True
        for i in range(epoch):
            batch_pf = product_fingerprint[i * FLAGS.batch_size:(i + 1) * FLAGS.batch_size]
            batch_rf = reaction_fingerprint[i * FLAGS.batch_size:(i + 1) * FLAGS.batch_size]
            batch_label = label[i * FLAGS.batch_size:(i + 1) * FLAGS.batch_size]
            if train_all_flag:
                _, accuracy, loss, similarity = \
                    self.sess.run([self.trainop, self.accuracy, self.similarity_loss, self.similarity],
                                  feed_dict={self.product_fingerprint: batch_pf,
                                             self.reaction_fingerprint: batch_rf,
                                             self.label: batch_label,
                                             self.drop_prob: np.reshape(0.3, (1, 1))})
            else:
                accuracy, loss, similarity = \
                    self.sess.run([self.accuracy, self.similarity_loss, self.similarity],
                                  feed_dict={self.product_fingerprint: batch_pf,
                                             self.reaction_fingerprint: batch_rf,
                                             self.label: batch_label,
                                             self.drop_prob: np.reshape(0.3, (1, 1))})
                if (similarity[0][0] <= 0.3 and batch_label[0][0] == 0) or \
                        (similarity[0][0] >= 0.7 and batch_label[0][0] == 1):
                    pass
                else:
                    self.sess.run(self.trainop, feed_dict={self.product_fingerprint: batch_pf,
                                                           self.reaction_fingerprint: batch_rf,
                                                           self.label: batch_label,
                                                           self.drop_prob: np.reshape(0.3, (1, 1))})

            print('cosine:', str(similarity[0][0]), "label:", str(batch_label[0][0]), "loss:", loss)
            step = self.sess.run(self.global_step)
            if step > 0 and step % 100 == 0:
                print("Train step {:6d}".format(step))
            if step > 0 and step % 10000 == 0:
                self.save(global_step=step)

    def predict(self, product_fingerprint, reaction_fingerprint, label):
        self.train_flag = False
        FLAGS.batch_size = 128
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
                                         self.drop_prob: np.reshape(0.0, (1, 1))})
            print("Accuracy:{0:10.5f},Total_Loss:{1:10.5f},Average_cosine:{2:10.5f}".format(accuracy, loss,
                                                                                            np.mean(similarity,
                                                                                                    axis=0)[0]))
            # print(str(similarity))

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
    net = SimilarityNet(True)
    pf = np.random.uniform(-1, 1, [1280, 2048]).astype(np.float32)
    rf = np.random.uniform(-1, 1, [1280, 2048]).astype(np.float32)
    label = np.ones([1280, 1], dtype=np.float32)

    net.train(pf, rf, label)
