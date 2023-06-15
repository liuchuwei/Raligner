import argparse
import os
import sys

# args define
parser = argparse.ArgumentParser(description='VAE for dimension reduction')
parser.add_argument('-i', '--input', required=True, help="intput file")
parser.add_argument('-o', '--output', required=True, help="output directory")
args = parser.parse_args(sys.argv[1:])
global FLAGS
FLAGS = args

import pandas as pd
from dataLoader import dataLoad
import tensorflow as tf
from vae_model import bulid_vae, run_vae, get_latent

# creat output directory
if not os.path.isdir(FLAGS.output):
    os.mkdir(FLAGS.output)


# read and preprocess data
exp, gene = dataLoad(FLAGS.input)
original_dim = len(gene)
input_exp = tf.constant(exp, dtype = tf.float32)

for i in range(100):
    # bulid model
    latent_dim = 100
    drop_rate = 0.2
    batch_size = 100
    epochs = 100
    learning_rate = 0.0005

    model = bulid_vae(original_dim=original_dim, latent_dim = latent_dim, drop_rate=drop_rate, learning_rate=learning_rate)

    # run model
    run_vae(input=input_exp, vae=model, epochs=epochs, batch_size=batch_size)

    # get latent
    z_latent = get_latent(input=input_exp, vae=model)

    # output
    out_path = FLAGS.output + "/latent_" + str(i) + ".csv"
    z_latent = pd.DataFrame(z_latent)
    z_latent.to_csv(out_path, index=False)