import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers, Input
from tensorflow.keras.layers import BatchNormalization, Activation

# Create a sampling layer
class Sampling(layers.Layer):
    """Uses (z_mean, z_log_var) to sample z, the vector encoding a digit."""

    def call(self, inputs):
        z_mean, z_log_var = inputs
        batch = tf.shape(z_mean)[0]
        dim = tf.shape(z_mean)[1]
        epsilon = tf.random.normal(shape=(batch, dim))
        return z_mean + tf.exp(0.5 * z_log_var) * epsilon

# encoder
def vae_encoder(original_dim, latent_dim):
    encoder_inputs = Input(shape=(original_dim,))
    z_mean_dense_linear = layers.Dense(latent_dim, name="z_mean", kernel_initializer='glorot_uniform')(encoder_inputs)
    z_mean_dense_batchnorm = BatchNormalization()(z_mean_dense_linear)
    z_mean_encoded = Activation('relu')(z_mean_dense_batchnorm)

    z_log_var_dense_linear = layers.Dense(latent_dim, name="z_log_var")(encoder_inputs)
    z_log_var_dense_batchnorm = BatchNormalization()(z_log_var_dense_linear)
    z_log_var_encoded = Activation('relu')(z_log_var_dense_batchnorm)

    z = Sampling()([z_mean_encoded, z_log_var_encoded])
    encoder = keras.Model(encoder_inputs, [z_mean_encoded, z_log_var_encoded, z], name="encoder")

    return encoder, z


# decoder
def vae_decoder(drop_rate, z, original_dim):
    drop_layer = layers.Dropout(rate=drop_rate, noise_shape=None)(z)
    decoder = layers.Dense(original_dim, kernel_initializer='glorot_uniform', activation='relu')
    reconstruct = decoder(drop_layer)

    return decoder


# Define the VAE as a Model with a custom train_step
class VAE(keras.Model):
    def __init__(self, encoder, decoder, **kwargs):
        super().__init__(**kwargs)
        self.encoder = encoder
        self.decoder = decoder
        self.total_loss_tracker = keras.metrics.Mean(name="total_loss")
        self.reconstruction_loss_tracker = keras.metrics.Mean(
            name="reconstruction_loss"
        )
        self.kl_loss_tracker = keras.metrics.Mean(name="kl_loss")

    @property
    def metrics(self):
        return [
            self.total_loss_tracker,
            self.reconstruction_loss_tracker,
            self.kl_loss_tracker,
        ]

    def train_step(self, data):
        with tf.GradientTape() as tape:
            z_mean, z_log_var, z = self.encoder(data)
            reconstruction = self.decoder(z)
            reconstruction_loss = tf.reduce_sum(
                keras.losses.binary_crossentropy(data, reconstruction)
            )

            kl_loss = -0.5 * (1 + z_log_var - tf.square(z_mean) - tf.exp(z_log_var))
            kl_loss = tf.reduce_mean(tf.reduce_sum(kl_loss, axis=1))
            total_loss = reconstruction_loss + kl_loss
        grads = tape.gradient(total_loss, self.trainable_weights)
        self.optimizer.apply_gradients(zip(grads, self.trainable_weights))
        self.total_loss_tracker.update_state(total_loss)
        self.reconstruction_loss_tracker.update_state(reconstruction_loss)
        self.kl_loss_tracker.update_state(kl_loss)
        return {
            "loss": self.total_loss_tracker.result(),
            "reconstruction_loss": self.reconstruction_loss_tracker.result(),
            "kl_loss": self.kl_loss_tracker.result(),
        }


def bulid_vae(original_dim, latent_dim, drop_rate, learning_rate):
    encoder, z = vae_encoder(original_dim, latent_dim)
    decoder = vae_decoder(drop_rate, z, original_dim)
    vae = VAE(encoder, decoder)
    vae.compile(optimizer=keras.optimizers.Adam(learning_rate=learning_rate))

    return vae


def run_vae(input, vae, epochs, batch_size):
    vae.fit(input, epochs=epochs, batch_size=batch_size)


def get_latent(input, vae):
    _, _, z_latent = vae.encoder(input)
    return z_latent