
__author__ = "Malvika Viswanathan, Leqi Yin, Yashwant Kurmi, Zhongliang Zu"
__maintainer__ = "Malvika Viswanathan, Zhongliang Zu"
__email__ = "zhongliang.zu@vumc.org"


import numpy as np
import tensorflow  as tf
import scipy.io
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt



index = [1,10,12,13,14,25,38,44,51,54,64]

'''
Please delete following lines of code for X_train and Y before training
'''
X_train = np.zeros((100000,11))
Y = np.zeros((100000,2))
n_features=11
'''
#Load training data from MATLAB file

mat = scipy.io.loadmat(...)
X = mat[...]
Y = mat[...]
R1 = mat[...]

# Scaling Y values to bring them in the same range
Y[:,0] = np.round(Y[:,0]*100, decimals=2)
Y[:,1] = Y[:,1]/100
n_features = len(X[index,0])
X_train = np.transpose(X[index,:]/np.transpose(R1))
'''
def huber_loss(y_true, y_pred, delta=0.1):
    error = y_true - y_pred
    is_small_error = tf.abs(error) <= delta
    small_error_loss = tf.square(error) / 2
    large_error_loss = delta * (tf.abs(error) - (delta / 2))
    return tf.where(is_small_error, small_error_loss, large_error_loss)

# Split the data into training, validation, and testing sets
Xtrain, X_val, Ytrain, Y_val = train_test_split(X_train, Y, test_size=0.1, shuffle=True, random_state=42)


# Define the neural network model
model1 = tf.keras.Sequential([
    tf.keras.layers.Dense(100, input_shape=(n_features,), activation='tanh'),
    tf.keras.layers.Dense(100, activation='tanh'),
    tf.keras.layers.Dense(50, activation='tanh'),
    tf.keras.layers.Dense(2, activation='linear')
])

model1.summary()

# Compile the model with the custom loss function and Adam optimizer
model1.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=0.00001,clipnorm=1.0), loss=huber_loss)

# Train the model
history = model1.fit(Xtrain, Ytrain, epochs=25, batch_size=64, validation_data=(X_val, Y_val))

# Plot training loss and validation loss
plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('model loss')
plt.ylabel('loss')
plt.xlabel('epoch')
plt.legend(['train', 'test'], loc='upper left')
plt.show()