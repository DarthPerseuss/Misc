# Runs a Boltzmann Machine to classify a random set of 
# uniformly distributed discrete numbers in {-1, 1}. 

import numpy as np
from numpy import random
import matplotlib.pyplot as plt

n_neurons = 10
t_steps = 500
learning_steps = 200
lr = 1

P = 500
init_weights = np.triu(random.uniform(-1, 1, [10, 10]))
W = (init_weights + init_weights.T) - np.diag(np.diag(init_weights))
theta = np.zeros([1, n_neurons])

states = np.zeros([10, t_steps])

data = np.random.randint(-1, 1, [n_neurons, P])


clamped_stat1 = np.zeros([10, 1])
clamped_stat2 = np.zeros([10, 10])

for j in range(10):
    clamped_stat1[j, 0] = np.sum(data[j, :]) / P
    for i in range(10):
        clamped_stat2[i, j] = np.sum(data[j, :] * data[i, :]) / P


print(clamped_stat1, clamped_stat2)
delta_w = []

x = random.randint(-1, 1, n_neurons)

for ii in range(learning_steps):
    for jj in range(t_steps):
        n_choose = random.choice(np.arange(10))
        x[n_choose] = np.tanh(np.sum(W * x) + theta[0, n_choose])
        ref = random.uniform()
        if ref > x[n_choose]:
            x[n_choose] = -1
        else:
            x[n_choose] = 1
        states[:, jj] = x

    free_stat1 = np.zeros([10, 1])
    free_stat2 = np.zeros([10, 10])
    for col in range(10):
        free_stat1[col, 0] = np.sum(states[col, :]) / t_steps
        for row in range(10):
            free_stat2[row, col] = np.sum(np.transpose(states[col, :] * states[row, :])) / t_steps

    weights = W

    for iii in range(10):
        theta[0, iii] = theta[0, iii] + lr * (clamped_stat1[iii, 0] - free_stat1[iii, 0])
        for jjj in range(10):
            print(weights[iii, jjj])
            if jjj != iii:
                weights[iii, jjj] = weights[iii, jjj] + lr * (clamped_stat2[iii, jjj] - free_stat2[iii, jjj])
            print(weights[iii, jjj])
    delta_w.append(np.sum(np.sum(W - weights)))

print(delta_w)
plt.figure(1)
plt.plot(delta_w, 'r-')
plt.xlabel('Iteration')
plt.ylabel('Change in weights')
plt.title('Convergence')
plt.show()
