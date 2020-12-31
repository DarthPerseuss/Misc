# Implementation of a simple single-layered perceptron

import numpy as np

n_dim = 50
n_inputs = 120  # this is 'p' in the exercise
lr = 0.1

# ins_big = ins[np.sum(ins, 2) >= 25]
# ins_small = ins[np.sum(ins, 2) < 25]


# Perceptron
indices = np.arange(n_inputs)


def percept(inputs, weights, iterations):

    # print('dot is {}'.format(np.dot(weights, inputs.T)))
    # print('sum is {}'.format(np.sum(inputs)))
    # print('mean is {}'.format(np.mean(inputs)))
    for i in range(iterations):
        ind = np.random.choice(indices, replace=True)
        if np.sum(inputs[ind]) >= n_dim/2 > np.sum(weights[ind] * inputs[ind]):
            weights[ind] += lr * np.mean(inputs[ind])

        if np.sum(inputs[ind]) < n_dim/2 <= np.sum(weights[ind] * inputs[ind]):
            weights[ind] -= lr * np.mean(inputs[ind])

    return weights


def predict(inputs, weights):
    predictions = []
    for ii in indices:
        if np.sum(inputs[ii]) >= n_dim / 2:
            if np.sum(weights[ii] * inputs[ii]) >= n_dim / 2:
                predictions.append(1)
            else:
                predictions.append(-1)
        else:
            if np.sum(weights[ii] * inputs[ii]) < n_dim / 2:
                predictions.append(1)
            else:
                predictions.append(-1)

    return predictions


# w = np.random.uniform(0.5, 1.5, n_dim)
w = np.zeros(n_dim)
ins = np.random.randint(0, 2, (n_inputs, n_dim))


updated_weights = percept(ins, w, 10000)
print(updated_weights)

results = predict(ins, updated_weights)

print(results)
