# -*- coding: utf-8 -*-
"""HW2.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1mO1005oYOyQtD_1gB6AXhvb0sCKKKC1a
"""

import numpy as np
import pandas as pd
from sklearn import linear_model
from sklearn import preprocessing
import copy

class Algorithm:
	
  def __init__(self, lambduh=0.1, max_iter=1000):
    self.lambduh = lambduh
    self.max_iter = max_iter
    
  def soft_threshold(self, a, lambduh):
    if a < -lambduh:
      return a+lambduh
    elif a > lambduh:
        return a-lambduh
    else:
      return 0
				   
  def min_beta_multivariate(self, x, y, beta, j):
    n = len(y)
    selector = [i for i in range(x.shape[1]) if i != j]
    norm_x_j = np.linalg.norm(x[:, j])
    a = x[:, j].dot(y[:, np.newaxis] - x[:, selector].dot(beta[:, np.newaxis][selector, :]))
    passin = self.lambduh*n/2
    res = self.soft_threshold(a, passin)
    return res/(norm_x_j**2)

  def generate_simulate_data(self):
    n = 100
    np.random.seed(123456789)
    X = np.random.normal(loc=2, scale=1, size=n)
    epsilon = np.random.normal(loc=0, scale=0.1, size=n)
    y = 1 + 2*X + 3*X**2 + 4*X**3 + epsilon
    predictors = np.vstack([X**i for i in range(1, 11)]).T
    scaler = preprocessing.StandardScaler()
    x = scaler.fit_transform(predictors)
    y = y-np.mean(y)
    return x, y

  def cycliccoorddescent(self, x, y, beta_init):
    """
		cycliccoorddescent that implements the cyclic coordinate descent algorithm. The cyclic 
		coordinate descent algorithm proceeds sequentially. At each iteration, the algorithm 
		increments the index j of the coordinate to minimize over. Then the algorithm performs 
		partial minimization with respect to the coordinate beta_j corresponding to that index. 
		After updating the coordinate beta_j , the algorithm proceeds to the next iteration. 
		The function takes as input the initial point, the initial step-size value, and the 
		maximum number of iterations. The stopping criterion is the maximum number of iterations.
		"""
    beta = copy.deepcopy(beta_init)
    beta_vals = beta
    d = np.size(x, 1)
    iter = 0
    while iter < self.max_iter:        
      for j in range(d):
        min_beta_j = self.min_beta_multivariate(x, y, beta, j)
        beta[j] = min_beta_j
      beta_vals = np.vstack((beta_vals, beta))
      iter += 1
      if iter % 100 == 0:
        print('Coordinate descent iteration', iter)
    return beta_vals

if __name__=='__main__':
    algorithm = Algorithm()
    x, y = algorithm.generate_simulate_data()
    beta_init = np.zeros(np.size(x, 1))
    betas_cyclic = algorithm.cycliccoorddescent(x, y, beta_init)

print(betas_cyclic[1000])

### comparing with sklearn.linear.model.Lasso
"""Generate simulated data sets and standarize data."""
n = 100
np.random.seed(123456789)
X = np.random.normal(loc=2, scale=1, size=n)
epsilon = np.random.normal(loc=0, scale=0.1, size=n)
y = 1 + 2*X + 3*X**2 + 4*X**3 + epsilon
predictors = np.vstack([X**i for i in range(1, 11)]).T
scaler = preprocessing.StandardScaler()
x = scaler.fit_transform(predictors)
y = y-np.mean(y)
""" fit model """
lasso_model = linear_model.Lasso(alpha=0.05, fit_intercept=False, max_iter=1000,selection='cyclic',normalize=False)
lasso_fit = lasso_model.fit(x, y)
beta_star = lasso_fit.coef_

print(beta_star)