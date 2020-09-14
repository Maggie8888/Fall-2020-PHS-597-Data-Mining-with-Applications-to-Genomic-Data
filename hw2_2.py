# -*- coding: utf-8 -*-
"""HW2_2.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1Tqr-3G0tu0rzu20sVkhh7OfEIxt3_jBs
"""

import numpy as np
import pandas as pd
from sklearn.cross_decomposition import PLSRegression
from sklearn import preprocessing
np.random.seed(123456789)

### define function for PLS ###
class Algorithm:

  def generate_simulate_data(self):
    n = 1000
    q = 3
    p = 10
    X = np.random.normal(size=n * p).reshape((n, p))
    B = np.array([[1, 2] + [0] * (p - 2)] * q).T
    Y = np.dot(X, B) + np.random.normal(size=n * q).reshape((n, q)) + 5
    scaler = preprocessing.StandardScaler()
    x = scaler.fit_transform(X)
    y = scaler.fit_transform(Y)
    return x, y


  def __init__(self,tol1=0.001, tol2=1):
    self.tol1 = 0.001
    self.tol2 = 1

  def colMaxSum(self,mat): 
    # Variable to store index of column 
    # with maximum
    idx = -1
    # Variable to store max Sum
    maxSum = -10**9
    # Traverse matrix column wise 
    for i in range(np.size(mat,1)): 
        Sum = 0
        # calculate Sum of column 
        for j in range(len(mat)): 
            Sum += mat[j][i]**2 
        # Update maxSum if it is less  
        # than current Sum 
        if (Sum > maxSum): 
          maxSum = Sum
          # store index 
          idx = i 
        # return result 
    return idx
  
 
  def pls(self, x, y):
    """
##PLS   Partial Least Squares Regrassion
##between the independent variables, X and dependent Y as
## X = T*P' + E;
## Y = U*Q' + F = T*B*Q' + F1;
##
## Inputs:
## X     data matrix of independent variables
## Y     data matrix of dependent variables
## tol1   the tolerant of convergence 
## tol2   the tolerant of convergence 
###
## Outputs:
## T     score matrix of X
## P     loading matrix of X
## U     score matrix of Y
## Q     loading matrix of Y
## B     matrix of regression coefficient
## W     weight matrix of X
## Using the PLS model, for new X1, Y1 can be predicted as
## Y1 = (X1*P)*B*Q' = X1*(P*B*Q')
## or
## Y1 = X1*(W*inv(P'*W)*inv(T'*T)*T'*Y)
##
## Without Y provided, the function will return the principal components as
## X = T*P' + E
##
		"""
    rx=len(x)
    cx=np.size(x,1)
    ry=len(y)
    cy=np.size(y,1)
    n=max(cx,cy)
    T=np.zeros((rx,n))
    P=np.zeros((cx,n))
    U=np.zeros((ry,n))
    Q=np.zeros((cy,n))
    B=np.zeros((n,n))
    W=P
    k=0
    """ iteration loop if residual is larger than specfied"""
    while np.linalg.norm(y,2)>self.tol2 and k<n and np.isnan(np.linalg.norm(y,2))==False:
      """ choose the column of x has the largest square of sum as t."""
      """ choose the column of y has the largest square of sum as u."""   
      tidx =  self.colMaxSum(x)
      uidx =  self.colMaxSum(y)
      t1 = x[:,tidx]
      u = y[:,uidx]
      t = np.zeros((rx,1))
      """ iteration for outer modeling until convergence"""

      while np.linalg.norm((t1-t),2) > self.tol1:
          w = np.dot(np.transpose(x),u)
          w = w/np.linalg.norm(w,2)
          t = t1
          t1 = np.dot(x,w)
          q = np.dot(np.transpose(y),t1)
          q = q/np.linalg.norm(q,2)
          u = np.dot(y,q)
          """ update p based on t """
      t=t1
      p=np.dot(np.transpose(x),t[:, np.newaxis])/(np.transpose(t)*t)
      pnorm=np.linalg.norm(p,2)
      p=p/pnorm
      t=pnorm*t
      w=pnorm*w
      """ regression and residuals """
      b = np.dot(np.transpose(u),t)/np.dot(np.transpose(t),t)
      x = x- np.dot(t,np.transpose(p))
      #y = y - b*np.dot(t,np.transpose(q))
      y = y - np.dot(t[:, np.newaxis],np.transpose(q[:, np.newaxis]))
      """ save iteration results to outputs: """
      #k=k+1
      T=pd.DataFrame(T)
      P=pd.DataFrame(P)
      U=pd.DataFrame(U)
      Q=pd.DataFrame(Q)
      W=pd.DataFrame(W)
      B=pd.DataFrame(B)
      T.iloc[:,k]=t
      P.iloc[:,k]=p
      U.iloc[:,k]=u
      Q.iloc[:,k]=q
      W.iloc[:,k]=w
      B.iloc[k,k]=b
      k=k+1
    T.iloc[:,k:]=[]
    P.iloc[:,k:]=[]
    U.iloc[:,k:]=[]
    Q.iloc[:,k:]=[]
    W.iloc[:,k:]=[]
    B=B.iloc[0:k,0:k]
    return T,P,U,Q,W,B

if __name__=='__main__':
    algorithm = Algorithm()
    x, y = algorithm.generate_simulate_data()
    T,P,U,Q,W,B = algorithm.pls(x, y)

print(B)
print(np.diagonal(B))
print(np.dot(B,np.transpose(Q)))

np.random.seed(123456789)
n = 1000
q = 3
p = 10
X = np.random.normal(size=n * p).reshape((n, p))
B = np.array([[1, 2] + [0] * (p - 2)] * q).T
Y = np.dot(X, B) + np.random.normal(size=n * q).reshape((n, q)) + 5
scaler = preprocessing.StandardScaler()
x = scaler.fit_transform(X)
y =scaler.fit_transform(Y)
pls2 = PLSRegression(n_components=3)
fit2= pls2.fit(x, y)
#print "Estimated B"
print(fit2.coef_)