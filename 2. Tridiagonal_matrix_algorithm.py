# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 20:40:14 2018

@author: stesn
"""

import numpy as np
from scipy.sparse import diags
np.set_printoptions(suppress = True)
diagonals = [[-14., 19., 18., 13., -16., 13., 16., 15., 20., 18.],\
             [2., 1., -3., 1., 6., 2., 3., 3., 3.],\
             [-1., 1., 4., 6., 6., 1., 5., -6., 6.]]
B = np.array([-30., 44., 22., -28., -72., 51., -96., -57., -54., 45])


A = diags(diagonals, [0, -1, 1]).todense()
n = A.shape[1]
P = [-A[0,1]/A[0,0]]
Q = [B[0]/A[0,0]]
for i in range(1, n-1):
    P.append(-A[i,i+1]/(A[i,i] + A[i,i-1]*P[i-1]))
for i in range(1, n):
    Q.append((B[i] - A[i,i-1]*Q[i-1])/(A[i,i] + A[i,i-1]*P[i-1]))
P.append(0)
X = np.zeros(shape=(n))
X[n-1] = Q[n-1]
for i in range(n-2,-1,-1):
    X[i] = Q[i] + P[i]*X[i+1]
print('X = ')
print(X)
print('A*X = ')
print(np.matmul(A,X))