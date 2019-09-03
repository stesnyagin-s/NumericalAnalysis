# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 18:10:50 2018

@author: stesn
"""

import numpy as np
np.set_printoptions(suppress = True)

def seidel(B, Beta):
    for i in range(n):
        B[i,:] = B[i,:]/XYZ[i]
        Beta[i] = Beta[i]/XYZ[i]
    X = Beta.copy()
    Beta_norm = np.linalg.norm(Beta, ord = np.inf)
    B_norm = np.linalg.norm(B, ord = np.inf)
    X_old = 0
    print('{:<2s} {:<11} {:<11} {:<11} {:<11} {:<11}'.
          format('i', 'x','y','z', 'eps', 'delta'))
    print('{:<2} {:<11.8f} {:<11.8f} {:<11.8f} {:<11.8f} {:<11.8f}'.
          format(0, X[0], X[1], X[2], 0.0, 0.0))
    for i in range(1,6):
        X_old = X.copy()
        for j in range(n):
            X[j] = B[j,:].dot(X) + Beta[j]
        B_norm = np.linalg.norm(B, ord = np.inf)
        delta = np.linalg.norm(X - X_old, ord = np.inf)
        eps = ((B_norm**(i))/(1 - B_norm))*(Beta_norm)
        print('{:<2} {:<11.8f} {:<11.8f} {:<11.8f} {:<11.8f} {:<11.8f}'.
              format(i, X[0], X[1], X[2], eps, delta))
    return X

def iteration(B, Beta):
    for i in range(n):
        B[i,:] = B[i,:]/XYZ[i]
        Beta[i] = Beta[i]/XYZ[i]
    X = Beta.copy()
    Beta_norm = np.linalg.norm(Beta, ord = np.inf)
    B_norm = np.linalg.norm(B, ord = np.inf)
    X_old = 0
    print('{:<2s} {:<11} {:<11} {:<11} {:<11} {:<11}'.
          format('i', 'x','y','z', 'eps', 'delta'))
    print('{:<2} {:<11.8f} {:<11.8f} {:<11.8f} {:<11.8f} {:<11.8f}'.
          format(0, X[0], X[1], X[2], 0.0, 0.0))
    for i in range(1,11):
        X_old = X.copy()
        X = B.dot(X) + Beta
        delta = np.linalg.norm(X - X_old, ord = np.inf)
        eps = ((B_norm**(i))/(1 - B_norm))*(Beta_norm)
        print('{:<2} {:<11.8f} {:<11.8f} {:<11.8f} {:<11.8f} {:<11.8f}'.
              format(i, X[0], X[1], X[2], eps, delta))
    return X

A = np.array([[-32., -10., -7.],
              [-5., -35., -5.],
              [-1., -5., -38.],])
Beta = np.array([-285., -5., -264.])



n = A.shape[1]
XYZ = np.empty(shape = (n))
B = A.copy()
for i in range(n):
    XYZ[i] = B[i,i]
    B[i,i] = 0
B = -B

print('Метод простых итераций:')
X = iteration(B.copy(), Beta.copy())
print('A*X =')
print(np.matmul(A,X))
print('')
print('Метод Зейделя:')
X = seidel(B.copy(), Beta.copy())
print('A*X =')
print(np.matmul(A,X))

    
