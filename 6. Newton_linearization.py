# -*- coding: utf-8 -*-
"""
Created on Sat Nov 17 12:04:36 2018

@author: stesn
"""
import numpy as np
def z(X):
    return 3*X[0]*X[1]

def df1_x(X):
    return -2

def df1_y(X):
    return 5

def df1_z(X):
    return 4

def df2_x(X):
    return 4*X[0]
    
def df2_y(X):
    return 6*X[1]

def df2_z(X):
    return 3

def df3_x(X):
    return 3*X[1]
    
def df3_y(X):
    return 3*X[0]

def df3_z(X):
    return -1


def f1(X):
    return -2*X[0] + 5*X[1] + 4*X[2]

def f2(X):
    return 2*X[0]*X[0] + 3*X[1]*X[1] + 3*X[2] - 1

def f3(X):
    return 3*X[0]*X[1] - X[2]

def findInverse(A):
    E = np.identity(n)
    AE = np.concatenate((A, E), axis = 1)
    for i in range(n):
        if AE[i,i] == 0:
            k = i + 1
            while k < n:
                if AE[i,k] != 0:
                    AE[[i,k],] = AE[[k,i],]
                k = k + 1
        if AE[i,i] != 0:
            AE[i,:] = AE[i,:]/AE[i,i]
            j = i + 1
            while j < n:
                AE[j,:] = AE[j,:] - AE[i,:]*(AE[j,i]/AE[i,i])
                j = j + 1
    for i in range(n-1,0, -1):
        AE[i,:] = AE[i,:]/AE[i,i]
        j = i - 1
        while(j >= 0):
          AE[j,:] = AE[j,:] - AE[i,:]*(AE[j,i]/AE[i,i])
          j = j - 1
    R = AE[:, range(n, 2*n)]
    return R

def solve(X):
    F = np.array([f1(X), f2(X), f3(X)])
    print(X)
    print('{:<2s} {:<11} {:<11} {:<11} {:<11} {:<11} {:<11}'.
      format('n', 'x_n','y_n','z_n', 'f1(X)', 'f2(X)', 'f3(X)'))
    print('{:<2} {:<11.8f} {:<11.8f} {:<11.8f} {:<11.8f} {:<11.8f} {:<11.8f}'.
      format(0, X[0], X[1], X[2], f1(X), f2(X), f3(X)))
      
    for i in range(1,5):
        J = np.array([[df1_x(X), df1_y(X), df1_z(X)],
                      [df2_x(X), df2_y(X), df2_z(X)],
                      [df3_x(X), df3_y(X), df3_z(X)]])
        J_inv = findInverse(J)
        X = X - J_inv.dot(F)
        F = np.array([f1(X), f2(X), f3(X)])
        print('{:<2} {:<11.8f} {:<11.8f} {:<11.8f} {:<11.8f} {:<11.8f} {:<11.8f}'.
          format(i, X[0], X[1], X[2], f1(X), f2(X), f3(X)))
    
n = 3    
X_1 = np.array([1.0,1.0, 1.0])
X_2 = np.array([-1.0,-1.0, 1.0])
solve(X_1)
print('')
solve(X_2)