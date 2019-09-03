# -*- coding: utf-8 -*-
"""
Created on Sun Nov 18 17:14:29 2018

@author: stesn
"""
import numpy as np
def solve(A,B):
    n = A.shape[1]
    P = [-A[0,1]/A[0,0]]
    Q = [B[0]/A[0,0]]
    for i in range(1, n-1):
        P.append(-A[i,i+1]/(A[i,i] + A[i,i-1]*P[i-1]))
    for i in range(1, n):
        Q.append((B[i] - A[i,i-1]*Q[i-1])/(A[i,i] + A[i,i-1]*P[i-1]))
    P.append(0)
    X = np.zeros(shape=(n,))
    X[n-1] = Q[n-1]
    for i in range(n-2,-1,-1):
        X[i] = Q[i] + P[i]*X[i+1]
    return X
#y'' - y' +4*y = x^2 - 2x + 4, -5y(1) + 4y'(1)=5, 4y(2) + 5y'(2) = 2
def K(x):
    return 1
def L(x):
    return -1
def M(x):
    return 4
def F(x):
    return x*x  - 2*x + 4
a = 1
b = 2
h = 0.2

R = 4
V = 5
S = -5
W = 4
T = 5
Z = 2
n = int((b-a)/h) + 1
x = np.empty(shape = (n,))
for i in range(n):
    x[i] = a + h*i
    
A = np.zeros(shape = (n,n))
for i in range(1, n-1):
    A[i,i] = M(x[i]) - 2*K(x[i])/(h**2)
for i in range(1, n-1):
    A[i,i-1] = K(x[i])/(h**2) - L(x[i])/(2*h)
for i in range(1, n-1):
    A[i,i+ 1] = K(x[i])/(h**2) + L(x[i])/(2*h)
A[0,0] = S - R/h
A[0,1] = R/h
A[n-1,n-1] = -V/h - W
A[n-1,n-2] = V/h

B = np.zeros(shape = (n,))
for i in range(1,n-1):
    B[i] = F(x[i])
B[0] = T
B[n-1] = -Z
print('A:\n',A)
print('B:\n',B)
Y = solve(A,B)
for i in range(n):
    print('x =',x[i], 'y =', Y[i])