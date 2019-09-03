# -*- coding: utf-8 -*-
"""
Created on Sat Nov 17 21:35:19 2018

@author: stesn

"""
import matplotlib.pyplot as plt
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

x = np.array([-4.0, -1.0, 1.0, 4.0, 6.0, 9.0])
y = np.array([-14.0, 3.0 , -5.0, -18.0, -15.0, 31.0])




n = x.size
H = np.empty(shape = (n-1,))
for i in range(n-1):
    H[i]=x[i+1] - x[i]
A = np.empty(shape = (n-2,n-2))
for i in range(n-2):
    A[i,i] = (H[i] + H[i+1])/3
for i in range(n-3):
    A[i,i+1] = H[i+1]/6
for i in range(1, n-2):
    A[i,i-1] = H[i]/6
B = np.empty(shape = (n-2,))
for i in range(n-2):
    B[i]= (y[i+2] - y[i+1])/H[i+1] -(y[i+1] - y[i])/H[i]
M = solve(A,B)
M = np.append(M, 0)
M = np.insert(M,0,0)
def S(value):
    i = 0
    for i in range(1,n):
        if((x[i-1] <= value) & (value <= x[i])):
            break
    S = M[i]*((value - x[i-1])**3)/(6*H[i-1])\
    + M[i-1]*((x[i] - value)**3)/(6*H[i-1])\
    + (y[i] - M[i]*(H[i-1]**2)/6)*(value - x[i-1])/H[i-1]\
    + (y[i-1] - M[i-1]*(H[i-1]**2)/6)*(x[i] - value)/H[i-1]
    return S
vfunc = np.vectorize(S)
test = np.empty(shape=[n-1])
for i in range(n-1):
    test[i] = ((x[i+1] + x[i])/2)
ans = vfunc(test)
for i in range(n-1):
    print('x =', test[i], 'y =',ans[i] )

test = np.arange(-4,9,0.02)
plt.figure(figsize=(10,5))
plt.plot(test,vfunc(test),'b',label = 'C-Spline')
plt.plot(x,y, 'go',label = 'Заданные точки')
plt.legend(loc='upper left')
plt.ylim(np.min(vfunc(test))-2,np.max(vfunc(test))+2)
plt.xlim(-5,10)
