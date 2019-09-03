# -*- coding: utf-8 -*-
"""
Created on Sat Nov 17 13:55:53 2018

@author: stesn
"""
import matplotlib.pyplot as plt

import numpy as np
def Lagrange(x,X,Y):
    n = X.size
    L = 0
    for i in range(n):
        p = 1
        for k in range(n):
            if(k != i):
                p*=(x - X[k])/(X[i] - X[k])
        L+=Y[i]*p
    return L

def Newton(x,X,Y):
    F = []
    def findDif(i,k):
        if(i!=k):
            dif = (findDif(i+1, k) - findDif(i, k-1))/(X[k] - X[i])
            if(i == 0):
                F.append(dif)
            return dif
        else:
            if(i == 0):
                F.append(Y[i])
            return Y[i]
    n = X.size
    findDif(0,n-1)
    P = F[0]
    m = 1
    for i in range(n-1):
        m *= x - X[i]
        P = P + m*F[i+1]
    return P

def Canon(X,Y):
    n = X.size
    A = np.zeros(shape= (n,n))
    B = np.zeros(shape = (n,1))
    for i in range(n):
        for j in range(n):
            A[j,i] = X[j]**(n-i-1)
        B[i] = Y[i]
    AB = np.concatenate((A, B), axis = 1)
    for i in range(n):
        if AB[i,i] == 0:
            k = i + 1
            while k < n:
                if AB[i,k] != 0:
                    AB[[i,k],] = AB[[k,i],]
                k = k + 1
        if AB[i,i] != 0:
            AB[i,:] = AB[i,:]/AB[i,i]
            j = i + 1
            while j < n:
                AB[j,:] = AB[j,:] - AB[i,:]*(AB[j,i]/AB[i,i])
                j = j + 1
    for i in range((n-1),-1, -1): 
        for j in range(i+1,n):
            AB[i,n] = AB[i,n] - AB[i, j]*AB[j,n]
    return AB[:, n]

x_1 = np.array([0.4, 1.8, 3.2, 5.1, 6.8])
y_1 = np.array([21.0, 21.7, 23.2, 21.1, 18.4])
x_2 = np.array([1.8, 3.2, 5.1, 6.8, 8.5])
y_2 = np.array([21.7, 23.2, 21.1, 18.4, 16.7])
x_3 = np.array([0.4, 1.8, 3.2, 5.1, 6.8, 8.5])
y_3 = np.array([21.0, 21.7, 23.2, 21.1, 18.4, 16.7])
x = 2.5

print("x =",x)
print('Многочлен Лагранжа:', Lagrange(x,x_1,y_1))
print('Многочлен Ньютона:', Newton(x,x_2,y_2))
coefs = Canon(x_3,y_3)
def canonf(x):
    P = 0
    for i in range(coefs.size):
        P = P + coefs[i]*((x**(coefs.size - i - 1)))
    return P
print("Канонический многочлен:", canonf(x))

x = np.arange(-5,12,0.1)
plt.figure(figsize=(10,5))
plt.plot(x,Lagrange(x,x_1,y_1),'y',label = 'Многочлен Лагранжа')
plt.plot(x,Newton(x,x_2,y_2),'b',label = 'Многочлен Ньютона')
plt.plot(x, canonf(x),'r',label = 'Канонический многочлен')
plt.plot(x_3,y_3, 'go',label = 'Заданные точки')
plt.legend(loc='upper right')
plt.ylim(15,25)
plt.xlim(-1,11)
