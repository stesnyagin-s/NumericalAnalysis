 # -*- coding: utf-8 -*-
"""
Created on Sat Nov 17 18:22:13 2018

@author: stesn
"""
import numpy as np
import matplotlib.pyplot as plt

def createMatrix(X,Y, k):
    A = np.empty(shape = (k+1,k+1))
    B = np.empty(shape = (k+1,1))
    for i in range(k+1):
        for j in range(k+1):
            A[i,j] = np.sum(X**(i+j))
        B[i] = np.sum(Y*X**i)
    return np.concatenate((A, B), axis = 1)   
        
def solve(AB):
    n = AB.shape[0]
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

x = np.array([-2.5, 0.2, 2.1, 5.0, 7.8, 10.1])
y = np.array([1.5, 0.0, -0.8, -11.2, -10.85, -11.85])



coef1 = solve(createMatrix(x,y,1))
def f1(X):
    return coef1[0] + coef1[1]*X
print('Невязка 1 = ', sum((f1(x)-y)**2))

coef2 = solve(createMatrix(x,y,2))
def f2(X):
    return coef2[0] + coef2[1]*X + coef2[2]*X*X
print('Невязка 2 = ', sum((f2(x)-y)**2))

coef3 = solve(createMatrix(x,y,3))
def f3(X):
    return coef3[0] + coef3[1]*X + coef3[2]*X*X + coef3[3]*X*X*X
print('Невязка 3 = ', sum((f3(x)-y)**2))


x_a = np.arange(-8,15,0.1)
plt.figure(figsize=(10,5))
plt.plot(x_a,f1(x_a),'k',label ='1 степень')
plt.plot(x_a,f2(x_a),'b',label ='2 степень')
plt.plot(x_a,f3(x_a),'r',label ='3 степень')
plt.plot(x,y, 'go', label='Заданные точки')
plt.legend(loc='upper right')
plt.ylim(-15,5)
plt.xlim(-5,13)
