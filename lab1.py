# -*- coding: utf-8 -*-
"""
Created on Sun Nov 18 19:31:58 2018

@author: stesn
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from math import exp, sin
pd.set_option('precision',7)
pd.set_option('display.expand_frame_repr', False)


a = 0.00001
l = np.pi
t = 0.0
T = 4

def solve(A,B):
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
    return X



def yavnaya_shema(p,h,tau):
    N = int(l/h)
    M = int(T/tau + 1)
    Y = np.zeros(shape = (N+1))
    MAE = 0
    for i in range(N+1):
        x = a + i*h
        Y[i] = sin(x)
        MAE = MAE + np.abs(U(x,0) - Y[i])
    Y_0 = Y.copy()
    
    for j in range(1, M+1): 
        t = j*tau
        for i in range(1,N):
            Y[i] = (tau*a/(h*h))*(Y_0[i+1] - 2*Y_0[i] + Y_0[i-1]) + Y_0[i]
        if(p == 1):
            Y[0] = Y[1] - h*exp(-a*t)
            Y[N] = Y[N-1] - h*exp(-a*t)
        elif(p == 2):
            Y[0] = (-exp(-a*t)*2*h + 4*Y[1] - Y[2])/3
            Y[N] = (-2*h*exp(-a*t) - Y[N-2] + 4*Y[N-1])/3
        elif(p == 3):
            Y[0] = (Y[1] + (h*h/(2*a*tau))*Y_0[0] - h*exp(-a*t))/(1 + h*h/(2*a*tau))
            Y[N] = (Y[N-1] + (h*h/(2*a*tau))*Y_0[N] - h*exp(-a*t))/(1 + h*h/(2*a*tau))
        for i in range(N+1):
            x = i*h
            MAE = MAE + np.abs(U(x,t) - Y[i])
        Y_0 = Y.copy()
    MAE = MAE/((N+1)*(M+1))
    return Y,MAE

def neyavnaya_shema(p,h,tau):
    N = int(l/h)
    M = int(T/tau + 1)
    MAE = 0
    Y = np.zeros(shape = (N+1))
    for i in range(N+1):
        x = a + i*h
        Y[i] = sin(x)
        MAE = MAE + np.abs(U(x,0) - Y[i])
    Y_0 = Y.copy()
    
    for j in range(1, M+1): 
        t = j*tau
        A = np.zeros(shape = (N+1, N+1))
        B = np.zeros(shape = (N+1,))
        if(p == 1):
            A[0,0] = -1
            A[0,1] = 1
            B[0] = h*exp(-a*t)
            A[N,N-1] = -1
            A[N,N] = 1
            B[N] = -h*exp(-a*t)
        elif(p == 2):
            A[0,0] = -3
            A[0,1] = 4
            A[0,2] = -1
            B[0] = 2*h*exp(-a*t)
            A[N,N-2] = 1
            A[N,N-1] = -4
            A[N,N] = 3
            B[N] = -2*h*exp(-a*t)
        elif(p == 3):
            A[0,0] = -(1 + h*h/(2*a*tau))
            A[0,1] = 1
            B[0] = h*exp(-a*t) - (h*h/(2*a*tau))*Y_0[0]
            A[N,N-1] = -1
            A[N,N] = 1 + h*h/(2*a*tau)
            B[N] = -h*exp(-a*t) + (h*h/(2*a*tau))*Y_0[N]
            
        for i in range(1,N):
            A[i,i-1] = -a*tau/(h*h)
            A[i,i] = 1 + 2*a*tau/(h*h)
            A[i,i+1] = -a*tau/(h*h)
            B[i] = Y_0[i]
            
        c = (A[0,2]/A[1,2])
        A[0,:] = A[0,:] - A[1,:]*c
        B[0] = B[0] - B[1]*c
        c = (A[N,N-2]/A[N-1,N-2])
        A[N,:] = A[N,:] - A[N-1,:]*c
        B[N] = B[N] - B[N-1]*c
        
        Y = solve(A,B)
        for i in range(N+1):
            x = i*h
            MAE = MAE + np.abs(U(x,t) - Y[i])
        Y_0 = Y.copy()
    MAE = MAE/((N+1)*(M+1))
    return Y, MAE

def сrank_nicolson(p,h,tau):
    N = int(l/h)
    M = int(T/tau + 1)
    Y = np.zeros(shape = (N+1))
    MAE = 0
    for i in range(N+1):
        x = a + i*h
        Y[i] = sin(x)
        MAE = MAE + np.abs(U(x,0) - Y[i])
    Y_0 = Y.copy()
    for j in range(1, M+1): 
        t = j*tau
        A = np.zeros(shape = (N+1, N+1))
        B = np.zeros(shape = (N+1,))
        if(p == 1):
            A[0,0] = -1
            A[0,1] = 1
            B[0] = h*exp(-a*t)
            A[N,N-1] = -1
            A[N,N] = 1
            B[N] = -h*exp(-a*t)
        elif(p == 2):
            A[0,0] = -3
            A[0,1] = 4
            A[0,2] = -1
            B[0] = 2*h*exp(-a*t)
            A[N,N-2] = 1
            A[N,N-1] = -4
            A[N,N] = 3
            B[N] = -2*h*exp(-a*t)
        elif(p == 3):
            A[0,0] = -(1 + h*h/(2*a*tau))
            A[0,1] = 1
            B[0] = h*exp(-a*t) - (h*h/(2*a*tau))*Y_0[0]
            A[N,N-1] = -1
            A[N,N] = 1 + h*h/(2*a*tau)
            B[N] = -h*exp(-a*t) + (h*h/(2*a*tau))*Y_0[N]

        for i in range(1,N):
            A[i,i-1] = -1
            A[i,i] = 2*(1 + h*h/(tau*a))
            A[i,i+1] = 1
            B[i] = 2*(h*h/(a*tau) - 1)*Y_0[i] + Y_0[i+1] + Y_0[i-1]
        if(A[0,2] != 0):      
            c = (A[0,2]/A[1,2])
            A[0,:] = A[0,:] - A[1,:]*c
            B[0] = B[0] - B[1]*c
        if(A[N,N-2] != 0):
            c = (A[N,N-2]/A[N-1,N-2])
            A[N,:] = A[N,:] - A[N-1,:]*c
            B[N] = B[N] - B[N-1]*c
        
        Y = solve(A,B)
        for i in range(N+1):
            x = i*h
            MAE = MAE + np.abs(U(x,t) - Y[i])
        Y_0 = Y.copy()
    MAE = MAE/((N+1)*(M+1))
    return Y, MAE


# p = 1 - двухточечная первый порядок точности
# p = 2 - трехточечная второй порядок точности
# p = 3 - двухточечная второй порядок точности
p = 3
h = l/10
tau = h*h/2
Y_1,MAE_1 = yavnaya_shema(p,h,tau)
Y_2,MAE_2 = neyavnaya_shema(p,h,tau)
Y_3,MAE_3 = сrank_nicolson(p,h,tau)



U = lambda x, t : np.exp(-a*t)*np.sin(x)



X = np.arange(0, l+0.01, h)
Y_true = U(X,T)
data = pd.DataFrame(data = {'x': X, 'Явная': Y_1, 'Неявная': Y_2,\
                     'Кранка-Николсона':Y_3, 'Истинное':Y_true})
print('MAE явная схема:', MAE_1)
print('MAE неявная схема:', MAE_2)
print('MAE кранка-николсона схема:', MAE_3)

print(data)
plt.figure(figsize=(10,5))
plt.plot(X,Y_1,label = 'Явная схема')
plt.plot(X,Y_2,label = 'Неявная схема')
plt.plot(X,Y_3,label = 'Схема Кранка-Николсона')
plt.plot(X,Y_true,label = 'Аналитическое решение')
plt.legend(loc='upper left')
plt.xlabel('x')
plt.ylabel('u')

h_list = []
MAE_list = []
for i in range(2, 25):
    h_list.append(l/i)
    Y_3,MAE = сrank_nicolson(p,l/i,tau)
    MAE_list.append(MAE)
    
plt.figure(figsize=(10,5))
plt.plot(h_list,MAE_list)
plt.xlabel('h')
plt.ylabel('MAE')
plt.title('Зависимость средней абсолютной ошибки схемы Кранка-Николсона от h')
tau_list = []
MAE_list = []
for i in range(3, 25):
    tau_list.append(T/i)
    Y_3,MAE = сrank_nicolson(p,h,T/i)
    MAE_list.append(MAE)
    
plt.figure(figsize=(10,5))
plt.plot(tau_list,MAE_list)
plt.xlabel('tau')
plt.ylabel('MAE')
plt.title('Зависимость средней абсолютной ошибки схемы Кранка-Николсона от tau')