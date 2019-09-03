# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 21:47:27 2019

@author: stesn
"""
 
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
from matplotlib import cm

def seidel(N,M):
    h_y = (b_y - a_y)/M
    h_x = (b_x - a_x)/N
    X = np.zeros(shape = (N+1,M+1))
    for i in range(N+1):
        x = a_x + h_x*i
        X[i,0] = np.sin(x)
        X[i,M] = np.e*np.sin(x)
    for j in range(1,M):
        y = a_y + h_y*j
        X[0,j]=-np.exp(y)*h_x
        X[N,j]=-np.exp(y)*h_x
    while True:
        X_old = X.copy()
        for j in range(1,M):
            for i in range(0,N+1):
                if(i == 0):
                    y = a_y + h_y*j
                    X[0,j] = X[1,j] - np.exp(y)*h_x
                elif(i == N):
                    y = a_y + h_y*j
                    X[N,j] = X[N-1,j] - np.exp(y)*h_x
                else:
                    X[i,j] = (1/(2*(1/(h_x*h_x) + 1/(h_y*h_y))))*((X[i+1,j]\
                     + X[i-1,j])/(h_x*h_x) + (X[i,j+1] + X[i,j-1])/(h_y*h_y))
        
        if(np.max(np.abs(X-X_old)) < 0.001):break
    return X

def iteration(N,M):
    h_y = (b_y - a_y)/N
    h_x = (b_x - a_x)/N
    X = np.empty(shape = (N+1,M+1))
    for i in range(N+1):
        x = a_x + h_x*i
        X[i,0] = np.sin(x)
        X[i,M] = np.e*np.sin(x)
    for i in range(1,N):
        for j in range(1,M):
            X[i,j] = X[i,0] + h_y*j*(X[i,M] - X[i,0])
    for j in range(1,M):
        y = a_y + h_y*j
        X[0,j] = X[1,j] - h_x*np.exp(y)
        X[N,j] = X[N-1,j] - h_x*np.exp(y)
    X_old = 0
    while True:
        X_old = X.copy()
        for i in range(1,N):
            for j in range(1,M):
                X[i,j] = (1/(2*(1/(h_x*h_x) + 1/(h_y*h_y))))*((X_old[i+1,j]\
                 + X_old[i-1,j])/(h_x*h_x) + (X_old[i,j+1] + X_old[i,j-1])/(h_y*h_y))
        for j in range(1,M):
            y = a_y + h_y*j
            X[0,j] = X_old[1,j] - np.exp(y)*h_x
            X[N,j] = X_old[N-1,j] - np.exp(y)*h_x
        if(np.max(np.abs(X-X_old)) < 0.001):break
    return X

def relaxation(N,M):
    w = 1.1
    h_y = (b_y - a_y)/M
    h_x = (b_x - a_x)/N
    X = np.zeros(shape = (N+1,M+1))
    for i in range(N+1):
        x = a_x + h_x*i
        X[i,0] = np.sin(x)
        X[i,M] = np.e*np.sin(x)
    for j in range(1,M):
        y = a_y + h_y*j
        X[0,j]=-np.exp(y)*h_x
        X[N,j]=-np.exp(y)*h_x
    while True:
        X_old = X.copy()
        for j in range(1,M):
            for i in range(0,N+1):
                if(i == 0):
                    y = a_y + h_y*j
                    X[0,j] = w*(X[1,j] - np.exp(y)*h_x) + (1 - w)*X[0,j]
                elif(i == N):
                    y = a_y + h_y*j
                    X[N,j] = w*(X[N-1,j] - np.exp(y)*h_x) + (1 - w)*X[N,j]
                else:
                    X[i,j] = w*((1/(2*(1/(h_x*h_x) + 1/(h_y*h_y))))*((X[i+1,j]\
                     + X[i-1,j])/(h_x*h_x) + (X[i,j+1] + X[i,j-1])/(h_y*h_y))) + (1 - w)*X[i,j]
        
        if(np.max(np.abs(X-X_old)) < 0.001):break
    return X

def MAE(X):
    N = X.shape[0]-1
    M = X.shape[1]-1
    h_x = (b_x - a_x)/N
    h_y = (b_y - a_y)/M
    sum = 0
    for i in range(N+1):
        for j in range(M+1):
            sum = sum + abs(X[i,j]-U(a_x + i*h_x, a_y + j*h_y))
    return sum/((N+1)*(M+1))
    

a_x = 0
b_x = np.pi
a_y = 0
b_y = 1
N = 30
M = 30

U = lambda x,y : np.sin(x)*np.exp(y)
X_1 = iteration(N,M)
X_2 = seidel(N,M)
X_3 = relaxation(N,M)
print("Метод простых итераций:")
print('MAE:', MAE(X_1))
print("Метод Зейделя:")
print('MAE:', MAE(X_2))
print("Метод релаксации:")
print('MAE:', MAE(X_3))

h_x_list = []
MAE_list = []
for i in range(4,50):
    h_x_list.append((b_x - a_x)/i)
    MAE_list.append(MAE(seidel(i,M)))
plt.figure(figsize=(10,5))
plt.plot(h_x_list,MAE_list)
plt.xlabel('h_x')
plt.ylabel('MAE')
plt.title('Зависимость средней абсолютной ошибки метода Зейделя от h_x')

h_y_list = []
MAE_list = []
for j in range(3,40):
    h_y_list.append((b_y - a_y)/j)
    MAE_list.append(MAE(seidel(N,j)))
plt.figure(figsize=(10,5))
plt.plot(h_y_list,MAE_list)
plt.xlabel('h_y')
plt.ylabel('MAE')
plt.title('Зависимость средней абсолютной ошибки метода Зейделя от h_y')

X = np.linspace(a_x,b_x,N+1)
Y = np.linspace(a_y,b_y,M+1)
X, Y = np.meshgrid(X, Y)
fig = plt.figure(figsize=(13,10))
ax = Axes3D(fig)
surf = ax.plot_wireframe(X, Y, X_3)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('U(x,y)')
ax.set_title('Метод Зейделя')
