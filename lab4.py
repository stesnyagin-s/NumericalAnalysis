# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 01:55:05 2019

@author: stesn
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
from matplotlib import cm


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

def per_napr(N_x,N_y,N_tau):
    h_x = (b_x - a_x)/N_x
    h_y = (b_y - a_y)/N_y
    tau = T/N_tau
    XY_old = np.empty(shape = (N_x + 1, N_y + 1))
    XY = np.zeros(shape = (N_x + 1, N_y + 1))
    for i in range(N_x + 1):
        for j in range(N_y + 1):
            x = a_x + i*h_x
            y = a_y + j*h_y
            XY_old[i,j] = np.cos(2*x)*np.cosh(y)
    for k in range(N_tau + 1):
        t = tau*(k - 1) + tau/2
        A = np.zeros(shape = (N_x + 1, N_x + 1))
        B = np.zeros(shape = (N_x + 1))
        for j in range(1,N_y):
            for i in range(1,N_x):
                A[i,i-1] = -a/(h_x*h_x)
                A[i,i] = (2*a/(h_x*h_x) + 2/tau)
                A[i,i+1] = -a/(h_x*h_x)
                B[i] = (a/(h_y*h_y))*(XY_old[i,j+1] - 2*XY_old[i,j] + XY_old[i,j-1])\
                    + (2/tau)*XY[i,j]
            y = a_y + j*h_y
            A[0,0] = 1
            B[0] = np.cosh(y)*np.exp(-3*a*t)
            A[N_x,N_x] = 1
            B[N_x] = 0
            XY[:,j] = solve(A,B)
        for i in range(N_x + 1):
            x = a_x + i*h_x
            XY[i,0] = np.cos(2*x)*np.exp(-3*a*t)
            XY[i,N_y] = h_y*(3/4)*np.cos(2*x)*np.exp(-3*a*t) + XY[i,N_y -1]
        XY_old = XY.copy()
        
        t = tau*k
        A = np.zeros(shape = (N_y + 1, N_y + 1))
        B = np.zeros(shape = (N_y + 1))
        for i in range(1,N_x):
            for j in range(1,N_y):
                A[j,j-1] = -a/(h_y*h_y)
                A[j,j] = (2*a/(h_y*h_y) + 2/tau)
                A[j,j+1] = -a/(h_y*h_y)
                B[j] = (a/(h_x*h_x))*(XY_old[i+1,j] - 2*XY_old[i,j] + XY_old[i-1,j])\
                    + (2/tau)*XY[i,j]
            x = a_x + i*h_x
            A[0,0] = 1
            B[0] = np.cos(2*x)*np.exp(-3*a*t)
            A[N_y-1,N_y] = -1
            A[N_y,N_y] = 1
            B[N_y] = h_y*(3/4)*np.cos(2*x)*np.exp(-3*a*t)
            XY[i,:] = solve(A,B)
        for j in range(N_y + 1):
            y = a_y + j*h_y
            XY[0,j] = np.cosh(y)*np.exp(-3*a*t)
            XY[N_x,j] = 0
        XY_old = XY.copy()
    return XY

def drob_shag(N_x,N_y,N_tau):
    h_x = (b_x - a_x)/N_x
    h_y = (b_y - a_y)/N_y
    tau = T/N_tau
    XY_old = np.empty(shape = (N_x + 1, N_y + 1))
    XY = np.zeros(shape = (N_x + 1, N_y + 1))
    for i in range(N_x + 1):
        for j in range(N_y + 1):
            x = a_x + i*h_x
            y = a_y + j*h_y
            XY_old[i,j] = np.cos(2*x)*np.cosh(y)
    for k in range(N_tau + 1):
        t = tau*(k - 1) + tau/2
        A = np.zeros(shape = (N_x + 1, N_x + 1))
        B = np.zeros(shape = (N_x + 1))
        for j in range(1,N_y):
            for i in range(1,N_x):
                A[i,i-1] = a/(h_x*h_x)
                A[i,i] = -(2*a/(h_x*h_x) + 1/tau)
                A[i,i+1] = a/(h_x*h_x)
                B[i] = -XY_old[i,j]/tau
            y = a_y + j*h_y
            A[0,0] = 1
            B[0] = np.cosh(y)*np.exp(-3*a*t)
            A[N_x,N_x] = 1
            B[N_x] = 0
            XY[:,j] = solve(A,B)
        for i in range(N_x + 1):
            x = a_x + i*h_x
            XY[i,0] = np.cos(2*x)*np.exp(-3*a*t)
            XY[i,N_y] = h_y*(3/4)*np.cos(2*x)*np.exp(-3*a*t) + XY[i,N_y -1]
        XY_old = XY.copy()
        
        t = tau*k
        A = np.zeros(shape = (N_y + 1, N_y + 1))
        B = np.zeros(shape = (N_y + 1))
        for i in range(1,N_x):
            for j in range(1,N_y):
                A[j,j-1] = a/(h_y*h_y)
                A[j,j] = -(2*a/(h_y*h_y) + 2/tau)
                A[j,j+1] = a/(h_y*h_y)
                B[j] = -XY_old[i,j]/tau
            x = a_x + i*h_x
            A[0,0] = 1
            B[0] = np.cos(2*x)*np.exp(-3*a*t)
            A[N_y-1,N_y] = -1
            A[N_y,N_y] = 1
            B[N_y] = h_y*(3/4)*np.cos(2*x)*np.exp(-3*a*t)
            XY[i,:] = solve(A,B)
        for j in range(N_y + 1):
            y = a_y + j*h_y
            XY[0,j] = np.cosh(y)*np.exp(-3*a*t)
            XY[N_x,j] = 0
        XY_old = XY.copy()
    return XY
    
def MAE(X):
    N = X.shape[0]-1
    M = X.shape[1]-1
    h_x = (b_x - a_x)/N
    h_y = (b_y - a_y)/M
    sum = 0
    for i in range(N+1):
        for j in range(M+1):
            sum = sum + abs(X[i,j]-U(a_x + i*h_x, a_y + j*h_y, T))
    return sum/((N+1)*(M+1))

a = 0.05
a_x = 0
b_x = np.pi/4
a_y = 0
b_y = np.log(2)
T = 10
U = lambda x,y,t: np.cos(2*x)*np.cosh(y)*np.exp(-3*a*t)
N_x = 50
N_y = 50
N_tau = 50
print('N_x:', N_x)
print('N_y:', N_y)
print('N_tau:', N_tau)
print('Метод переменных направлений, MAE:', MAE(per_napr(N_x,N_y,N_tau)))
print('Метод дробных шагов, MAE:', MAE(drob_shag(N_x,N_y,N_tau)))

#h_x_list = []
#MAE_list = []
#for i in range(4,50):
#    h_x_list.append((b_x - a_x)/i)
#    MAE_list.append(MAE(per_napr(i,N_y,N_tau)))
#plt.figure(figsize=(10,5))
#plt.plot(h_x_list,MAE_list)
#plt.xlabel('h_x')
#plt.ylabel('MAE')
#plt.title('Зависимость средней абсолютной метода переменных направлений ошибки от h_x')
#
#h_y_list = []
#MAE_list = []
#for j in range(3,50):
#    h_y_list.append((b_y - a_y)/j)
#    MAE_list.append(MAE(drob_shag(N_x,j,N_tau)))
#plt.figure(figsize=(10,5))
#plt.plot(h_y_list,MAE_list)
#plt.xlabel('h_y')
#plt.ylabel('MAE')
#plt.title('Зависимость средней абсолютной ошибки метода переменных направлений от h_y')
#
#tau_list = []
#MAE_list = []
#for j in range(3,50):
#    tau_list.append((b_y - a_y)/j)
#    MAE_list.append(MAE(drob_shag(N_x,N_y,j)))
#plt.figure(figsize=(10,5))
#plt.plot(tau_list,MAE_list)
#plt.xlabel('tau')
#plt.ylabel('MAE')
#plt.title('Зависимость средней абсолютной ошибки  метода переменных направлений от tau')

X = np.linspace(a_x,b_x,N_x +1)
Y = np.linspace(a_y,b_y,N_y +1)
X, Y = np.meshgrid(X, Y)
fig = plt.figure(figsize=(13,10))
ax = Axes3D(fig)
surf = ax.plot_wireframe(X, Y, per_napr(N_x,N_y,N_tau))
ax.view_init(elev=25., azim=110,)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('U(x,y,10)')
ax.set_title('Метод переменных направлений')
