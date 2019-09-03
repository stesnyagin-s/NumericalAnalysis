# -*- coding: utf-8 -*-
"""
Created on Sun Nov 18 19:31:58 2018

@author: stesn
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
pd.set_option('precision',8)
pd.set_option('display.expand_frame_repr', False)

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



def yavnaya_shema(p):
    Y = np.zeros(shape = (N+1))
    for i in range(N+1):
        Y[i] = f_3(a + i*h)
    Y_old = Y.copy()
    
    for j in range(1, M+1): 
        t = j*tau
        for i in range(1,N):
            Y[i] = (alpha_1*tau/(h*h) - alpha_2*tau/(2*h))*Y_old[i-1]\
            + (tau*alpha_3 + 1 - 2*alpha_1*tau/(h*h))*Y_old[i]\
            + (alpha_1*tau/(h*h) + alpha_2*tau/(2*h))*Y_old[i+1]\
            + tau*f(a + i*h,(j-1)*tau)
        if(p == 1):
            Y[0] = (h*f_1(j*tau) - phi_1*Y[1])/(h*phi_2 - phi_1)
            Y[N] = (h*f_2(j*tau) + phi_4*Y[N-1])/(h*phi_5 + phi_4)
        elif(p == 2):
            Y[0] = (2*h*f_1(t) - 4*phi_1*Y[1] + phi_1*Y[2])/(2*h*phi_2 - 3*phi_1)
            Y[N] = (2*h*f_2(t) + 4*phi_4*Y[N-1] - phi_4*Y[N-2])/(2*h*phi_5 + 3*phi_4)
        elif(p == 3):
            Y[0] = (Y[1]*(2*alpha_1/h) + (h/tau)*Y_old[0] + h*f(a,t) - ((2*alpha_1 - alpha_2*h)/phi_1)*f_1(t))/(2*alpha_1/h + h/tau - alpha_3*h - (phi_2/phi_1)*(2*alpha_1 - alpha_2*h))
            Y[N] = (Y[N-1]*(2*alpha_1/h) + (h/tau)*Y_old[N] + h*f(b,t) + ((2*alpha_1 + alpha_2*h)/phi_4)*f_2(t))/(2*alpha_1/h + h/tau - alpha_3*h + (phi_5/phi_4)*(2*alpha_1 + alpha_2*h))
        Y_old = Y.copy()
    return Y

def neyavnaya_shema(p):
    Y = np.zeros(shape = (N+1))
    for i in range(N+1):
        Y[i] = f_3(a + i*h)
    Y_old = Y.copy()
    
    for j in range(1, M+1): 
        t = j*tau
        A = np.zeros(shape = (N+1, N+1))
        B = np.zeros(shape = (N+1,))
        if(p == 1):
            A[0,0] = (h*phi_2 - phi_1)
            A[0,1] = phi_1
            B[0] = h*f_1(t)
            A[N,N-1] = -phi_4
            A[N,N] = h*phi_5 + phi_4
            B[N] = h*f_2(t)
        elif(p == 2):
            A[0,0] = 2*h*phi_2 - 3*phi_1
            A[0,1] = 4*phi_1
            A[0,2] = -phi_1
            B[0] = 2*h*f_1(t)
            A[N,N-2] = phi_4
            A[N,N-1] = -4*phi_4
            A[N,N] = 2*h*phi_5 + 3*phi_4
            B[N] = 2*h*f_2(t)
        elif(p == 3):
            A[0,0] = 2*alpha_1/h + h/tau - alpha_3*h\
                    - (phi_2/phi_1)*(2*alpha_1 - alpha_2*h)
            A[0,1] = -2*alpha_1/h
            B[0] = (h/tau)*Y_old[0] + h*f(a,t) - ((2*alpha_1 - alpha_2*h)/phi_1)*f_1(t)
            A[N,N-1] = -2*alpha_1/h
            A[N,N] = 2*alpha_1/h + h/tau - alpha_3*h\
                    + (phi_5/phi_4)*(2*alpha_1 + alpha_2*h)
            B[N] = (h/tau)*Y_old[N] + h*f(b,t) + ((2*alpha_1 + alpha_2*h)/phi_4)*f_2(t)
            
        for i in range(1,N):
            A[i,i-1] = (2*alpha_1 - h*alpha_2)*(tau/(2*h*h))
            A[i,i] = alpha_3*tau - 2*alpha_1*tau/(h*h) - 1
            A[i,i+1] = (2*alpha_1 + h*alpha_2)*tau/(2*h*h)
            B[i] = -Y_old[i] - tau*f(a + i*h, j*tau)
            
        c = (A[0,2]/A[1,2])
        A[0,:] = A[0,:] - A[1,:]*c
        B[0] = B[0] - B[1]*c
        c = (A[N,N-2]/A[N-1,N-2])
        A[N,:] = A[N,:] - A[N-1,:]*c
        B[N] = B[N] - B[N-1]*c
        
        Y = solve(A,B)
        Y_old = Y.copy()
    return Y

def сrank_nicolson(p):
    theta = 0.5
    
    Y = np.zeros(shape = (N+1))
    for i in range(N+1):
        Y[i] = f_3(a + i*h)
    Y_old = Y.copy()
    
    for j in range(1, M+1): 
        t = j*tau
        A = np.zeros(shape = (N+1, N+1))
        B = np.zeros(shape = (N+1,))
        
        if(p == 1):
            A[0,0] = (h*phi_2 - phi_1)
            A[0,1] = phi_1
            B[0] = h*f_1(t)
            A[N,N-1] = -phi_4
            A[N,N] = h*phi_5 + phi_4
            B[N] = h*f_2(t)
        elif(p == 2):
            A[0,0] = 2*h*phi_2 - 3*phi_1
            A[0,1] = 4*phi_1
            A[0,2] = -phi_1
            B[0] = 2*h*f_1(t)
            A[N,N-2] = phi_4
            A[N,N-1] = -4*phi_4
            A[N,N] = 2*h*phi_5 + 3*phi_4
            B[N] = 2*h*f_2(t)
        elif(p == 3):
            A[0,0] = 2*alpha_1/h + h/tau - alpha_3*h\
                    - (phi_2/phi_1)*(2*alpha_1 - alpha_2*h)
            A[0,1] = -2*alpha_1/h
            B[0] = (h/tau)*Y_old[0] + h*f(a,t) - ((2*alpha_1 - alpha_2*h)/phi_1)*f_1(t)
            A[N,N-1] = -2*alpha_1/h
            A[N,N] = 2*alpha_1/h + h/tau - alpha_3*h\
                    + (phi_5/phi_4)*(2*alpha_1 + alpha_2*h)
            B[N] = (h/tau)*Y_old[N] + h*f(b,t) + ((2*alpha_1 + alpha_2*h)/phi_4)*f_2(t)

        for i in range(1,N):
            A[i,i-1] = (2*alpha_1 - h*alpha_2)*(tau*theta)/(2*h*h)
            A[i,i] = alpha_3*theta*tau - 2*alpha_1*theta*tau/(h*h) - 1
            A[i,i+1] = (2*alpha_1 + h*alpha_2)*(tau*theta)/(2*h*h)
            B[i] = -theta*tau*f(a + i*h, j*tau) - Y_old[i]\
              - (1 - theta)*(alpha_1*tau*(Y_old[i+1] - 2*Y_old[i] + Y[i-1])/(h*h)\
                 + alpha_2*tau*(Y_old[i+1] - Y_old[i-1])/(2*h)\
                 + alpha_3*tau*Y_old[i] + tau*f(a + i*h, (i-1)*tau))
              
        c = (A[0,2]/A[1,2])
        A[0,:] = A[0,:] - A[1,:]*c
        B[0] = B[0] - B[1]*c
        c = (A[N,N-2]/A[N-1,N-2])
        A[N,:] = A[N,:] - A[N-1,:]*c
        B[N] = B[N] - B[N-1]*c
        
        Y = solve(A,B)
        Y_old = Y.copy()
    return Y

def f(x,t):
    return (x+2)/(t + 11)
def f_1(t):
    return 10. + t/2.
def f_2(t):
    return 20. + t*t/4. + 0.2*t
def f_3(x):
    return -1. + 4.*np.sin(np.pi/2*x )

alpha_1 = 5.
alpha_2 = -5.
alpha_3 = -1.5
phi_1 = 2.
phi_2 = -2.
phi_4 = 4.
phi_5 = -4.

a = -1.
b = 3.
h = 0.8
tau = 0.025
t = 0.0
T = 0.1
N = int((b - a)/h)
M = int(T/tau)

# p = 1 - двухточечная первый порядок точности
# p = 2 - трехточечная второй порядок точности
# p = 3 - двухточечная второй порядок точности

p = 3
Y_1 = yavnaya_shema(p)
Y_2 = neyavnaya_shema(p)
Y_3 = сrank_nicolson(p)
X = np.arange(a, b+0.01, h)
data = pd.DataFrame(data = {'x': X, 'Явная схема': Y_1, 'Неявная схема': Y_2,\
                     'Схема Кранка-Николсона':Y_3})
print(data)
plt.figure(figsize=(10,5))
plt.plot(X,Y_1,label = 'Явная схема')
plt.plot(X,Y_2,label = 'Неявная схема')
plt.plot(X,Y_3,label = 'Схема Кранка-Николсона')
plt.legend(loc='upper left')
plt.xlabel('x')
plt.ylabel('y(0.1, x)')
