# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 12:57:17 2018

@author: stesn
"""

import numpy as np


def Determinant(A):
    n = A.shape[1]
    D = 1
    for i in range(n):
        if A[i,i] == 0:
            k = i + 1
            while k < n:
                if A[i,k] != 0:
                    A[[i,k],] = A[[k,i],]
                k = k + 1
        if A[i,i] != 0:
            D = D*A[i,i]
            A[i,:] = A[i,:]/A[i,i]
            j = i + 1
            while j < n:
                A[j,:] = A[j,:] - A[i,:]*(A[j,i]/A[i,i])
                j = j + 1
    return D

def Trapizoidal(h):
    n = int((b-a)/h)
    sum = 0
    x = a
    for i in range(n-1):
        x = x + h
        sum = sum + f(x)
    T = (h/2.0)*(f(a) + 2.0*sum + f(b))
    return T

def Simpson(h):
    n = int((b-a)/h)
    oddsum = 0
    evensum = 0
    x = a
    for i in range(1, n):
        x = x + h
        if (i%2 == 0):
            evensum += f(x)
        else:
            oddsum += f(x)
    S = (h/3.0)*(f(a) + 4.0*oddsum + 2*evensum + f(b))
    return S

def rectangle(h):
    n = int((b-a)/h)
    S = 0
    x = a
    for i in range(n):
        S = S + h*f((x+x+h)/2)
        x = x + h
    return S

def f(x):
  return (-4*x + 3)/(x*x + 4)

a = -7.
b = -3.
h_1 = 1.
h_2 = 0.5
h_3 = 0.25


p = 1
z_1 = rectangle(h_1)
z_2 = rectangle(h_2)
z_3 = rectangle(h_3)
D_1 = Determinant(np.array([[z_1, h_1**p, h_1**(p+1)],
                            [z_2, h_2**p, h_2**(p+1)],
                            [z_3, h_3**p, h_3**(p+1)]]))
D_2 = Determinant(np.array([[1, h_1**p, h_1**(p+1)],
                            [1, h_2**p, h_2**(p+1)],
                            [1, h_3**p, h_3**(p+1)]]))
z = D_1/D_2
print('Формула прямоугольников')
print('h_1=',h_1, ', z_1 =',z_1)
print('h_2=',h_2, ', z_2 =',z_2)
print('h_3=',h_3, ', z_3 =',z_3)
print('Уточнение по формуле Рунге-Ромберта', z)


p = 2
z_1 = Trapizoidal(h_1)
z_2 = Trapizoidal(h_2)
z_3 = Trapizoidal(h_3)
D_1 = Determinant(np.array([[z_1, h_1**p, h_1**(p+1)],
                            [z_2, h_2**p, h_2**(p+1)],
                            [z_3, h_3**p, h_3**(p+1)]]))
D_2 = Determinant(np.array([[1, h_1**p, h_1**(p+1)],
                            [1, h_2**p, h_2**(p+1)],
                            [1, h_3**p, h_3**(p+1)]]))
z = D_1/D_2
print('')
print('Формула трапеций')
print('h_1=',h_1, ', z_1 =',z_1)
print('h_2=',h_2, ', z_2 =',z_2)
print('h_3=',h_3, ', z_3 =',z_3)
print('Уточнение по формуле Рунге:', z)

p = 4
z_1 = Simpson(h_1)
z_2 = Simpson(h_2)
R = h_2/h_1
z = z_1 + (z_1 - z_2)/(R**p - 1)
print('')
print('Формула Симпсона')
print('h_1=',h_1, ', z_1 =',z_1)
print('h_2=',h_2, ', z_2 =',z_2)
print('Уточнение по формуле Рунге-Ромберта', z)
#print('Integral = {:.7f}'.format(T))
