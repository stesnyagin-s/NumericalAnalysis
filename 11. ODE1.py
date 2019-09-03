# -*- coding: utf-8 -*-
"""
Created on Sun Nov 18 13:47:06 2018

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

def euler(h):
    print('')
    print('h =',h)
    print('{:<2s} {:<11} {:<11}'.
      format('i', 'x','y'))
    n = int((b - a)/h)
    i = 0
    x = x_0
    y = y_0
    print('{:<2} {:<11.8f} {:<11.8f}'.
          format(i, x, y))
    for i in range(1,n+1):
        y += h*f(x,y)
        x += h
        print('{:<2} {:<11.8f} {:<11.8f}'.
          format(i, x, y))
    print('y(8) =',y)
    return y

def runge_kutta(h):
    print('')
    print('h =',h)
    print('{:<2s} {:<11} {:<15} {:<15} {:<15} {:<15} {:<15}'.
          format('i', 'x','y', 'K_1', 'K_2', 'K_3', 'K_4'))
    n = int((b - a)/h)
    i = 0
    x = x_0
    y = y_0
    K_1 = f(x,y)
    K_2 = f(x + h/2, y + (h/2)*K_1)
    K_3 = f(x + h/2, y + (h/2)*K_2)
    K_4 = f(x + h, y + h*K_3)
    print('{:<2} {:<11.8f} {:<15.8f} {:<15.8f} {:<15.8f} {:<15.8f} {:<15.8f}'.
          format(i, x, y, K_1, K_2, K_3, K_4))
    for i in range(1, n+1):
        y += (h/6)*(K_1 + 2*K_2 + 2*K_3 + K_4)
        x += h
        K_1 = f(x,y)
        K_2 = f(x + h/2, y + (h/2)*K_1)
        K_3 = f(x + h/2, y + (h/2)*K_2)
        K_4 = f(x + h, y + h*K_3)
        print('{:<2} {:<11.8f} {:<15.8f} {:<15.8f} {:<15.8f} {:<15.8f} {:<15.8f}'.
              format(i, x, y, K_1, K_2, K_3, K_4))
    print('y(8) =',y)
    return y

def f(x,y):
    return y - 3.*x*x + 5

a = 4.
b = 8.
h_1 = 1.
h_2 = 0.5
h_3 = 0.4
x_0 = 4.
y_0 = -4.


print('Метод Эйлера:')
p = 1
z_1 = euler(h_1)
z_2 = euler(h_2)
z_3 = euler(h_3)
D_1 = Determinant(np.array([[z_1, h_1**p, h_1**(p+1)],
                            [z_2, h_2**p, h_2**(p+1)],
                            [z_3, h_3**p, h_3**(p+1)]]))
D_2 = Determinant(np.array([[1, h_1**p, h_1**(p+1)],
                            [1, h_2**p, h_2**(p+1)],
                            [1, h_3**p, h_3**(p+1)]]))
z = D_1/D_2

print('\nУточнение y(8) по формуле Рунге:',z)

print('\n\nМетод Рунге-Кутты:')
z_1 = runge_kutta(h_1)
print()
z_2 = runge_kutta(h_2)
R = h_2/h_1
p = 4
z_pp = z_1 + (z_1 - z_2)/(R**p - 1)
print('\nУточнение y(8) по формуле Рунге-Ромберта:', z_pp)