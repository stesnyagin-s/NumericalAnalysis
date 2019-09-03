# -*- coding: utf-8 -*-
"""
Created on Sun Nov 18 15:07:17 2018

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
    print('{:<2s} {:<11} {:<11} {:<11}'.
      format('i', 'x','y', 'z'))
    i = 0
    n = int((b - a)/h)
    x = x_0
    y = y_0
    z = z_0
    x_old = x
    y_old = y
    z_old = z
    print('{:<2} {:<11.8f} {:<11.8f} {:<11.8f}'.
          format(i, x, y, z))
    for i in range(1,n+1):
        x = x_old + h
        y = y_old + h*z_old
        z = z_old + h*f(x_old, y_old, z_old)
        print('{:<2} {:<11.8f} {:<11.8f} {:<11.8f}'.
          format(i, x, y, z))
        x_old = x
        y_old = y
        z_old = z
    print('y(1) =',y)
    return y

def runge_kutta(h):
    print('')
    print('h =',h)
    print('{:<2} {:<12} {:<12} {:<12} {:<6} {:<6} {:<6} {:<6} {:<6} {:<6} {:<6} {:<6}'.
              format('i', 'x', 'y', 'z', 'K_1y', 'K_1z', 'K_2y', 'K_2z','K_3y', 'K_3z','K_4y', 'K_4z'))
    i = 0
    n = int((b - a)/h)
    x = x_0
    y = y_0
    z = z_0
    K_1z = f(x,y,z)
    K_1y = z
    K_2z = f(x + h/2, y + (h/2)*K_1y,z + (h/2)*K_1z)
    K_2y = z + (h/2)*K_1z
    K_3z = f(x + h/2, y + (h/2)*K_2y,z + (h/2)*K_2z)
    K_3y = z + (h/2)*K_2z
    K_4z = f(x + h, y + h*K_3y,z + h*K_3z)
    K_4y = z + h*K_3z
    
    print('{:<2} {:<12.8f} {:<12.8f} {:<12.8f} {:<6.3f} {:<6.3f} {:<6.3f} {:<6.3f} {:<6.3f} {:<6.3f} {:<6.3f} {:<6.3f}'.
              format(i, x, y, z, K_1y, K_1z, K_2y, K_2z,K_3y, K_3z,K_4y, K_4z))
    for i in range(1, n+1):
        y = y + (h/6)*(K_1y + 2*K_2y + 2*K_3y + K_4y)
        z = z + (h/6)*(K_1z + 2*K_2z + 2*K_3z + K_4z)
        x = x + h
        K_1z = f(x,y,z)
        K_1y = z
        K_2z = f(x + h/2, y + (h/2)*K_1y,z + (h/2)*K_1z)
        K_2y = z + (h/2)*K_1z
        K_3z = f(x + h/2, y + (h/2)*K_2y,z + (h/2)*K_2z)
        K_3y = z + (h/2)*K_2z
        K_4z = f(x + h, y + h*K_3y,z + h*K_3z)
        K_4y = z + h*K_3z
        print('{:<2} {:<12.8f} {:<12.8f} {:<12.8f} {:<7.3f} {:<6.3f} {:<6.3f} {:<6.3f} {:<6.3f} {:<6.3f} {:<6.3f} {:<6.3f}'.
              format(i, x, y, z, K_1y, K_1z, K_2y, K_2z,K_3y, K_3z,K_4y, K_4z))
    print('y(1)', y)
    return y
    
def f(x,y,z):
    return -z + 2.*y - 2.*x - 5.
a = -1.
b = 1.
h_1 = 0.5
h_2 = 0.25
h_3 = 0.2
x_0 = -1.
y_0 = -1.
z_0 = -4.


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

print('\nУточнение y(1) по формуле Рунге:',z)

print('\n\nМетод Рунге-Кутты:')
z_1 = runge_kutta(h_1)
print()
z_2 = runge_kutta(h_2)
R = h_2/h_1
p = 4
z_pp = z_1 + (z_1 - z_2)/(R**p - 1)
print('\nУточнение y(1) по формуле Рунге-Ромберта:', z_pp)