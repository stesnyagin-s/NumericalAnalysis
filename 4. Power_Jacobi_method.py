# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 21:31:50 2018

@author: stesn
"""

import numpy as np
np.set_printoptions(suppress = True, precision = 8)

def power_method():
    W = np.empty(shape = (n))
    W.fill(1/np.sqrt(n))
    print('W_0:\n',W)
    V = 0.0
    p = 0.0
    for i in range(1,11):
        V = A.dot(W)
        p = W.dot(V)
        W = V/np.linalg.norm(V)
        print('\ni =', i)
        print('V:\n', V)
        print('p=', p)
        print('W:\n', W)
    
def jacobi(A):
    V = np.identity(n)
    max = -np.inf
    k=m=0
    c = 1
    while (c <= 6):
        phi = 0
        for i in range(n):
            for j in range(i+1, n):
                if(np.abs(A[i,j]) > max):
                    max = np.abs(A[i,j])
                    k = i
                    m = j
        if(A[k,k] == A[m,m]):
            if(A[k,m] > 0):
                phi = np.pi/4
            else:
                phi = -np.pi/4
        else:
            phi = 0.5*np.arctan(2*A[k,m]/(A[k,k] - A[m,m]))
        print('\ni =',c)
        print('a_max =', A[k,m])
        print('phi = ', phi)
        H = np.identity(n)
        H[k,k] = np.cos(phi)
        H[m,k] = np.sin(phi)
        H[k,m] = -np.sin(phi)
        H[m,m] = np.cos(phi)
        print('H:\n', H)
        V = V.dot(H)
        A = np.matmul(np.matmul(H.transpose(), A),H)
        print('A:\n',A)
        c = c + 1
        max = -np.inf
    print('V:\n',V)
    V[:,0] = V[:,0]/np.linalg.norm(V[:,0], ord = np.inf)
    V[:,1] = V[:,1]/np.linalg.norm(V[:,1], ord = np.inf)
    V[:,2] = V[:,2]/np.linalg.norm(V[:,2], ord = np.inf)
    print('\nV_norn:\n',V)
    return A,V

A = np.array([[15., -11., 3.],
               [-11., 15., 8.],
               [3., 8., -20.]])
    
    
n = A.shape[1]
print('Степенной метод:')
power_method()
print('')
print('Метод Якоби:')
p, W = jacobi(A.copy())