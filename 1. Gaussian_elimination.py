# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 13:15:22 2018

@author: stesn
"""
import numpy as np
np.set_printoptions(suppress = True)
A = np.matrix([[-1., -5., -4., -6., 4],
               [-2., -10., -4., -2., -4.],
               [-3., -15., -12., -4., 5.],
               [10., -6., -11., 5., 4.],
               [-5., 2., -3., -4., 2.]])
B = np.matrix([[-25.],
               [18.],
               [-54.],
               [-111.],
               [20.]])


n = A.shape[1]
Determinant = 1
AB = np.concatenate((A, B), axis = 1)
for i in range(n):
    if AB[i,i] == 0:
        k = i + 1
        while k < n:
            if AB[i,k] != 0:
                AB[[i,k],] = AB[[k,i],]
            k = k + 1
    if AB[i,i] != 0:
        Determinant = Determinant*AB[i,i]
        AB[i,:] = AB[i,:]/AB[i,i]
        j = i + 1
        while j < n:
            AB[j,:] = AB[j,:] - AB[i,:]*(AB[j,i]/AB[i,i])
            j = j + 1
for i in range((n-1),-1, -1): 
    for j in range(i+1,n):
        AB[i,n] = AB[i,n] - AB[i, j]*AB[j,n]
        

E = np.identity(n)
AE = np.concatenate((A, E), axis = 1)
for i in range(n):
    if AE[i,i] == 0:
        k = i + 1
        while k < n:
            if AE[i,k] != 0:
                AE[[i,k],] = AE[[k,i],]
            k = k + 1
    if AE[i,i] != 0:
        Determinant = Determinant*AB[i,i]
        AE[i,:] = AE[i,:]/AE[i,i]
        j = i + 1
        while j < n:
            AE[j,:] = AE[j,:] - AE[i,:]*(AE[j,i]/AE[i,i])
            j = j + 1
for i in range(n-1,0, -1):
    AE[i,:] = AE[i,:]/AE[i,i]
    j = i - 1
    while(j >= 0):
      AE[j,:] = AE[j,:] - AE[i,:]*(AE[j,i]/AE[i,i])
      j = j - 1
print("Определитель:")
print(Determinant)
X = AB[:, n]
print('X = ')
print(X)
R = AE[:, range(n, 2*n)]
print("Обратная матрица:")
print(R)
print("A*X = ")
print(np.matmul(A,X))
print("A*A^-1 = ")
print(np.matmul(A,R))
