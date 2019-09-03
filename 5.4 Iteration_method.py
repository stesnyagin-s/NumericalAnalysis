# -*- coding: utf-8 -*-
"""
Created on Sat Nov 17 11:06:18 2018

@author: stesn
"""
def f(x):
    return 3*x*x*x + 4*x + 16
def g(x):
    return -((4*x + 16)/3)**(1/3)
def dg(x):
    return -(4)/(3*3**(1/3)*((4*x + 16))**(2/3))
a = -2.
b = -1.


eps = 0.01
x = 0.0
x_old = float('Inf')
x_old_old = 0.0
i = 1
print('{:<2s} {:<11} {:<11} {:<11} {:<11}'.
      format('i', 'x','f(x)','g(x)', 'dg(x)'))

while((i <=5) or ((((x - x_old)**2)/abs(2*x_old - x - x_old_old)) > eps)\
      or (abs(f(x)) > eps)):
    x_old_old = x_old
    x_old = x
    
    if(abs(dg(a)) < 1):
        x =  a
        print('{:<2} {:<11.8f} {:<11.8f} {:<11.8f} {:<11.8f}'.
              format(i, x, f(x), g(x), dg(x)))
        a = g(x)
        
    elif(abs(dg(b)) < 1):
        x = b
        print('{:<2} {:<11.8f} {:<11.8f} {:<11.8f} {:<11.8f}'.
              format(i, x, f(x), g(x), dg(x)))
        b = g(x)
    i = i + 1