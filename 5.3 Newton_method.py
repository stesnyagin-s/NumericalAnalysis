# -*- coding: utf-8 -*-
"""
Created on Sat Nov 17 10:47:57 2018

@author: stesn
"""
def f(x):
    return 0.5*x*x*x + 10.*x - 20.
def df(x):
    return 1.5*x*x + 10.
def ddf(x):
    return 3.*x
a = 1.
b = 2.


eps = 0.01
i = 0
x = 0
x_old = 0
if (f(a)*ddf(a) > 0):
    x_old = b
    x = a
elif (f(b)*ddf(b) > 0):
    x_old = a
    x = b
print('{:<2s} {:<11} {:<11} {:<11} {:<11} {:<11}'.
      format('n', 'a_(n-1)','f(a_(n-1))','df(a_(n-1))', 'ddf', 'a_n'))
while True:
    x_old = x
    x = x_old - f(x_old)/df(x_old)
    print('{:<2} {:<11.8f} {:<11.8f} {:<11.8f} {:<11.8f} {:<11.8f}'.
    format(i, x_old, f(x_old), df(x_old), ddf(x_old), x))
    if(not ((abs(x - x_old) > eps) or (abs(f(x)) > eps) or (i <= 5))):break
    i = i + 1