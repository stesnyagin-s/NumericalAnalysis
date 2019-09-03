# -*- coding: utf-8 -*-
"""
Created on Sat Nov 17 10:47:57 2018

@author: stesn
"""
def f(x):
  return 0.5*x*x*x + 5*x -4
a = 0.
b = 1.


eps = 0.01
i = 1
toRight = False
c = (a*f(b) - b*f(a))/(f(b) - f(a))
c_old = a
if (f(c)*f(b) < 0):
    toRight = True
print('{:<2s} {:<11} {:<11} {:<11} {:<11} {:<11} {:<11}'.
      format('i', 'a','c','b', 'f(a)', 'f(c)', 'f(b)'))
while True:
    c_old = c
    c = (a*f(b) - b*f(a))/(f(b) - f(a))
    print('{:<2} {:<11.8f} {:<11.8f} {:<11.8f} {:<11.8f} {:<11.8f} {:<11.8f}'.
    format(i, a, c, b, f(a), f(c), f(b)))
    if(not ((abs(c - c_old) > eps) or (abs(f(c_old)) > eps) or (i <= 5))):break
    if toRight:
        a = c
    else:
        b = c   
    i = i + 1