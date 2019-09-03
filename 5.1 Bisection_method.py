# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 11:24:23 2018

@author: stesn
"""

def f(x):
  return 0.5*x*x*x + 4.*x + 2.2
a = -1.
b = 0.


eps = 0.01
i = 1
print('{:<2s} {:<11} {:<11} {:<11} {:<11} {:<11} {:<11}'.
      format('i', 'a','c','b', 'f(a)', 'f(c)', 'f(b)'))
while True:
    c = (a + b)/2.
    print('{:<2} {:<11.8f} {:<11.8f} {:<11.8f} {:<11.8f} {:<11.8f} {:<11.8f}'.
    format(i, a, c, b, f(a), f(c), f(b)))
    if((b - a) < 2*eps):
        break
    if f(a)*f(c) < 0:
        b = c
    else:
        a = c
    i = i + 1