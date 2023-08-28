# -*- coding: utf-8 -*-
"""
some test problems from regu

@author: Admin
"""
import numpy as np

def gravity_model(n,a,b,d):
    dx = 1/n
    dxl = (b-a)/n
    x = dx * (np.arange(n) + 0.5) 
    xl = a + dxl * (np.arange(n) + 0.5) 
    X,XL = np.meshgrid(x,xl)
    A = dx * (d/(d**2 + (X-XL)**2)**(3/2)) 
    return A, x

def density_pieciwise(n):
    nn = int(n/3)
    rhox = np.ones(n)
    rhox[:nn] = 2
    return rhox

def density_sin(x, fp):
    rhox = np.sin(2*np.pi*fp*x)
    return rhox

