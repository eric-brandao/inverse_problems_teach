# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 15:37:18 2023

@author: Admin
"""

import numpy as np # arrays
import matplotlib.pyplot as plt # plots
plt.rcParams.update({'font.size': 12})
import regu_test_problems as tp
import sys
sys.path.append('D:/Work/dev/insitu_sim_python/insitu')
import lcurve_functions as lc
#%%
n = 64
a = 0
b = 1
d = 0.25
#x = np.linspace(0, 1, n)



A,x = tp.gravity_model(n = n,a = a,b = b,d = d)

rhox = tp.density_sin(x, 0.5) + 0.5*tp.density_sin(x, 1.0)
b = A @ rhox
n = np.random.normal(loc = 0.0, scale = 0.001, size = len(b))
bn = b + n
#%%
[U,s,V] = lc.csvd(A)

#%%
#lc.plot_picard(U,s,b)
#lc.plot_picard(U,s,bn, noise_norm = np.linalg.norm(n))

#%%
rhox_k = lc.tsvd(U,s,V,bn,10)
rhox_tau = lc.ssvd(U,s,V,bn, np.linalg.norm(n))

plt.figure()
plt.plot(x, rhox_k)
plt.plot(x, rhox_tau)
plt.plot(x, rhox)
plt.ylim((-0.1,1.5))

#%%
# Model
A = np.array([[1.0, 1.0]])

# Measurement
b = 1.0
n = np.random.normal(loc = 0.0, scale = 0.11, size = 1)
bn = b + n

x = lc.cvx_solver(A, bn, np.linalg.norm(n))

#%%
lam_opt = lc.l_cuve(U, s, bn, plotit = True)