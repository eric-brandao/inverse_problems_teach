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

rhox = tp.density_sin(x, 0.5) + 0.5*tp.density_sin(x, 1.0)# + 0.5*tp.density_sin(x, 2.0)
b = A @ rhox
n = np.random.normal(loc = 0.0, scale = 0.1, size = len(b))
bn = b + n
snr = 20*np.log10(np.linalg.norm(n)/np.linalg.norm(bn))
print('SNR = {}'.format(snr))
#%%
[U,s,V] = lc.csvd(A)

#%%
#lc.plot_picard(U,s,b)
#lc.plot_picard(U,s,bn, noise_norm = np.linalg.norm(n))

#%%
rhox_k = lc.tsvd(U,s,V,bn,7)
rhox_tau = lc.ssvd(U,s,V,bn, 0.5*np.linalg.norm(n))

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
lam_l = lc.l_cuve(U, s, bn, plotit = True, plot_curv = True)
# %%

lam_gcv = lc.gcv_lambda(U, s, bn, print_gcvfun = True)

#%%
x_delta, lam_dp = lc.discrep(U, s, V.T, bn, np.linalg.norm(n), x_0=None)

#%%
# lam_ncp, dist, reg_param = 
lam_gcv, dist = lc.ncp(U, s, bn)