# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 08:45:15 2023

@author: Admin
"""
import numpy as np # arrays
import matplotlib.pyplot as plt # plots
plt.rcParams.update({'font.size': 10})
import cvxpy as cp
import sklearn
import sys
sys.path.append('D:/Work/dev/insitu_sim_python/insitu')
import lcurve_functions as lc

from controlsair import compare_alpha, compare_spk
from controlsair import AlgControls, AirProperties, load_cfg, sph2cart
from sources import Source
# from receivers import Receiver
import receivers
from field_free import FreeField
from decompositionclass import Decomposition
from decomposition_ev_ig import DecompositionEv2
#%%
air = AirProperties(temperature = 20)
controls = AlgControls(c0 = air.c0, freq_vec = [1000])

receivers = Receiver()
receivers.random_3d_array(x_len=0.6, y_len=0.8, z_len=0.25, n_total = 290, zr = 0.1)
# receivers.double_planar_array(x_len = 0.3, n_x = 10, y_len = 0.3, n_y = 10, 
#                               zr = 0.01, dz = 0.01)
# receivers.brick_array(x_len = 0.3, n_x = 10, y_len = 0.3, n_y = 10, z_len = 0.1, n_z = 3, 
#                       zr = 0.1)
# receivers.planar_array(x_len = 0.3, n_x = 10, y_len = 0.3, n_y = 10, zr = 0.1)
#%%

field1 = FreeField(air, controls, receivers)
field1.planewave_ff(theta = np.deg2rad(-45), phi = np.deg2rad(30))

# field2 = FreeField(air, controls, receivers)
# field2.planewave_ff(theta = np.deg2rad(45), phi = np.deg2rad(30))

field = FreeField(air, controls, receivers)
field.planewave_ff(theta = np.deg2rad(45), phi = np.deg2rad(60))
field.pres_s[0] += field1.pres_s[0]
field.add_noise(snr = 30)
# field.plot_pres()
field.plot_scene(vsam_size=2)

#%%
# theta = np.deg2rad(-45)
# phi = np.deg2rad(30)
# s_coord1 = sph2cart(2, np.pi/2-theta, phi)
# theta = np.deg2rad(45)
# phi = np.deg2rad(60)
# s_coord2 = sph2cart(2, np.pi/2-theta, phi)
# source = Source(coord = s_coord1)
# source.add_sources(coord = s_coord2)

# field = FreeField(air, controls, receivers)
# field.monopole_ff(sources = source)
# field.pres_s[0] += field.pres_s[1]
# field.add_noise(snr = 30)
# field.plot_scene(vsam_size=2)
#%%

ff_ded = Decomposition(field.pres_s[0], controls = controls, 
                       receivers=receivers, regu_par = 'GCV')
ff_ded.wavenum_dir(n_waves = 600, plot = True)

ff_ded.pk_tikhonov(plot_l = True, method = 'Tikhonov')
ff_ded.pk_interpolate()
#%%
ff_ded.plot_pk_sphere(freq=2000, db=True, dinrange=12, travel=False)
ff_ded.plot_pk_map(freq=2000, db=True, dinrange=12)

#%%
ff_ded = Decomposition(field.pres_s[0], controls = controls, 
                       receivers=receivers, regu_par = 'GCV')
ff_ded.wavenum_dir(n_waves = 162, plot = False)
ff_ded.pk_cs(snr=30, headroom = 0)

#%%
recs = receivers.Receiver()
recs.hemispherical_array(radius = 1, n_rec_target = 642)

#%%
fig = plt.figure()
ax = plt.axes(projection ="3d")
ax.scatter(recs.coord[:,0], recs.coord[:,1], recs.coord[:,2])
#ax.scatter(recs.pdir_all[:,0], recs.pdir_all[:,1], recs.pdir_all[:,2])

#%%
pk = DecompositionEv2()
pk.prop_dir(n_waves = 642, plot = True)