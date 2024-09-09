import sys
sys.path.append('D:/Work/dev/insitu_sim_python/insitu')
import decomp_quad_EU as dqdt_eu
import decomp_quad as dqdt
import decomp_quad_v2 as dqdt_v2
#from decomp2waves import Decomposition_2W  # Plane waves
from decomp2mono import Decomposition_2M  # Monopoles
from controlsair import AirProperties, AlgControls, sph2cart  # Controls
from material import PorousAbsorber  # Material
from field_bemflush import BEMFlushSq  # Field
from field_inf_nlr import NLRInfSph  # Field Inf NLR
from receivers import Receiver  # Receivers
import matplotlib.pyplot as plt  # Plot
from sources import Source  # Source
import numpy as np  # Numpy library
plt.rcParams.update({'font.size': 14})

import lcurve_functions as lc

#%% NLR Field ref
r = 0.3
theta = 0 
phi = 0 

main_folder = 'D:/Work/UFSM/Disciplinas/Problemas Inversos/10_DCISM/saved_fields/'

f_name_field_nlr_zs = 'NLR_zs_d10cm_el0d_az0d_r30cm_resist5k_z__0'

field_nlr_zs = NLRInfSph()  # NLR field
field_nlr_zs.load(path=main_folder, filename=f_name_field_nlr_zs)

# NLR
zs_nlr = field_nlr_zs.pres_s[0]/field_nlr_zs.uz_s[0]  # Zs = P / U
zs_nlr_mean = np.zeros(len(field_nlr_zs.controls.freq), dtype=complex)  # Mean of Zs
alpha_nlr = np.zeros(len(field_nlr_zs.controls.freq))  # Absorption coefficient

for jf in range(0, len(field_nlr_zs.controls.freq)):  # Absorption coefficient estimation
    zs_nlr_mean[jf] = np.mean(zs_nlr[:, jf])
    alpha_nlr[jf] = 1 - (np.abs(np.divide((zs_nlr_mean[jf] * np.cos(theta) - 1),
                                          (zs_nlr_mean[jf] * np.cos(theta) + 1)))) ** 2
#%% Field to decomp

f_name_field_nlr = 'NLR_Inf_d10cm_el0d_az0d_r30cm_resist5k_r3d_128_mics_30x30x5cm'

field_nlr = NLRInfSph()  # NLR field
field_nlr.load(path=main_folder, filename=f_name_field_nlr)
snr = 40
field_nlr.add_noise(snr=snr, uncorr=False)  # array


#%%
ng = 30
a = 0
b = 30
retraction = 0

source_coord = sph2cart(r, np.pi / 2 - theta, phi)  # Source coordinates
#source = Source(coord=source_coord)

decomp_qdt_v2 = dqdt_v2.Decomposition_QDT(p_mtx=field_nlr.pres_s[0], controls=field_nlr.controls,
                               receivers=field_nlr.receivers, source_coord=source_coord, quad_order=ng,
                                            a = a, b = b, retraction = retraction, image_source_on = True,
                                            regu_par = 'l-curve')
decomp_qdt_v2.gauss_legendre_sampling()
#decomp_qdt_v2.gauss_laguerre_sampling()
#decomp_qdt_v2.uniform_sampling()
#decomp_qdt_v2.mid_point_sampling()

#r, zr, r1, rq = decomp_qdt_v2.get_rec_parameters(field_nlr.receivers)
#iq = decomp_qdt_v2.kernel_p(18, rq)
decomp_qdt_v2.pk_tikhonov(plot_l=False, method='Tikhonov')
#decomp_qdt_v2.least_squares_pk()
# decomp_qdt_v2.reconstruct_p(field_nlr_zs.receivers)


#r, zr, r1, rq = decomp_qdt_v2.get_rec_parameters(field_nlr_zs.receivers)
#dqdt_v2.kernel_uz(18, rq, decomp_qdt_v2.hs + zr - 1j * decomp_qdt_v2.q)
# decomp_qdt_v2.reconstruct_uz(field_nlr_zs.receivers)
decomp_qdt_v2.zs(Lx=0.1, n_x=21, Ly=0.1, n_y=21, theta=[np.deg2rad(0)], avgZs=True);  # Zs
nmse_qdt_ps = lc.nmse_freq(decomp_qdt_v2.p_recon , field_nlr_zs.pres_s[0])
nmse_qdt_uzs = lc.nmse_freq(decomp_qdt_v2.uz_recon , field_nlr_zs.uz_s[0])

# pk = np.array(decomp_qdt_v2.pk)
# alpha_pk = np.zeros(15)
# for jf in np.arange(15):
#     Qref = np.sum(np.abs(pk[jf,1])**2)/np.abs(pk[jf,0])**2
#     alpha_pk[jf] = 1-(np.abs(Qref))

#%%
# decomp_qdt = dqdt.Decomposition_QDT(p_mtx=field_nlr.pres_s[0], controls=field_nlr.controls,
#                                receivers=field_nlr.receivers, source_coord=source_coord, quad_order=ng)

# decomp_qdt.pk_tikhonov(plot_l=True, method='Tikhonov', a=a, b=b, retraction=retraction)
# decomp_qdt.zs(a=a, b=b, Lx=0.1, n_x=21, Ly=0.1, n_y=21, theta=[theta], avgZs=True, retraction=retraction);  # Zs
# nmse_qdt_ps = lc.nmse_freq(decomp_qdt.p_s , field_nlr_zs.pres_s[0])
# nmse_qdt_uzs = lc.nmse_freq(decomp_qdt.uz_s , field_nlr_zs.uz_s[0])
#%%
decomp_mono = Decomposition_2M(p_mtx=field_nlr.pres_s[0], controls=field_nlr.controls,
                               receivers=field_nlr.receivers, source_coord=source_coord)
decomp_mono.pk_tikhonov(plot_l=False, method='Tikhonov')
decomp_mono.zs(Lx=0.1, n_x=21, Ly=0.1, n_y=21, theta=[theta], avgZs=True); # Zs

#%%
plt.figure()
plt.semilogx(field_nlr.controls.freq, field_nlr.material.alpha, label = 'Miki')
plt.semilogx(field_nlr_zs.controls.freq, alpha_nlr, label = 'NLR')
plt.semilogx(decomp_qdt_v2.controls.freq, decomp_qdt_v2.alpha[0,:], label = 'QDT')
# plt.semilogx(decomp_qdt_v2.controls.freq, alpha_pk, label = 'QDT Q')
plt.semilogx(decomp_mono.controls.freq, decomp_mono.alpha[0,:], label = '2M')
plt.legend()
plt.xlabel('Frequency [Hz]')
plt.ylabel(r'$\alpha$ [-]')
plt.grid()
plt.ylim((-0.2, 1.2))
plt.tight_layout()

#%%
plt.figure()
plt.loglog(field_nlr.controls.freq, nmse_qdt_ps, label = 'NMSE ps QDT')
plt.loglog(field_nlr.controls.freq, nmse_qdt_uzs, label = 'NMSE uzs QDT')
plt.legend()
plt.xlabel('Frequency [Hz]')
plt.ylabel(r'$\alpha$ [-]')
plt.grid()
#plt.ylim((-0.2, 1.2))
plt.tight_layout()

#%%
theta_deg = np.array([0, 45])

mat = PorousAbsorber(air = field_nlr.air, controls = field_nlr.controls)
mat.miki(resistivity = field_nlr.material.resistivity)
alpha_pw_ref = np.zeros((len(theta_deg), len(field_nlr.controls.k0)))
for jthe in np.arange(len(theta_deg)):
    _,_,alpha_pw_ref[jthe,:] = mat.layer_over_rigid(thickness = field_nlr.material.thickness, 
                                 theta = np.deg2rad(theta_deg[jthe]))

plt.figure()
plt.semilogx(mat.freq, alpha_pw_ref[0,:], label = 'Miki - {}deg'.format(theta_deg[0]))
plt.semilogx(mat.freq, alpha_pw_ref[1,:], label = 'Miki - {}deg'.format(theta_deg[1]))
plt.semilogx(decomp_qdt_v2.controls.freq, decomp_qdt_v2.alpha[0,:], label = 'QDT')
#plt.semilogx(mat.freq, alpha_pw_ref[2,:], label = 'Miki - {}deg'.format(theta_deg[2]))
plt.legend()
plt.xlabel('Frequency [Hz]')
plt.ylabel(r'$\alpha$ [-]')
plt.grid()
plt.ylim((-0.2, 1.2))
#%% Compute plane wave ref.coef
weights_vec = np.zeros(len(decomp_qdt_v2.weights)+1)
weights_vec[0] = 1
weights_vec[1:] = decomp_qdt_v2.weights

Rp_s = np.zeros(len(field_nlr.controls.k0), dtype = 'complex')
Rp = np.zeros(len(field_nlr.controls.k0), dtype = 'complex')
for jf, k0 in enumerate(field_nlr.controls.k0):
    kernel_vec = np.zeros(len(decomp_qdt_v2.weights)+1, dtype = 'complex')
    kernel_vec[0] = 1
    kernel_vec[1:] = np.exp(-1j*decomp_qdt_v2.q * k0 * np.cos(np.deg2rad(theta_deg[0])))
    s_vec = decomp_qdt_v2.pk[1:, jf]/decomp_qdt_v2.pk[0,jf]
    fun_vec = s_vec * kernel_vec
    #print(np.linalg.norm(fun_vec-s_vec))
    Rp_s[jf] = np.dot(weights_vec, s_vec) #np.dot(weights_vec, fun_vec)
    Rp[jf] = np.dot(weights_vec, fun_vec)
alpha_pw_rec_s = 1 - (np.abs(Rp_s))**2
alpha_pw_rec = 1 - (np.abs(Rp))**2
#%%
plt.figure()
plt.semilogx(mat.freq, alpha_pw_ref[0,:], label = 'Miki - {}deg'.format(theta_deg[0]))
plt.semilogx(mat.freq, alpha_pw_ref[1,:], label = 'Miki - {}deg'.format(theta_deg[1]))
plt.semilogx(field_nlr.controls.freq, alpha_pw_rec_s, label = '{}deg'.format(theta_deg[0]))
plt.semilogx(field_nlr.controls.freq, alpha_pw_rec, label = '{}deg'.format(theta_deg[0]))
plt.semilogx(field_nlr.controls.freq, 1 - (np.abs(decomp_qdt_v2.pk[1, :]/decomp_qdt_v2.pk[0,:]))**2, label = '{}deg'.format(theta_deg[0]))
plt.legend()
plt.xlabel('Frequency [Hz]')
plt.ylabel(r'$\alpha$ [-]')
plt.grid()
plt.ylim((-0.2, 1.2))
plt.tight_layout()

#%%
Vp_recon_surf = decomp_qdt_v2.p_recon_ref/decomp_qdt_v2.p_recon_inc
Vp_recon_mean = np.mean(Vp_recon_surf, axis = 0)
alpha_recon_mean = 1 - (np.abs(Vp_recon_mean))**2

plt.figure()
plt.semilogx(mat.freq, alpha_pw_ref[0,:], label = 'Miki - {}deg'.format(theta_deg[0]))
plt.semilogx(mat.freq, alpha_pw_ref[1,:], label = 'Miki - {}deg'.format(theta_deg[1]))
plt.semilogx(field_nlr.controls.freq, alpha_recon_mean, label = 'Vp_recon')
plt.legend()
plt.xlabel('Frequency [Hz]')
plt.ylabel(r'$\alpha$ [-]')
plt.grid()
plt.ylim((-0.2, 1.2))
plt.tight_layout()


# # --------------------------------------------------------------------------------------------------------------------------
# #                                       Main folder from where to save/load simulations
# # -------------------------------------------------------------------------------------------------------------------------

# main_folder = 'D:/Work/UFSM/Disciplinas/Problemas Inversos/10_DCISM/saved_fields/'#'C:\\Workspace\\Python\\Decomp_Quad\\saved_fields\\'

# # --------------------------------------------------------------------------------------------------------------------------
# #                                                  Configure test parameters
# # --------------------------------------------------------------------------------------------------------------------------

# r = 0.2  # Source distance to sample's center [m]
# Lx = 0.6  # Sample size (x) [m]
# Ly = 0.6  # Sample size (y) [m]
# d = 0.1  # Sample thickness [m]
# resist = 5000  # Flow resistivity [Ns/m^4]

# ng = 199  # Quadrature order
# p_hz = 0.02  # Reconstruction point height
# retraction = 0.01  # Source's retraction for the reconstruction
# a = 0  # lower bond of the integral
# b = 70  # upper bond of the integral

# # Array parameters
# n_x = 8  # number of receivers in the x direction
# n_y = 8  # number of receivers in the y direction
# x_len = 0.4  # length of the x direction
# y_len = 0.4  # length of the y direction

# # Elevation and azimuth [deg] -> [rad]
# theta_d = 0
# theta = np.deg2rad(theta_d)
# phi_d = 0
# phi = np.deg2rad(phi_d)

# # --------------------------------------------------------------------------------------------------------------------------
# #                                                     AIR AND CONTROLS
# # --------------------------------------------------------------------------------------------------------------------------

# air = AirProperties(c0=343.0, rho0=1.21)
# controls = AlgControls(c0=air.c0, freq_vec=[100, 125, 160, 200, 250, 315, 400, 500,
#                                             630, 800, 1000, 1250, 1600, 2000])

# # --------------------------------------------------------------------------------------------------------------------------
# #                                                         MATERIAL
# # --------------------------------------------------------------------------------------------------------------------------

# material = PorousAbsorber(air=air, controls=controls)
# material.miki(resistivity=resist)
# material.layer_over_rigid(thickness=d, theta=theta)
# # material.plot_absorption()

# # --------------------------------------------------------------------------------------------------------------------------
# #                                                          SOURCE
# # --------------------------------------------------------------------------------------------------------------------------

# source_coord = sph2cart(r, np.pi / 2 - theta, phi)  # Source coordinates
# source = Source(coord=source_coord)

# # --------------------------------------------------------------------------------------------------------------------------
# #                                                         RECEIVERS
# # --------------------------------------------------------------------------------------------------------------------------

# bfs_array = 'r3d'  # Array codename for BEM field
# nlr_array = 'r3d'  # Array codename for NLR calculation
# n_mics = 128  # Number of mics
# dimension = '30x30x5'  # Array dimension

# receivers_r3d = Receiver()  # small volume random array
# receivers_r3d.random_3d_array(x_len=0.3, y_len=0.3, z_len=0.05, zr=0.013, n_total=n_mics)

# receivers_5mic = Receiver()  # used for reconstruction tests
# receivers_5mic.line_array(startat=0.00, line_len=p_hz+0.02, n_rec=5, direction='z')


# # --------------------------------------------------------------------------------------------------------------------------
# #                                                     BEM FIELD - ARRAY
# # --------------------------------------------------------------------------------------------------------------------------

# # field = BEMFlushSq(air=air, controls=controls, material=material, sources=source, receivers=receivers_r3d, n_gauss=36)
# # field.generate_mesh(Lx=Lx, Ly=Ly, Nel_per_wavelenth=6)  # Generate a mesh
# # field.plotly_scene(renderer="browser")  # Plot and see the scene
# # field.psurf()
# # field.p_fps()
# # field.uz_fps()
# f_name_field = f'bfs_Lx{int(Lx*100)}cm_Ly{int(Ly*100)}cm_d{int(d*100)}cm_el{int(theta_d)}d_az{int(phi_d)}d_r' \
#                f'{int(r*100)}cm_resist{int(resist/1e3)}k_{bfs_array}_{n_mics}_mics_{dimension}cm'
# # field.save(path=main_folder, filename=f_name_field)

# # --------------------------------------------------------------------------------------------------------------------------
# #                                                      FIELD - 5 MIC
# # --------------------------------------------------------------------------------------------------------------------------

# # field_5mic = BEMFlushSq(air=air, controls=controls, material=material, sources=source, receivers=receivers_5mic,
# #                         n_gauss=36)
# # field_5mic.generate_mesh(Lx=Lx, Ly=Ly, Nel_per_wavelenth=6)  # Generate a mesh
# # field_5mic.plotly_scene(renderer="browser")  # Plot and see the scene
# # field_5mic.psurf()
# # field_5mic.p_fps()
# # field_5mic.uz_fps()
# f_name_field_5mic = f'BEM_5mic_Lx{int(Lx*100)}cm_Ly{int(Ly*100)}cm_d{int(d * 100)}cm_el{int(theta_d)}d_az{int(phi_d)}' \
#                     f'd_r{int(r*100)}cm_resist{int(resist/1e3)}k'
# # field_5mic.save(path=main_folder, filename=f_name_field_5mic)

# # --------------------------------------------------------------------------------------------------------------------------
# #                                                      FIELD - NLR
# # --------------------------------------------------------------------------------------------------------------------------

# # field_nlr = NLRInfSph(air=air, controls=controls, material=material, sources=source, receivers=receivers_r3d)
# # field_nlr.p_nlr()
# # field_nlr.uz_nlr()
# f_name_field_nlr = f'NLR_Inf_d{int(d*100)}cm_el{int(theta_d)}d_az{int(phi_d)}d_r'\
#                    f'{int(r*100)}cm_resist{int(resist/1e3)}k_{nlr_array}_{n_mics}_mics_{dimension}cm'
# # field_nlr.save(path=main_folder, filename=f_name_field_nlr)

# # --------------------------------------------------------------------------------------------------------------------------
# #                                                  FIELD - NLR ZS - Z = 0
# # --------------------------------------------------------------------------------------------------------------------------

# receivers_nlr_zs = Receiver()
# receivers_nlr_zs.planar_array(x_len=0.1, n_x=21, y_len=0.1, n_y=21, zr=0.0)

# # field_nlr_zs = NLRInfSph(air=air, controls=controls, material=material, sources=source, receivers=receivers_nlr_zs)
# # field_nlr_zs.p_nlr()
# # field_nlr_zs.uz_nlr()
# f_name_field_nlr_zs = f'NLR_zs_d{int(d*100)}cm_el{int(theta_d)}d_az{int(phi_d)}d_r'\
#                    f'{int(r*100)}cm_resist{int(resist/1e3)}k_z__0'
# # field_nlr_zs.save(path=main_folder, filename=f_name_field_nlr_zs)

# # --------------------------------------------------------------------------------------------------------------------------
# #                                                        Loading BEM
# # --------------------------------------------------------------------------------------------------------------------------

# field = BEMFlushSq()  # array field
# field.load(path=main_folder, filename=f_name_field)
# # field.plotly_scene(renderer="browser")  # Plot and see the scene

# field_5mic = BEMFlushSq()  # reconstruction test field
# field_5mic.load(path=main_folder, filename=f_name_field_5mic)

# field_nlr = NLRInfSph()  # NLR field
# field_nlr.load(path=main_folder, filename=f_name_field_nlr)

# field_nlr_zs = NLRInfSph()  # NLR field
# field_nlr_zs.load(path=main_folder, filename=f_name_field_nlr_zs)  # Surface impedance NLR

# # --------------------------------------------------------------------------------------------------------------------------
# #                                                           NOISE
# # --------------------------------------------------------------------------------------------------------------------------

# snr = 40
# field.add_noise(snr=snr, uncorr=False)  # array
# field_5mic.add_noise(snr=snr, uncorr=False)  # 5 mics
# field_nlr.add_noise(snr=snr, uncorr=False)  # array - NLR
# field_nlr_zs.add_noise(snr=snr, uncorr=False)  # array - NLR ZS

# # --------------------------------------------------------------------------------------------------------------------------
# #                                              DECOMPOSITION - QUADRATURE
# # --------------------------------------------------------------------------------------------------------------------------

# decomp_qdt = Decomposition_QDT(p_mtx=field.pres_s[0], controls=controls, material=material,
#                                receivers=receivers_r3d, source_coord=source_coord, quad_order=ng)

# decomp_qdt.pk_tikhonov(plot_l=False, method='Ridge', a=a, b=b, retraction=retraction)


# # --------------------------------------------------------------------------------------------------------------------------
# #                                              DECOMPOSITION - PLANE WAVES
# # --------------------------------------------------------------------------------------------------------------------------

# decomp_plane = Decomposition_2W(p_mtx=field.pres_s[0], controls=controls, material=material,
#                                 receivers=receivers_r3d, source_coord=source_coord)
# decomp_plane.pk_tikhonov(plot_l=False, method='Ridge')

# # --------------------------------------------------------------------------------------------------------------------------
# #                                               DECOMPOSITION - MONOPOLES
# # --------------------------------------------------------------------------------------------------------------------------

# decomp_mono = Decomposition_2M(p_mtx=field.pres_s[0], controls=controls, material=material,
#                                receivers=receivers_r3d, source_coord=source_coord)
# decomp_mono.pk_tikhonov(plot_l=False, method='Ridge')

# # --------------------------------------------------------------------------------------------------------------------------
# #                                     PRESSURE, PARTICLE VELOCITY AND SURFACE IMPEDANCE
# # --------------------------------------------------------------------------------------------------------------------------

# # Quadrature
# decomp_qdt.reconstruct_p(a=a, b=b, hz=p_hz, n_pts=1, retraction=retraction)   # Quadrature
# decomp_qdt.reconstruct_uz(a=a, b=b, hz=p_hz, n_pts=1, retraction=retraction)  # Particle velocity
# decomp_qdt.zs(a=a, b=b, Lx=0.1, n_x=21, Ly=0.1, n_y=21, theta=[theta], avgZs=True, retraction=retraction)  # Zs

# # Plane Waves
# decomp_plane.reconstruct_pu(hz=p_hz, n_pts=1)  # Pressure and particle velocity
# decomp_plane.zs(Lx=0.1, n_x=21, Ly=0.1, n_y=21, theta=[theta], avgZs=True)  # Zs

# # Monopoles
# decomp_mono.reconstruct_pu(hz=p_hz, n_pts=1)  # Pressure and particle velocity
# decomp_mono.zs(Ly=0.1, n_y=21, theta=[theta], avgZs=True)  # Zs

# # NLR
# zs_nlr = field_nlr_zs.pres_s[0]/field_nlr_zs.uz_s[0]  # Zs = P / U

# zs_nlr_mean = np.zeros(len(controls.freq), dtype=complex)  # Mean of Zs
# alpha_nlr = np.zeros(len(controls.freq))  # Absorption coefficient

# for jf in range(0, len(controls.freq)):  # Absorption coefficient estimation
#     zs_nlr_mean[jf] = np.mean(zs_nlr[:, jf])
#     alpha_nlr[jf] = 1 - (np.abs(np.divide((zs_nlr_mean[jf] * np.cos(theta) - 1),
#                                           (zs_nlr_mean[jf] * np.cos(theta) + 1)))) ** 2

# # --------------------------------------------------------------------------------------------------------------------------
# #                                              PLOT - RECONSTRUCTED PRESSURE
# # --------------------------------------------------------------------------------------------------------------------------

# # plt.figure()
# plt.title('Reconstructed pressure at $p = %.2f$ [m]' % p_hz)  # field_5mic.pres_s[00][2, :] = pressure at 0.02m (p_hz)
# plt.semilogx(controls.freq, np.real(field_5mic.pres_s[00][2, :]), '-', color='black',  # BEM (Real)
#              linewidth=3, markersize=5, alpha=0.8, label=r'$P_\mathrm{Real}$ (BEM)')
# plt.semilogx(controls.freq, np.imag(field_5mic.pres_s[00][2, :]), '-', color='gray',  # BEM (Imaginary)
#              linewidth=3, markersize=5, alpha=0.8, label=r'$P_\mathrm{Imag}$ (BEM)')
# plt.semilogx(controls.freq, np.real(decomp_qdt.p_recon[0, :]), '-.', color='darkviolet',  # QUAD (Real)
#              linewidth=2, markersize=5, alpha=0.8, label=r'$P_\mathrm{re, Real}$ (QBDCIM)')
# plt.semilogx(controls.freq, np.imag(decomp_qdt.p_recon[0, :]), '-.', color='greenyellow',  # QUAD (Real)
#              linewidth=2, markersize=5, alpha=0.8, label=r'$P_\mathrm{re, Imag}$ (QBDCIM)')
# plt.semilogx(controls.freq, np.real(decomp_plane.p_recon[0, :]), '--', color='darkcyan',   # PLANE WAVE (Real)
#              linewidth=2, markersize=5, alpha=0.8, label=r'$P_\mathrm{re, Real}$ (2PW)')
# plt.semilogx(controls.freq, np.imag(decomp_plane.p_recon[0, :]), '--', color='orange',  # PLANE WAVE (Imag)
#              linewidth=2, markersize=5, alpha=0.8, label=r'$P_\mathrm{re, Imag}$ (2PW)')
# plt.semilogx(controls.freq, np.real(decomp_mono.p_recon[0, :]), '-o', color='blue',  # MONOPOLE (Real)
#              linewidth=2, markersize=5, alpha=0.8, label=r'$P_\mathrm{re, Real}$ (2MP)')
# plt.semilogx(controls.freq, np.imag(decomp_mono.p_recon[0, :]), '-o', color='red',  # MONOPOLE (Imag)
#              linewidth=2, markersize=5, alpha=0.8, label=r'$P_\mathrm{re, Imag}$ (2MP)')
# # plt.semilogx(controls.freq, np.real(field_nlr.pres_s[00][0, :]), '-o', color='teal',  # NLR (Real)
# #              linewidth=3, markersize=5, alpha=0.8, label=r'$P_\mathrm{NLR, Real}$')
# # plt.semilogx(controls.freq, np.imag(field_nlr.pres_s[00][0, :]), '-o', color='coral',  # NLR (Imag)
# #              linewidth=3, markersize=5, alpha=0.8, label=r'$P_\mathrm{NLR, Imag}$')
# plt.xlabel('Frequency [Hz]')
# plt.ylabel('Amplitude [Pa]')
# plt.xticks([50, 100, 500, 1000, 2000, 4000, 8000, 10000],
#            ['50', '100', '500', '1k', '2k', '4k', '8k', '10k'])
# plt.grid(linestyle='--', which='both')
# plt.legend(loc='best')
# plt.xlim((80, 3000))
# # plt.show()

# # --------------------------------------------------------------------------------------------------------------------------
# #                                               PLOT - RECONSTRUCTED UZ
# # --------------------------------------------------------------------------------------------------------------------------

# plt.figure()
# plt.title('Reconstructed particle velocity at $p = %.2f$ [m] ' % p_hz + '($n_g = %.0f$,' % ng +
#           ' $a = %i$ ' % a + 'e $b = %i$)' % b)
# plt.semilogx(controls.freq, np.real(field_5mic.uz_s[00][2, :]), '-', color='black',  # BEM (Real)
#              linewidth=3, markersize=5, alpha=0.8, label=r'$U_\mathrm{z, Real}$ (BEM)')
# plt.semilogx(controls.freq, np.imag(field_5mic.uz_s[00][2, :]), '-', color='gray',  # BEM (Imag)
#              linewidth=3, markersize=5, alpha=0.8, label=r'$U_\mathrm{z, Imag}$ (BEM)')
# plt.semilogx(controls.freq, np.real(decomp_qdt.uz_recon[0, :]), '-.', color='darkviolet',  # QUAD Reconstructed (Real)
#              linewidth=2, markersize=5, alpha=0.8, label=r'$U_\mathrm{z, re, Real}$ (QBDCIM)')
# plt.semilogx(controls.freq, np.imag(decomp_qdt.uz_recon[0, :]), '-.', color='greenyellow',  # QUAD Reconstructed  (Imag)
#              linewidth=2, markersize=5, alpha=0.8, label=r'$U_\mathrm{z,re, Imag}$ (QBDCIM)')
# plt.semilogx(controls.freq, np.real(decomp_plane.uz_recon[0, :]), '--', color='darkcyan',
#              linewidth=2, markersize=5, alpha=0.8, label=r'$U_\mathrm{z, re, Real}$ (2PW)')  # PLANE WAVE
# plt.semilogx(controls.freq, np.imag(decomp_plane.uz_recon[0, :]), '--', color='orange',
#              linewidth=2, markersize=5, alpha=0.8, label=r'$U_\mathrm{z,re, Imag}$ (2PW)')  # PLANE WAVE
# plt.semilogx(controls.freq, np.real(decomp_mono.uz_recon[0, :]), '-o', color='blue',
#              linewidth=2, markersize=5, alpha=0.8, label=r'$U_\mathrm{z, re, Real}$ (2MP)')  # MONOPOLE (Real)
# plt.semilogx(controls.freq, np.imag(decomp_mono.uz_recon[0, :]), '-o', color='red',
#              linewidth=2, markersize=5, alpha=0.8, label=r'$U_\mathrm{z,re, Imag}$ (2MP)')  # MONOPOLE (Imag)
# plt.semilogx(controls.freq, np.real(field_nlr.uz_s[00][0, :]), '--', color='black',  # NLR (Real)
#              linewidth=3, markersize=5, alpha=0.8, label=r'$U_\mathrm{z, Real}$ (NLR)')
# plt.semilogx(controls.freq, np.imag(field_nlr.uz_s[00][0, :]), '--', color='gray',  # NLR (Imag)
#              linewidth=3, markersize=5, alpha=0.8, label=r'$U_\mathrm{z, Imag}$ (NLR)')
# plt.xlabel('Frequency [Hz]')
# plt.ylabel('Amplitude [m/s]')
# plt.xticks([50, 100, 500, 1000, 2000, 4000, 8000, 10000],
#            ['50', '100', '500', '1k', '2k', '4k', '8k', '10k'])
# plt.grid(linestyle='--', which='both')
# plt.legend(loc='best', ncol=2)
# plt.xlim((80, 3000))
# # plt.show()

# # --------------------------------------------------------------------------------------------------------------------------
# #                                             INCIDENT PARTICLE VELOCITY
# # --------------------------------------------------------------------------------------------------------------------------

# plt.figure()
# plt.title('Incident particle velocity at $p = %.2f$ [m]' % p_hz)
# plt.semilogx(controls.freq, np.real(decomp_mono.uz_recon_inc[0, :]), '-o', color='blue',  # MONOPOLE (Real)
#              linewidth=2, markersize=5, alpha=0.8, label=r'$U_\mathrm{z, re, Real}$ (2MP)')
# plt.semilogx(controls.freq, np.imag(decomp_mono.uz_recon_inc[0, :]), '-o', color='red',  # MONOPOLE (Imag)
#              linewidth=2, markersize=5, alpha=0.8, label=r'$U_\mathrm{z,re, Imag}$ (2MP)')
# plt.semilogx(controls.freq, np.real(decomp_qdt.uz_recon_inc[0, :]), '-.', color='darkviolet',  # QUAD (Real)
#              linewidth=2, markersize=5, alpha=0.8, label=r'$U_\mathrm{z, re, Real}$ (QBDCIM)')
# plt.semilogx(controls.freq, np.imag(decomp_qdt.uz_recon_inc[0, :]), '-.', color='greenyellow',  # QUAD  (Imag)
#              linewidth=2, markersize=5, alpha=0.8, label=r'$U_\mathrm{z,re, Imag}$ (QBDCIM)')
# plt.semilogx(controls.freq, np.real(decomp_plane.uz_recon_inc[0, :]), '--', color='darkcyan',  # PLANE WAVE (Real)
#              linewidth=2, markersize=5, alpha=0.8, label=r'$U_\mathrm{z, re, Real}$ (2PW)')
# plt.semilogx(controls.freq, np.imag(decomp_plane.uz_recon_inc[0, :]), '--', color='orange',  # PLANE WAVE (Imag)
#              linewidth=2, markersize=5, alpha=0.8, label=r'$U_\mathrm{z,re, Imag}$ (2PW)')
# plt.xlabel('Frequency [Hz]')
# plt.ylabel('Amplitude [m/s]')
# plt.xticks([50, 100, 500, 1000, 2000, 4000, 8000, 10000],
#            ['50', '100', '500', '1k', '2k', '4k', '8k', '10k'])
# plt.grid(linestyle='--', which='both')
# plt.legend(loc='best', ncols=3)
# plt.xlim((80, 3000))
# plt.show()

# # --------------------------------------------------------------------------------------------------------------------------
# #                                             REFLECTED PARTICLE VELOCITY
# # --------------------------------------------------------------------------------------------------------------------------

# plt.figure()
# plt.title('Reflected particle velocity at $p = %.2f$ [m]' % p_hz)
# plt.semilogx(controls.freq, np.real(decomp_mono.uz_recon_ref[0, :]), '-o', color='blue',  # MONOPOLE (Real)
#              linewidth=2, markersize=5, alpha=0.8, label=r'$U_\mathrm{z, re, Real}$ (2MP)')
# plt.semilogx(controls.freq, np.imag(decomp_mono.uz_recon_ref[0, :]), '-o', color='red',  # MONOPOLE (Imag)
#              linewidth=2, markersize=5, alpha=0.8, label=r'$U_\mathrm{z,re, Imag}$ (2MP)')
# plt.semilogx(controls.freq, np.real(decomp_qdt.uz_recon_ref[0, :]), '-.', color='darkviolet',  # QUAD (Real)
#              linewidth=2, markersize=5, alpha=0.8, label=r'$U_\mathrm{z, re, Real}$ (QBDCIM)')
# plt.semilogx(controls.freq, np.imag(decomp_qdt.uz_recon_ref[0, :]), '-.', color='greenyellow',  # QUAD  (Imag)
#              linewidth=2, markersize=5, alpha=0.8, label=r'$U_\mathrm{z,re, Imag}$ (QBDCIM)')
# plt.xlabel('Frequency [Hz]')
# plt.ylabel('Amplitude [m/s]')
# plt.xticks([50, 100, 500, 1000, 2000, 4000, 8000, 10000],
#            ['50', '100', '500', '1k', '2k', '4k', '8k', '10k'])
# plt.grid(linestyle='--', which='both')
# plt.legend(loc='best')
# plt.xlim((80, 3000))
# # plt.show()

# # --------------------------------------------------------------------------------------------------------------------------
# #                                                   PLOT - ABSORPTION
# # --------------------------------------------------------------------------------------------------------------------------


# plt.figure(figsize=(7, 5))
# # plt.title('Absorption coefficient ($n_g = %.0f$,' % ng + ' $a = %i$ ' % a + 'e $b = %i$)' % b)
# plt.semilogx(material.freq, material.alpha, '--k', linewidth=3, alpha=1.0, label='Miki')
# plt.semilogx(decomp_qdt.material.freq, alpha_nlr, '-o', color='gold',
#              linewidth=2, markersize=5, alpha=0.8, label='NLR')
# plt.semilogx(decomp_mono.material.freq, decomp_mono.alpha[0], '-o', color='blue',
#              linewidth=2, markersize=5, alpha=0.8, label='2MP')
# plt.semilogx(decomp_plane.material.freq, decomp_plane.alpha[0], '-o', color='c',
#              linewidth=2, markersize=5, alpha=0.8, label='2PW')
# plt.semilogx(decomp_qdt.material.freq, decomp_qdt.alpha[0], '-o', color='violet',
#              linewidth=2, markersize=5, alpha=0.8, label='QBDCIM')
# plt.legend(loc='best')
# plt.xticks([50, 100, 500, 1000, 2000, 4000, 8000, 10000],
#            ['50', '100', '500', '1k', '2k', '4k', '8k', '10k'])
# plt.grid(linestyle='--', which='both')
# plt.xlabel('Frequência [Hz]')
# plt.ylabel(r'$\alpha$ [-]')
# plt.xlim((80, 3000))
# plt.ylim((-0.2, 1.2))
# plt.tight_layout()
# # plt.show()

# # --------------------------------------------------------------------------------------------------------------------------
# #                                                   PLOT - IMPEDANCE
# # --------------------------------------------------------------------------------------------------------------------------

# plt.figure(figsize=(7, 5))
# # plt.title('Reconstructed surface impedance')
# plt.semilogx(material.freq, np.real(material.Zs)/(343*1.21), '-', color='black',  linewidth=3,
#              alpha=1.0, label='Miki (Real)')
# plt.semilogx(material.freq, np.imag(material.Zs)/(343*1.21), '-', color='gray', linewidth=3,
#              alpha=1.0, label='Miki (Imag)')
# plt.semilogx(decomp_qdt.material.freq, np.real(zs_nlr_mean), '-', color='b',
#              linewidth=2, markersize=5, alpha=0.8, label='NLR (Real)')
# plt.semilogx(decomp_qdt.material.freq, np.imag(zs_nlr_mean), '-', color='r',
#              linewidth=2, markersize=5, alpha=0.8, label='NLR (Imag)')
# # plt.semilogx(decomp_plane.material.freq, np.real(decomp_plane.Zs), '--', color='green',
# #              linewidth=2, markersize=5, alpha=0.8, label='2PW Real')
# # plt.semilogx(decomp_plane.material.freq, np.imag(decomp_plane.Zs), '--', color='orange',
# #              linewidth=2, markersize=5, alpha=0.8, label='2PW Imag')
# plt.semilogx(decomp_qdt.material.freq, np.real(decomp_qdt.Zs), '--', color='b',
#              linewidth=2, markersize=5, alpha=0.8, label='QBDCIM (Real)')
# plt.semilogx(decomp_qdt.material.freq, np.imag(decomp_qdt.Zs), '--', color='r',
#              linewidth=2, markersize=5, alpha=0.8, label='QBDCIM (Imag)')

# # plt.semilogx(decomp_qdt.material.freq, np.real(decomp_qdt.Zs_grad), '--', color='b',
# #              linewidth=2, markersize=5, alpha=0.8, label='QDT Real (Gradient)')
# # plt.semilogx(decomp_qdt.material.freq, np.imag(decomp_qdt.Zs_grad), '--', color='r',
# #              linewidth=2, markersize=5, alpha=0.8, label='QDT Imag (Gradient)')
# plt.legend(loc='best')
# plt.xticks([50, 100, 500, 1000, 2000, 4000, 8000, 10000],
#            ['50', '100', '500', '1k', '2k', '4k', '8k', '10k'])
# plt.grid(linestyle='--', which='both')
# plt.xlabel('Frequência [Hz]')
# plt.ylabel(r'$\alpha$ [-]')
# plt.xlim((80, 3000))
# plt.tight_layout()
# plt.show()
