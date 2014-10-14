# -*- coding: utf-8 -*-
"""
Compute the dispersion curves for a spring mass
lattice with hexagonal cells of side a.


@author: Nicolas Guarin-Zapata
@date: October 13, 2014
"""
import numpy as np
import scipy.linalg as LA
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

def disp_rel(kx_vec, ky_vec):
    m, n = kx_vec.shape
    omega = np.zeros((4, m, n))
    for i in range(m):
        for j in range(n):
            kx = kx_vec[i,j]
            ky = ky_vec[i,j]
            
            K = np.array([
                [-3, 0, (1j*(5*np.sin((np.sqrt(3)*ky + kx)/2) - 
                3*np.sin((np.sqrt(3)*ky - kx)/2)-4*np.sin(kx)) + 
                5*np.cos((np.sqrt(3)*ky + kx)/2) + 
                3*np.cos((np.sqrt(3)*ky - kx)/2) + 4*np.cos(kx))/4, 0],
                [0, -2*np.sqrt(3), 0, 
                 (1j*(5*np.sqrt(3)*np.sin((np.sqrt(3)*ky + kx)/2) -
                 3**1.5*np.sin((np.sqrt(3)*ky-kx)/2)) + 
                 5*np.sqrt(3)*np.cos((np.sqrt(3)*ky + kx)/2) + 
                 3**1.5*np.cos((np.sqrt(3)*ky-kx)/2))/4],
                [-(1j*(5*np.sin((np.sqrt(3)*ky + kx)/2) -
                3*np.sin((np.sqrt(3)*ky - kx)/2)-4*np.sin(kx)) -
                5*np.cos((np.sqrt(3)*ky + kx)/2)-
                3*np.cos((np.sqrt(3)*ky - kx)/2)-4*np.cos(kx))/4, 0, -3, 0],
                [0, -(1j*(5*np.sqrt(3)*np.sin((np.sqrt(3)*ky + kx)/2)-
                3**1.5*np.sin((np.sqrt(3)*ky - kx)/2))-
                5*np.sqrt(3)*np.cos((np.sqrt(3)*ky + kx)/2)-
                3**1.5*np.cos((np.sqrt(3)*ky - kx)/2))/4, 0, -2*np.sqrt(3)]])
                
            M = -np.array([[1., 0., 0., 0],
                          [0., 1., 0., 0],
                          [0., 0., 1., 0],
                          [0., 0., 0., 1]])
            vals = LA.eigvalsh(-K, -M)
            omega[:,i,j] = vals
        
    return omega

kx_vec, ky_vec = np.mgrid[-np.pi: np.pi: 201j, -np.pi: np.pi: 201j]
omega = disp_rel(kx_vec, ky_vec)
omega = np.abs(omega)

plt.close('all')
plt.figure()
plt.contourf(kx_vec, ky_vec, np.sqrt(omega[0,:,:] + omega[1,:,:]), cmap='gray')
plt.axis('image')
plt.xlabel(r"$k_x a/\pi$", size=18)
plt.ylabel(r"$k_y a/\pi$", size=18)
plt.colorbar()
plt.figure()
plt.contourf(kx_vec, ky_vec, np.sqrt(omega[2,:,:] + omega[3,:,:]), cmap='gray')
plt.axis('image')
plt.xlabel(r"$k_x a/\pi$", size=18)
plt.ylabel(r"$k_y a/\pi$", size=18)
plt.colorbar()
plt.show()
