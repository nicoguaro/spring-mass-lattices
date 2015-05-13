# -*- coding: utf-8 -*-
"""
Compute the dispersion curves for a spring mass
lattice with hexagonal cells of side a.


@author: Nicolas Guarin-Zapata
@date: October 13, 2014
"""
from __future__ import division
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
            
            km = 0.5*(np.sqrt(3)*ky - kx)
            kp = 0.5*(np.sqrt(3)*ky + kx)
            skx = np.sin(kx)
            skp = np.sin(kp)
            skm = np.sin(km)
            ckx = np.cos(kx)
            ckp = np.cos(kp)
            ckm = np.cos(km)
            sq3 = np.sqrt(3)
            
            K = np.array([
                [4, 0, (1j*(8*skx -5*skp + 3*skm) - 8*ckx - 5*ckp - 3*ckm)/4,
                 0],
                [0, 2*sq3, 0,  
                 -(1j*(5*sq3*skp - 3**1.5*skm) + 5*sq3*ckp + 3**1.5*ckm)/4],
                [-(1j*(8*skx -5*skp + 3*skm) + 8*ckx + 5*ckp + 3*ckm)/4, 0,
                 4, 0],
                [0, (1j*(5*sq3*skp - 3**1.5*skm) -  5*sq3*ckp - 3**1.5*ckm)/4,
                 0, 2*sq3]
                ])
                
            M = np.array([[1., 0., 0., 0],
                          [0., 1., 0., 0],
                          [0., 0., 1., 0],
                          [0., 0., 0., 1]])
            vals = LA.eigvalsh(K, M)
            omega[:,i,j] = vals
        
    return omega
    
def freqs(kx_vec, ky_vec):
    m, n = kx_vec.shape
    omega = np.zeros((4, m, n))
    km = 0.5*(np.sqrt(3)*ky_vec - kx_vec)
    kp = 0.5*(np.sqrt(3)*ky_vec + kx_vec)
    skx = np.sin(kx_vec)
    skp = np.sin(kp)
    skm = np.sin(km)
    ckx = np.cos(kx_vec)
    ckp = np.cos(kp)
    ckm = np.cos(km)
    sq3 = np.sqrt(3)
    
    w1 = -(sq3*np.sqrt(25*skp**2 - 30*skm*skp+9*skm**2 + 25*ckp**2 + 
        30*ckm*ckp + 9*ckm**2)-8*sq3)/4
    w2 = (np.sqrt(64*skx**2 + (48*skm - 80*skp)*skx + 25*skp**2 - 30*skm*skp + 
        9*skm**2 + 64*ckx**2 + (80*ckp + 48*ckm)*ckx + 25*ckp**2 + 30*ckm*ckp
        + 9*ckm**2) - 16)/4
    w3 = (sq3*np.sqrt(25*skp**2 - 30*skm*skp + 9*skm**2 + 25*ckp**2 +
        30*ckm*ckp + 9*ckm**2) + 8*sq3)/4
    w4 = (np.sqrt(64*skx**2 + (48*skm - 80*skp)*skx + 25*skp**2 - 30*skm*skp + 
        9*skm**2 + 64*ckx**2 + (80*ckp + 48*ckm)*ckx + 25*ckp**2 + 30*ckm*ckp
        + 9*ckm**2) + 16)
        
    omega[0,:,:] = w1
    omega[1,:,:] = w2
    omega[2,:,:] = w3
    omega[3,:,:] = w4
    
    return omega
        
    
    

kx_vec, ky_vec = np.mgrid[-np.pi: np.pi:201j, -np.pi: np.pi: 201j]
omega = disp_rel(kx_vec, ky_vec)
omega = np.abs(omega)

#plt.close('all')
#plt.figure()
#plt.contourf(kx_vec, ky_vec, np.sqrt(omega[0,:,:] + omega[1,:,:]), cmap='hot')
#plt.colorbar()
#plt.contour(kx_vec, ky_vec, np.sqrt(omega[0,:,:] + omega[1,:,:]), colors='k')
#plt.axis('image')
#plt.xlabel(r"$k_x a/\pi$", size=18)
#plt.ylabel(r"$k_y a/\pi$", size=18)
#plt.figure()
#plt.contourf(kx_vec, ky_vec, np.sqrt(omega[2,:,:] + omega[3,:,:]), cmap='hot')
#plt.colorbar()
#plt.contour(kx_vec, ky_vec, np.sqrt(omega[2,:,:] + omega[3,:,:]), colors='k')
#plt.axis('image')
#plt.xlabel(r"$k_x a/\pi$", size=18)
#plt.ylabel(r"$k_y a/\pi$", size=18)
#plt.show()


omega2 = freqs(kx_vec, ky_vec)
omega2 = np.abs(omega2)
from mayavi import mlab
#mlab.mesh(kx_vec, ky_vec, np.sqrt(omega2[0,:,:]), colormap='hot')
#mlab.mesh(kx_vec, ky_vec, np.sqrt(omega2[1,:,:]), colormap='hot')
#mlab.mesh(kx_vec, ky_vec, np.sqrt(omega2[2,:,:]), colormap='hot')
#mlab.mesh(kx_vec, ky_vec, np.sqrt(omega2[3,:,:]), colormap='hot')
#mlab.mesh(kx_vec, ky_vec, np.sqrt(omega2[0,:,:] + omega2[1,:,:]) -
#    np.sqrt(omega[0,:,:] + omega[1,:,:]), colormap='hot')
#mlab.mesh(kx_vec, ky_vec, np.sqrt(omega[2,:,:] + omega[3,:,:]), colormap='cool')
#mlab.mesh(kx_vec, ky_vec, np.sqrt(omega2[0,:,:] + omega2[1,:,:]), colormap='hot')
#mlab.mesh(kx_vec, ky_vec, np.sqrt(omega2[2,:,:] + omega2[3,:,:]), colormap='hot')
#mlab.axes()
mlab.show()

