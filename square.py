# -*- coding: utf-8 -*-
"""
Compute the dispersion curves for a spring mass
lattice with square cells of sided $a$, $b$ and
spring constant $c$.

@author: Nicolas Guarin-Zapata
@date: October 14, 2014
"""
import numpy as np
import scipy.linalg as LA
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

def disp_rel(kx_vec, ky_vec):
    m, n = kx_vec.shape
    omega = np.zeros((2, m, n))
    for i in range(m):
        for j in range(n):
            kx = kx_vec[i,j]
            ky = ky_vec[i,j]

            val1 = 8*np.sin(0.5*kx)**2
            val2 = 8*np.sin(0.5*ky)**2
            omega[:,i,j] = [val1, val2]
        
    return omega

kx_vec, ky_vec = np.mgrid[-np.pi: np.pi: 201j, -np.pi: np.pi: 201j]
omega = disp_rel(kx_vec, ky_vec)
omega = np.abs(omega)

plt.close('all')
plt.figure()
plt.contourf(kx_vec, ky_vec, np.sqrt(omega[0,:,:] + omega[1,:,:]), cmap='hot')
plt.colorbar()
plt.contour(kx_vec, ky_vec, np.sqrt(omega[0,:,:] + omega[1,:,:]), colors='k')
plt.axis('image')
plt.xlabel(r"$k_x a/\pi$", size=18)
plt.ylabel(r"$k_y a/\pi$", size=18)
plt.savefig("Notes/img/square-disp.pdf", bbox="tight")
plt.show()
