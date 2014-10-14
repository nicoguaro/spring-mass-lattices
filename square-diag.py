# -*- coding: utf-8 -*-
"""
Compute the dispersion curves for a spring mass
lattice with square cells of side $a$ and $b$
and spring  cosntants $c_1$, $c_2$, and $c_d$
for horizontal, vertical and diagonal springs.

@author: Nicolas Guarin-Zapata
@date: October 14, 2014
"""
import numpy as np
import scipy.linalg as LA
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

def disp_rel(kx_vec, ky_vec, w1, w2, wd):
    m, n = kx_vec.shape
    omega = np.zeros((2, m, n))
    for i in range(m):
        for j in range(n):
            kx = kx_vec[i,j]
            ky = ky_vec[i,j]

            val1 = 8*w1**2*np.sin(0.5*kx)**2 + 4*wd**2*(np.sin(0.5*kx + 0.5*ky)**2 +
                                                        np.sin(0.5*kx - 0.5*ky)**2)

            val2 = 8*w2**2*np.sin(0.5*ky)**2 + 4*wd**2*(np.sin(0.5*kx + 0.5*ky)**2 +
                                                        np.sin(0.5*kx - 0.5*ky)**2)

            omega[:,i,j] = [val1, val2]
        
    return omega

m = 1.  # Mass
c1 = 1.  # Horizontal spring
c2 = 1.  # Vertical spring
cd =  np.sqrt(2)/2 # Diagonal spring
w1 = np.sqrt(c1/m)
w2 = np.sqrt(c2/m)
wd = np.sqrt(cd/m)

kx_vec, ky_vec = np.mgrid[-np.pi: np.pi: 201j, -np.pi: np.pi: 201j]
omega = disp_rel(kx_vec, ky_vec, w1, w2, wd)
omega = np.abs(omega)

plt.close('all')
plt.figure()
plt.contourf(kx_vec, ky_vec, np.sqrt(omega[0,:,:] + omega[1,:,:]), cmap='hot')
plt.colorbar()
plt.contour(kx_vec, ky_vec, np.sqrt(omega[0,:,:] + omega[1,:,:]), colors='k')
plt.axis('image')
plt.xlabel(r"$k_x a/\pi$", size=18)
plt.ylabel(r"$k_y a/\pi$", size=18)
plt.savefig("Notes/img/square-diag-cd=%g.pdf"%cd, bbox="tight")
plt.show()
