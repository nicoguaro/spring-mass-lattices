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

def disp_rel(kx_vec, ky_vec, m1, m2, cd, cs):
    m, n = kx_vec.shape
    omega = np.zeros((4, m, n))
    for i in range(m):
        for j in range(n):
            kx = kx_vec[i,j]
            ky = ky_vec[i,j]
            
            val1 = -(np.sqrt((4*np.cos(kx)**2 - 8*np.cos(kx) + 4)*m2**2*cs**2 + ((8*cd - 8*cd*np.cos(kx))*m2**2 + (8*cd*np.cos(kx) - 8*cd)*m1*m2)*cs + 4*cd**2*
m2**2 + (cd**2*np.sin(ky + kx)**2 + (2*cd**2*np.sin(ky) + 2*cd**2*np.sin(kx))*np.sin(ky + kx) + cd**2*np.cos(ky + kx)**2 + 
(2*cd**2*np.cos(ky) + 2*cd**2*np.cos(kx) + 2*cd**2)*np.cos(ky + kx) + cd**2*np.sin(ky)**2 + 2*cd**2*np.sin(kx)*np.sin(ky) + cd**2*np.cos(ky)**2 + 
(2*cd**2*np.cos(kx) + 2*cd**2)*np.cos(ky) + cd**2*np.sin(kx)**2 + cd**2*np.cos(kx)**2 + 2*cd**2*np.cos(kx) - 7*cd**2)*m1*m2 + 4*cd**2*m1**2) + (2*np.cos(kx) - 2)*m2*
cs - 2*cd*m2 - 2*cd*m1)/(m1*m2)

            val2 = -(np.sqrt((4*np.cos(ky)**2 - 8*np.cos(ky) + 4)*m2**2*cs**2 + ((8*cd - 8*cd*np.cos(ky))*m2**2 + (8*cd*np.cos(ky) - 8*cd)*m1*m2)*cs + 4*cd**2*
m2**2 + (cd**2*np.sin(ky + kx)**2 + (2*cd**2*np.sin(ky) + 2*cd**2*np.sin(kx))*np.sin(ky + kx) + cd**2*np.cos(ky + kx)**2 + 
(2*cd**2*np.cos(ky) + 2*cd**2*np.cos(kx) + 2*cd**2)*np.cos(ky + kx) + cd**2*np.sin(ky)**2 + 2*cd**2*np.sin(kx)*np.sin(ky) + cd**2*np.cos(ky)**2 + 
(2*cd**2*np.cos(kx) + 2*cd**2)*np.cos(ky) + cd**2*np.sin(kx)**2 + cd**2*np.cos(kx)**2 + 2*cd**2*np.cos(kx) - 7*cd**2)*m1*m2 + 4*cd**2*m1**2) + 2*np.cos(ky)*m2*cs - 
2*m2*cs - 2*cd*m2 - 2*cd*m1)/(m1*m2)
            
            val3 = (np.sqrt((4*np.cos(kx)**2 - 8*np.cos(kx) + 4)*m2**2*cs**2 + ((8*cd - 8*cd*np.cos(kx))*m2**2 + (8*cd*np.cos(kx) - 8*cd)*m1*m2)*cs + 4*cd**2*
m2**2 + (cd**2*np.sin(ky + kx)**2 + (2*cd**2*np.sin(ky) + 2*cd**2*np.sin(kx))*np.sin(ky + kx) + cd**2*np.cos(ky + kx)**2 + 
(2*cd**2*np.cos(ky) + 2*cd**2*np.cos(kx) + 2*cd**2)*np.cos(ky + kx) + cd**2*np.sin(ky)**2 + 2*cd**2*np.sin(kx)*np.sin(ky) + cd**2*np.cos(ky)**2 + 
(2*cd**2*np.cos(kx) + 2*cd**2)*np.cos(ky) + cd**2*np.sin(kx)**2 + cd**2*np.cos(kx)**2 + 2*cd**2*np.cos(kx) - 7*cd**2)*m1*m2 + 4*cd**2*m1**2) + (2 - 2*np.cos(kx))*m2*
cs + 2*cd*m2 + 2*cd*m1)/(m1*m2)

            val4 = (np.sqrt((4*np.cos(ky)**2 - 8*np.cos(ky) + 4)*m2**2*cs**2 + ((8*cd - 8*cd*np.cos(ky))*m2**2 + (8*cd*np.cos(ky) - 8*cd)*m1*m2)*cs + 4*cd**2*
m2**2 + (cd**2*np.sin(ky + kx)**2 + (2*cd**2*np.sin(ky) + 2*cd**2*np.sin(kx))*np.sin(ky + kx) + cd**2*np.cos(ky + kx)**2 + 
(2*cd**2*np.cos(ky) + 2*cd**2*np.cos(kx) + 2*cd**2)*np.cos(ky + kx) + cd**2*np.sin(ky)**2 + 2*cd**2*np.sin(kx)*np.sin(ky) + cd**2*np.cos(ky)**2 + 
(2*cd**2*np.cos(kx) + 2*cd**2)*np.cos(ky) + cd**2*np.sin(kx)**2 + cd**2*np.cos(kx)**2 + 2*cd**2*np.cos(kx) - 7*cd**2)*m1*m2 + 4*cd**2*m1**2) + (2 - 2*np.cos(ky))*m2*
cs + 2*cd*m2 + 2*cd*m1)/(m1*m2)
            
            omega[:,i,j] = [val1, val2, val3, val4]
        
    return omega

m1 = 1.
m2 = 1.
cs = 1.
cd = 10.
kx_vec, ky_vec = np.mgrid[-np.pi: np.pi: 201j, -np.pi: np.pi: 201j]
omega = disp_rel(kx_vec, ky_vec, m1, m2, cd, cs)
omega = np.abs(omega)

plt.close('all')
# First pair of modes
plt.figure()
plt.contourf(kx_vec, ky_vec, np.sqrt(omega[0,:,:] + omega[1,:,:]), cmap='hot')
plt.colorbar()
plt.contour(kx_vec, ky_vec, np.sqrt(omega[0,:,:] + omega[1,:,:]), colors='k')
plt.axis('image')
plt.xlabel(r"$k_x a/\pi$", size=18)
plt.ylabel(r"$k_y a/\pi$", size=18)
plt.savefig("Notes/img/square-bcc-disp1-c=%g-m=%g.pdf"%(m2/m1, cd/cs), bbox="tight")
# Second pair of modes
plt.figure()
plt.contourf(kx_vec, ky_vec, np.sqrt(omega[2,:,:] + omega[3,:,:]), cmap='hot')
plt.colorbar()
plt.contour(kx_vec, ky_vec, np.sqrt(omega[2,:,:] + omega[3,:,:]), colors='k')
plt.axis('image')
plt.xlabel(r"$k_x a/\pi$", size=18)
plt.ylabel(r"$k_y a/\pi$", size=18)
plt.savefig("Notes/img/square-bcc-disp2-c=%g-m=%g.pdf"%(m2/m1, cd/cs), bbox="tight")
plt.show()
