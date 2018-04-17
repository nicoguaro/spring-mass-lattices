# -*- coding: utf-8 -*-
"""
Compute the dispersion curves for a spring mass
lattice with square cells of sided $a$, $b$ and
spring constant $c$.

@author: Nicolas Guarin-Zapata
"""
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 8

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

kx_vec, ky_vec = np.mgrid[-pi: pi: 201j, -pi: pi: 201j]
omega = disp_rel(kx_vec, ky_vec)
omega = np.abs(omega)
omega = np.sqrt(omega[0,:,:] + omega[1,:,:])

plt.close('all')
plt.figure(figsize=(4, 4))
plt.contourf(kx_vec/pi, ky_vec/pi, omega, cmap='hot')
plt.colorbar()
plt.contour(kx_vec/pi, ky_vec/pi, omega, colors='k')
plt.axis('image')
plt.xlabel(r"$k_x a/\pi$")
plt.ylabel(r"$k_y a/\pi$")
plt.savefig("Notes/img/square-disp.pdf", bbox_inches="tight")


#%% Plot in the contour of the irreducible Brillouin zone
plt.figure(figsize=(4, 2.2))
k = np.sqrt(kx_vec**2 + ky_vec**2)
xaxis = range(0,301)
# Vertical lines
ymax = np.max(omega)
plt.plot([0,0], [0,ymax], 'gray')
plt.plot([100,100], [0,ymax], 'gray')
plt.plot([200,200], [0,ymax], 'gray')
# Parts of the plot
plt.plot(xaxis[0:101], omega[100::,100], 'k')
plt.plot(xaxis[100:201], omega[-1,100::], 'k')
plt.plot(xaxis[200:301], omega.diagonal()[0:101], 'k')
plt.xticks([0, 100, 200, 300], [r"$\Gamma$", r"$X$", r"$M$", r"$\Gamma$"])
plt.ylabel(r"$\Omega$")
plt.savefig("Notes/img/square-disp-irreducible.pdf", bbox_inches="tight")

plt.show()
