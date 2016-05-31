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
cd = .1
kx_vec, ky_vec = np.mgrid[-np.pi: np.pi: 201j, -np.pi: np.pi: 201j]
omega = disp_rel(kx_vec, ky_vec, m1, m2, cd, cs)
omega = np.abs(omega)
omega1 = np.sqrt(omega[0,:,:] + omega[1,:,:])
omega2 = np.sqrt(omega[2,:,:] + omega[3,:,:])

plt.close('all')
# First pair of modes
plt.figure(figsize=(3, 3))
plt.contourf(kx_vec/pi, ky_vec/pi, omega1, cmap='hot')
plt.colorbar()
plt.contour(kx_vec/pi, ky_vec/pi, omega1, colors='k')
plt.axis('image')
plt.xlabel(r"$k_x a/\pi$")
plt.ylabel(r"$k_y a/\pi$")
plt.savefig("Notes/img/square-bcc-disp1-c=%g-m=%g.pdf"%(cd/cs, m2/m1),
            bbox_inches="tight")
# Second pair of modes
plt.figure(figsize=(3, 3))
plt.contourf(kx_vec/pi, ky_vec/pi, omega2, cmap='hot')
plt.colorbar()
plt.contour(kx_vec/pi, ky_vec/pi, omega2, colors='k')
plt.axis('image')
plt.xlabel(r"$k_x a/\pi$")
plt.ylabel(r"$k_y a/\pi$")
#plt.savefig("Notes/img/square-bcc-disp2-c=%g-m=%g.pdf"%(cd/cs, m2/m1),
#            bbox_inches="tight")


#%% Plot in the contour of the irreducible Brillouin zone
plt.figure(figsize=(3, 1.9))
k = np.sqrt(kx_vec**2 + ky_vec**2)
xaxis = range(0,301)
# Vertical lines
ymax = 14
plt.plot([0,0], [0,ymax], 'gray')
plt.plot([100,100], [0,ymax], 'gray')
plt.plot([200,200], [0,ymax], 'gray')
# Parts of the plot
plt.plot(xaxis[0:101], omega1[100::,100], 'k')
plt.plot(xaxis[100:201], omega1[-1,100::], 'k')
plt.plot(xaxis[200:301], omega1.diagonal()[0:101], 'k')
plt.plot(xaxis[0:101], omega2[100::,100], 'k')
plt.plot(xaxis[100:201], omega2[-1,100::], 'k')
plt.plot(xaxis[200:301], omega2.diagonal()[0:101], 'k')
plt.xticks([0, 100, 200, 300], [r"$\Gamma$", r"$X$", r"$M$", r"$\Gamma$"])
plt.ylabel(r"$\Omega$")
plt.xlabel(r"$k$-space position")
plt.ylim(0, ymax)
#plt.savefig("Notes/img/square-bcc-irreducible-c=%g-m=%g.pdf"%(cd/cs, m2/m1),
#            bbox_inches="tight")
plt.show()
