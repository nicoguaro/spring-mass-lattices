# -*- coding: utf-8 -*-
"""
Compute the dispersion curves for a spring mass
lattice with hexagonal cells of side a.


@author: Nicolas Guarin-Zapata
"""
from __future__ import division
import numpy as np
from numpy import sqrt, cos, sin, pi
import scipy.linalg as LA
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 8

def disp_rel(kx_vec, ky_vec):
    m, n = kx_vec.shape
    omega = np.zeros((4, m, n))
    for i in range(m):
        for j in range(n):
            kx = kx_vec[i,j]
            ky = ky_vec[i,j]

            km = 0.5*(sqrt(3)*ky - kx)
            kp = 0.5*(sqrt(3)*ky + kx)
            skx = sin(kx)
            skp = sin(kp)
            skm = sin(km)
            ckx = cos(kx)
            ckp = cos(kp)
            ckm = cos(km)
            sq3 = sqrt(3)

            K = np.array([
                [4, 0, 2*1j*skx - 1j*skp + 1j*skm - 2*ckx - ckp - ckm, 0],
                [0, 2*sq3, 0, -sq3*1j*skp + sq3*1j*skm - sq3*ckp - sq3*ckm],
                [-2*1j*skx + 1j*skp-1j*skm - 2*ckx - ckp - ckm, 0, 4, 0],
                [0, sq3*1j*skp - sq3*1j*skm - sq3*ckp - sq3*ckm, 0, 2*sq3]])



            M = np.array([[1., 0., 0., 0],
                          [0., 1., 0., 0],
                          [0., 0., 1., 0],
                          [0., 0., 0., 1]])
            vals = LA.eigvalsh(K, M)
            omega[:,i,j] = vals

    return omega

def freqs(kx_vec, ky_vec):
    """Analytical solution for the hexagon"""
    m, n = kx_vec.shape
    omega = np.zeros((4, m, n))
    km = 0.5*(sqrt(3)*ky_vec - kx_vec)
    kp = 0.5*(sqrt(3)*ky_vec + kx_vec)
    skx = sin(kx_vec)
    skp = sin(kp)
    skm = sin(km)
    ckx = cos(kx_vec)
    ckp = cos(kp)
    ckm = cos(km)
    sq3 = sqrt(3)

    w1 = 2*sq3 - sq3*sqrt(skp**2 - 2*skm*skp + skm**2 + ckp**2 + 2*ckm*ckp +
        ckm**2)
    w2 = 2*sq3 + sq3*sqrt(skp**2 - 2*skm*skp + skm**2 + ckp**2 + 2*ckm*ckp +
        ckm**2)
    w3 = 4 - sqrt(4*skx**2 + (4*skm - 4*skp)*skx + skp**2 - 2*skm*skp +
        skm**2 + 4*ckx**2 + (4*ckp + 4*ckm)*ckx + ckp**2 + 2*ckm*ckp + ckm**2)
    w4 = 4 + sqrt(4*skx**2 + (4*skm - 4*skp)*skx + skp**2 - 2*skm*skp +
        skm**2 + 4*ckx**2 + (4*ckp + 4*ckm)*ckx + ckp**2 + 2*ckm*ckp + ckm**2)

    omega[0,:,:] = w1
    omega[1,:,:] = w3
    omega[2,:,:] = w2
    omega[3,:,:] = w4

    return omega




kx_vec, ky_vec = np.mgrid[-np.pi: np.pi:201j, -np.pi: np.pi: 201j]
omega = disp_rel(kx_vec, ky_vec)
omega = np.abs(omega)
omega_low = np.sqrt(omega[0,:,:] + omega[1,:,:])
omega_high = np.sqrt(omega[2,:,:] + omega[3,:,:])
omega_max = np.max(omega_high)

plt.close('all')
plt.figure(figsize=(3, 3))
plt.contourf(kx_vec/pi, ky_vec/pi, omega_low, cmap='hot')
plt.colorbar()
plt.contour(kx_vec/pi, ky_vec/pi, omega_low, colors='k')
plt.axis('image')
plt.xlabel(r"$k_x a/\pi$")
plt.ylabel(r"$k_y a/\pi$")
plt.savefig("Notes/img/hexagon-disp1.pdf",
            bbox_inches="tight")
plt.figure(figsize=(3, 3))
plt.contourf(kx_vec/pi, ky_vec/pi, omega_high, cmap='hot')
plt.colorbar()
plt.contour(kx_vec/pi, ky_vec/pi, omega_high, colors='k')
plt.axis('image')
plt.xlabel(r"$k_x a/\pi$")
plt.ylabel(r"$k_y a/\pi$")
plt.savefig("Notes/img/hexagon-disp2.pdf",
            bbox_inches="tight")

#%% Plot in the contour of the irreducible Brillouin zone
plt.figure(figsize=(4, 2.2))
ky_max = pi/np.sqrt(3)
ky1 = np.linspace(0, ky_max, 50)
ky2 = ky_max*np.ones((50))
ky3 = np.linspace(ky_max, 0, 50)
ky_vec = np.zeros((150, 1))
ky_vec[0:50,0] = ky1
ky_vec[50:100,0] = ky2
ky_vec[100:150,0] = ky3
kx1 = np.zeros((50))
kx2 = np.linspace(0, pi/2)
kx3 = np.linspace(pi/2, 0, 50)
kx_vec = np.zeros((150, 1))
kx_vec[0:50,0] = kx1
kx_vec[50:100,0] = kx2
kx_vec[100:150,0] = kx3
omega = disp_rel(kx_vec, ky_vec)
xaxis = np.linspace(0, 150, 150)
omega_low = np.sqrt(omega[0,:,:] + omega[1,:,:])
omega_high = np.sqrt(omega[2,:,:] + omega[3,:,:])
# Vertical lines
ymax = 5
plt.plot([0,0], [0,ymax], 'gray')
plt.plot([50,50], [0,ymax], 'gray')
plt.plot([100,100], [0,ymax], 'gray')
# Parts of the plot
plt.plot(xaxis, omega_low, 'k')
plt.plot(xaxis, omega_high, 'k')
plt.xticks([0, 50, 100, 150], [r"$\Gamma$", r"$M$", r"$K$", r"$\Gamma$"])
plt.xlim(0, 150)
plt.ylabel(r"$\Omega$")
plt.xlabel(r"$k$-space position")
plt.savefig("Notes/img/hexagon-irreducible.pdf", bbox_inches="tight")
plt.show()


omega2 = freqs(kx_vec, ky_vec)
omega2 = np.abs(omega2)
from mayavi import mlab
mlab.mesh(kx_vec, ky_vec, np.sqrt(omega2[0,:,:] + omega2[1,:,:]), colormap='hot')
mlab.mesh(kx_vec, ky_vec, np.sqrt(omega2[2,:,:] + omega2[3,:,:]), colormap='hot')
#mlab.axes()
mlab.show()

