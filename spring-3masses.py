# -*- coding: utf-8 -*-
"""
Compute the dispersion curves for a spring mass
lattice formed with three different masses.

The parameters are:
    - mu_1: ratio of masses m2/m1
    - mu_2: ratio of masses m3/m1
    - ka: dimensionless wavenumber, a is the original separation
      between springs.
    - omega: is the dimensionless frequency, divided by omega_0
      omega_0**2 = c/m1.

@author: Nicolas Guarin-Zapata
@date: October 9, 2014
"""
import numpy as np
import scipy.linalg as LA
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

def disp_rel(kvec, mu_1, mu_2):
    n = kvec.shape[0]
    omega = np.zeros((3, n))
    for j,ka in enumerate(kvec):
        K = np.array([[-2., 1., np.exp(-1j*ka)],
                      [1., -2., 1.],
                      [np.exp(1j*ka), 1, -2]])
        M = -np.array([[1., 0., 0.],
                      [0., mu_1, 0.],
                      [0., 0., mu_2]])
        vals = LA.eigvalsh(-K, -M)
        omega[:,j] = vals
        
    return omega

mu_1 = 2.
mu_2 = 10.
kvec = np.linspace(0, np.pi, 501)
omega = np.sqrt(disp_rel(kvec, mu_1, mu_2))

plt.figure()
plt.plot(kvec/np.pi, omega[0,:])
plt.plot(kvec/np.pi, omega[1,:])
plt.plot(kvec/np.pi, omega[2,:])
plt.xlabel(r"$ka/\pi$", size=18)
plt.ylabel(r"$\omega/\omega_0$", size=18)
plt.ylim(0,2.5)
plt.show()