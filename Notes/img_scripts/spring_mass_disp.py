# -*- coding: utf-8 -*-
"""
Generate the plot for the dispersion curve

@author: Nicolas Guarin-Zapata
"""
from __future__ import division
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

#%% Dispersion curve
ka = np.linspace(-2*pi, 2*pi, 101)
omega = 2*np.abs(np.sin(0.5*ka))

plt.plot([-1, -1],[0, 2.5], c='gray', lw=2)
plt.plot([1, 1],[0, 2.5], 'gray', lw=2)
plt.plot([-1/pi, -1/pi],[0, 2.5], '--', c='gray', lw=2)
plt.plot([1/pi, 1/pi],[0, 2.5], '--', c='gray', lw=2)
plt.plot(ka/pi, omega, 'k', lw=2)
plt.xlabel(r"$ka/\pi$", size=18)
plt.ylabel(r"$\Omega$", size=18)
plt.ylim(0, 2.5)
plt.savefig("../img/spring-mass-plot.pdf")


#%% Speed curves
ka = np.linspace(1e-6, pi, 51)
vp = 2/ka*np.abs(np.sin(0.5*ka))
vg = np.cos(0.5*ka)*np.sign(np.sin(0.5*ka))

plt.figure()
plt.plot([1/pi, 1/pi],[0, 2.5], '--', c='gray', lw=2)
plt.plot(ka/pi, vp, 'k', lw=2, label="Phase speed")
plt.plot(ka/pi, vg, 'b', lw=2, label="Group speed")
plt.xlabel(r"$ka/\pi$", size=18)
plt.ylabel(r"Speed - $\Omega/(ka)$", size=18)
plt.ylim(0, 1.5)
plt.legend()
plt.savefig("../img/spring-mass-speeds.pdf")

plt.show()
