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
rcParams['font.size'] = 8

#%% Dispersion curve
ka = np.linspace(-2*pi, 2*pi, 101)
omega = 2*np.abs(np.sin(0.5*ka))

plt.figure(figsize=(4, 2.2))
plt.plot([-1, -1],[0, 2.5], c='gray')
plt.plot([1, 1],[0, 2.5], 'gray')
plt.plot([-1/pi, -1/pi],[0, 2.5], '--', c='gray')
plt.plot([1/pi, 1/pi],[0, 2.5], '--', c='gray')
plt.plot(ka/pi, omega, 'k')
plt.xlabel(r"$ka/\pi$")
plt.ylabel(r"$\Omega$")
plt.ylim(0, 2.5)
plt.savefig("../img/spring-mass-plot.pdf", bbox_inches="tight")


#%% Speed curves
ka = np.linspace(1e-6, pi, 51)
vp = 2/ka*np.abs(np.sin(0.5*ka))
vg = np.cos(0.5*ka)*np.sign(np.sin(0.5*ka))

plt.figure(figsize=(4, 2.2))
plt.plot([1/pi, 1/pi],[0, 2.5], '--', c='gray')
plt.plot(ka/pi, vp, 'k', label="Phase speed")
plt.plot(ka/pi, vg, 'b', label="Group speed")
plt.xlabel(r"$ka/\pi$")
plt.ylabel(r"Speed - $\Omega/(ka)$")
plt.ylim(0, 1.5)
plt.legend(fontsize=8)
plt.savefig("../img/spring-mass-speeds.pdf", bbox_inches="tight")

plt.show()
