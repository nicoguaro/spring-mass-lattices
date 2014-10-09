#!/usr/bin/python
"""
Dispersion for spring mass system
"""
import numpy as np
import matplotlib.pyplot as plt

a = 1
c = 10
m1 = 1
m2 = 2
pi = np.pi

k = np.linspace(-2*pi, 2*pi, 200)

omega2_p = c*( 1. + np.exp(1.0j*k*a) )/(np.exp(1./2.*1.0j*k*a))/np.sqrt(m1*m2)
omega2_m = -omega2_p
omega_p = np.abs(np.sqrt(omega2_p))
omega_m = np.abs(np.sqrt(omega2_m))



plt.plot(k/pi, omega_p, k/pi, omega_m)
plt.show()


