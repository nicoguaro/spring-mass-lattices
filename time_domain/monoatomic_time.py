# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 13:22:05 2017

@author: nicoguaro
"""
from __future__ import division, print_function
import numpy as np
from numpy import (linspace, exp, sqrt, pi, arcsin, outer)
from numpy.fft.fftpack import ifft
from numpy.fft.helper import fftshift
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import rcParams


rcParams['animation.convert_path'] = "C:\Program Files\ImageMagick-6.9.3\convert.exe"

nx = 5001
nfreq = 2048
x = linspace(-20, 20, nx)
omega = linspace(0, 4, nfreq)
kappa = 2*arcsin(abs(omega/2) + 0j)
ricker_amp = sqrt(pi)*4*omega**2*exp(-omega**2)
TX = exp(1j*outer(ricker_amp*kappa, x))
X = ifft(TX, axis=0)
magX = np.real(X)


def update(k):
    ax.cla()
    plt.plot(x, magX[k, :])
    plt.ylim(-1, 2)


nframes = int(nfreq/2)
fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot(111)
ani = animation.FuncAnimation(fig, update, frames=range(nframes), blit=False,
                              repeat=False)

ani.save("monoatomic_time.gif", writer='imagemagick',
         savefig_kwargs={'delay': 6})
#plt.show()