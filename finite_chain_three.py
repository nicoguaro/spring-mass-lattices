"""
Compute the spectrum response for a chain of $N$
unit cells composed of three different masses
three equal springs.


@author: Nicolas Guarin-Zapata
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse.linalg import dsolve
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 8


def n_springs(N, mu_1, mu_2, f_vec):
    """
    Compute the ratio $u_N/f_0$ for a chain
    of $N$ springs. The range of frequencies
    is given by $[f_{\min}, f_{\max}]$. $\mu_1$ is
    the ratio of masses $\mu_1 = m_2/m_1$, $\mu_2$ is
    the ratio of masses $\mu_2 = m_3/m_1$.
    """
    nfreq = len(f_vec)
    uf_ratio = np.zeros(nfreq)
    force = np.zeros(3*N + 1)
    force[0] = -1.0
    for k, freq in enumerate(f_vec):
        d_main = np.zeros(3*N + 1)
        d_main[1:-1:3] = -2. + mu_1*freq**2
        d_main[2:-1:3] = -2. + mu_2*freq**2
        d_main[3:-1:3] = -2. + freq**2
        d_main[0] = -1.0 + freq**2
        d_main[-1] = -1.0  + freq**2
        d_low = np.ones(3*N + 1)
        data = [d_low, d_main, d_low]
        diags = [-1,0,1]
        A = sparse.spdiags(data, diags, 3*N+1, 3*N+1, format='csc')


        x = dsolve.spsolve(A, force)
        uf_ratio[k] = x[-1]


    return uf_ratio

N_vec = [1, 5, 10, 50]
npts = 10000
mu_1 = 2.
mu_2 = 3.
f_vec = np.linspace(1e-3, 2.5, npts)
plt.figure(figsize=(3, 1.9))
for N in N_vec:
    uf_ratio = n_springs(N, mu_1, mu_2, f_vec)
    plt.semilogy(f_vec, abs(uf_ratio), label="N=%i"%N)

plt.ylim(1e-9, 1e3)
plt.yticks(np.logspace(-9, 3, 5))
plt.grid(True, color="gray", alpha=0.3)
plt.xlabel(r"$\Omega$")
plt.ylabel(r"$u_N/\hat{f}_0$")
plt.legend(loc="best", fontsize=6)
plt.savefig("Notes/img/three_finite-m1=%g-m2=%g.pdf"%(mu_1, mu_2),
            bbox_inches="tight")

plt.show()