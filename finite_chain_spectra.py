
"""
Compute the spectrum response for a chain of $N$
mass and springs


@author: Nicolas Guarin-Zapata
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse.linalg import dsolve
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 8


def n_springs(N, f_vec):
    """
    Compute the ratio $u_N/f_0$ for a chain
    of $N$ springs. The range of frequencies
    is given by $[f_{\min}, f_{\max}]$.
    """
    nfreq = len(f_vec)
    uf_ratio = np.zeros(nfreq)
    for k, freq in enumerate(f_vec):
        d_main = -2.0*np.ones(N) + freq**2
        d_main[0] = -1.0 + freq**2
        d_main[-1] = -1.0  + freq**2
        d_low = np.ones(N)
        data = [d_low, d_main, d_low]
        diags = [-1,0,1]
        A = sparse.spdiags(data, diags, N, N, format='csc')
        force = np.zeros(N)
        force[0] = -1.0

        x = dsolve.spsolve(A, force)
        uf_ratio[k] = x[-1]


    return uf_ratio

plt.figure(figsize=(4, 2.2))
N_vec = [1, 5, 10,  50]
npts = 10000
f_vec = np.linspace(1e-3, 5, npts)
for N in N_vec:
    uf_ratio = n_springs(N, f_vec)
    plt.semilogy(f_vec, abs(uf_ratio), label="N=%i"%N)
    plt.ylim(1e-12, 1e3)
    plt.yticks(np.logspace(-12, 3, 6))

plt.grid(True, color="gray", alpha=0.3)
plt.xlabel(r"$\Omega$")
plt.ylabel(r"$u_N/\hat{f}_0$")
plt.legend(loc="best", fontsize=8)
plt.savefig("Notes/img/single_finite.svg")
plt.savefig("Notes/img/single_finite.pdf", bbox_inches="tight")

plt.show()