
"""
Compute the spectrum response for a chain of $N$
mass and springs


@author: Nicolas Guarin-Zapata
@date: January 23, 2015
"""

import numpy as np

def n_springs(N, f_vec):
    """
    Compute the ratio $u_N/f_0$ for a chain
    of $N$ springs. The range of frequencies
    is given by $[f_{\min}, f_{\max}]$.
    """
    npts = len(f_vec)
    A_mat = np.zeros((N, N))
    
    #  Fill the matrices
    
    # Make a loop to compute the solution for every frequency
    
    return 2*f_vec # Return solution
    
    
f_vec = np.linspace(1e-6, 5)
uf_ratio = n_springs(10, f_vec)