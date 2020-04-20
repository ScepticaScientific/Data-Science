#!/home/ubuntu/miniconda/bin/python
#
# This code implements the canonical coherence analysis (CCA) of multivariate data. The computation is performed using
# an accurate non-parametric Fourier spectrum estimation different from the standard Python function scipy.signal.csd().
# For details on CCA, please refer to [1].
#
# INPUT:
# xx      - a multivariate time series of second increments of the physical observables 'x(t) = (x_1(t), ..., x_N(t))'
#
# OUTPUT:
# evt     - the total coherence coefficient over the frequency range
# ev      - the array of partial coherence coefficients over the frequency range
# freq    - the frequency range over which the coherence spectrum is computed
#
# REFERENCES:
# [1] A.A. Lyubushin, Data Analysis of Systems of Geophysical and Ecological Monitoring, Nauka, Moscow, 2007.
# [2] D.M. Filatov and A.A. Lyubushin, Physica A: Stat. Mech. Appl., 527 (2019) 121309. DOI: 10.1016/j.physa.2019.121309.
#
# The end user is granted perpetual permission to reproduce, adapt, and/or distribute this code, provided that
# an appropriate link is given to the original repository it was downloaded from.

import scipy.linalg as sl
import numpy as np
from numpy import linalg

import sys
sys.path.append('../utils')
from getPowerSpectrum import getPowerSpectrum

def getCanonicalCoherence(xx):
    ## Initialisation
    N = xx.shape[2 - 1]     # Number of variates in the vector time series

    # Auxiliary computation (of the frequency range)
    [_, freq] = getPowerSpectrum(xx[:, 1 - 1], xx[:, 1 - 1])
    freq_len = np.int(len(freq))

    ## Computing
    ev = np.zeros((freq_len, N))
    for ic in range(1, N + 1):
        # We split the original N-variate signal into an (N - 1)-variate and a single-variate ones
        x = xx[:, list(range(1 - 1, ic - 1)) + list(range(ic + 1 - 1, N))]      # (N - 1)-variate (i.e. vector) time series
        y = xx[:, range(ic - 1, ic)]                                            # Single-variate (i.e. scalar) time series

        # We compute the spectral matrices (see [1, p. 111]) ...
        # ... for the first signal, ...
        Sxx = np.zeros((N - 1, N - 1, freq_len), dtype = np.complex128)
        for i in range(1, N):
            for j in range(1, N):
                [Sxx[i - 1, j - 1, :], _] = getPowerSpectrum(x[:, i - 1], x[:, j - 1])

        # ... for the mixture of the first and second signals, ...
        Sxy = np.zeros((N - 1, 1, freq_len), dtype = np.complex128)
        for i in range(1, N):
            [Sxy[i - 1, 1 - 1, :], _] = getPowerSpectrum(x[:, i - 1], y[:, 1 - 1])

        Syx = np.zeros((1, N - 1, freq_len), dtype = np.complex128)
        for j in range(1, N):
            [Syx[1 - 1, j - 1, :], _] = getPowerSpectrum(y[:, 1 - 1], x[:, j - 1])

        # ... and for the second signal
        Syy = np.zeros((1, 1, freq_len), dtype = np.complex128)
        [Syy[1 - 1, 1 - 1, :], _] = getPowerSpectrum(y[:, 1 - 1], y[:, 1 - 1])

        # We compute the matrix whose eigenvalues are the measures of coherence and explicitly determine
        # the maximum value of the coherence at each frequency
        for i in range(1, freq_len + 1):
            # We implement formula (2.3.1) from [1]
            m1 = linalg.inv(sl.sqrtm(Sxx[:, :, i - 1]))
            m2 = Sxy[:, :, i - 1]
            m3 = linalg.inv(Syy[:, :, i - 1])
            m4 = Syx[:, :, i - 1]
            m5 = linalg.inv(sl.sqrtm(Sxx[:, :, i - 1]))
            U = np.matmul(np.matmul(np.matmul(np.matmul(m1, m2), m3), m4), m5)

            ev[i - 1, ic - 1] = max(np.real(linalg.eig(U)[1 - 1]))

        # Normalisation
        ev[:, ic - 1] = ev[:, ic - 1] / max(np.abs(ev[:, ic - 1]))

    # Finally we compute the total coherence frequency-wise
    evt = np.sqrt(np.prod(ev, axis = 2 - 1))

    return [evt, ev, freq]
