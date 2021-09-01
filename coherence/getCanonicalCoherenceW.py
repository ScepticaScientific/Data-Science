#!/home/ubuntu/miniconda/bin/python
#
# This code implements canonical coherence analysis (CCA) of multivariate data. The computation is performed using
# a Morlet wavelet spectrum estimation. For details on CCA, please refer to [1-2].
#
# At the input:
#   - 'ddx' is a multivariate stationary time series (usually the second derivative of a physical observable 'x(t)')
#     'ddx(t) = (ddx_1(t), ..., ddx_N(t))'; the variates are provided column-wise
#   - 'fs' is the sampling rate (optional)
#   - 'isInfo' is the flag prescribing to output a message for each variate processed (optional).
#
# Besides, for an array of time moments this function allows to determine the left and right borders, out of which
# the daughter wavelet, while moving in time at a fixed scale, does not involve those time moments' data samples
# for the computation of CWT:
#   - 'timesOfInterest' is a row vector of the time moments (these are usually extreme events reflected in the data)
#     for which the left and right borders are wanted (optional)
#   - 'energyThreshold' is the probability with which the daughter wavelet does not affect the time moments' data
#     samples (optional).
#
# At the output:
#   - 'evt' is the total coherence spectrum (the coherence coefficient between 0.0 and 1.0)
#   - 'ev' is the array of the partial coherence coefficients
#   - 'freqS' is the scale-related frequency range over which the coherence has been computed
#   - 'coi' is the cone of influence
#   - 'timeBorders' are borders of the energy cones out of which the above-specified time moments' data samples are
#     not affected by the daughter wavelet (optional).
#
# REFERENCES:
# [1] A.A. Lyubushin, Data Analysis of Systems of Geophysical and Ecological Monitoring, Nauka, Moscow, 2007.
# [2] A.A. Lyubushin, in: Complexity of Seismic Time Series: Measurement and Applications, Elsevier, Amsterdam, 2018.
#                     DOI: 10.1016/B978-0-12-813138-1.00006-7.
#
# The end user is granted perpetual permission to reproduce, adapt, and/or distribute this code, provided that
# an appropriate link is given to the original repository it was downloaded from.

import numpy as np
from numpy import linalg
import scipy.special as sspec

import sys
sys.path.append('../utils')
from getPowerSpectrumW import getPowerSpectrumW

def getCanonicalCoherenceW(ddx, fs = 1.0, timesOfInterest = [], energyThreshold = [], isInfo = False, isSmoothing = True):
    ## Initialisation
    N = ddx.shape[2 - 1]     # Number of variates in the vector time series

    # Auxiliary computation (of the scale-related frequency range and the cone of influence)
    [_, freqS, coi] = getPowerSpectrumW(ddx[:, [1 - 1, 1 - 1]], fs, isSmoothing)
    freqS_len = len(freqS)
    freqT_len = ddx.shape[1 - 1]

    ## Computing
    ev = np.zeros((freqS_len, freqT_len, N))
    for ic in range(1, N + 1):
        # We split the original N-variate signal into an (N - 1)-variate and a single-variate ones
        x = ddx[:, list(range(1 - 1, ic - 1)) + list(range(ic + 1 - 1, N))]      # (N - 1)-variate (i.e. vector) time series
        y = ddx[:, range(ic - 1, ic)]                                            # Single-variate (i.e. scalar) time series

        # We compute the spectral matrices (see [1, p. 111]) ...
        # ... for the first signal, ...
        Sxx = np.zeros((N - 1, N - 1, freqS_len, freqT_len), dtype = np.complex128)
        for i in range(1, N):
            for j in range(1, N):
                [Sxx[i - 1, j - 1, :, :], _, _] = getPowerSpectrumW(x[:, [i - 1, j - 1]], fs, isSmoothing)

        # ... for the mixtures of the first and second signals, ...
        Sxy = np.zeros((N - 1, 1, freqS_len, freqT_len), dtype = np.complex128)
        for i in range(1, N):
            [Sxy[i - 1, 1 - 1, :, :], _, _] = getPowerSpectrumW(np.asarray([x[:, i - 1], y[:, 1 - 1]]).transpose(), fs, isSmoothing)

        Syx = np.zeros((1, N - 1, freqS_len, freqT_len), dtype = np.complex128)
        for j in range(1, N):
            [Syx[1 - 1, j - 1, :, :], _, _] = getPowerSpectrumW(np.asarray([y[:, 1 - 1], x[:, j - 1]]).transpose(), fs, isSmoothing)

        # ... and for the second signal
        Syy = np.zeros((1, 1, freqS_len, freqT_len), dtype = np.complex128)
        [Syy[1 - 1, 1 - 1, :, :], _, _] = getPowerSpectrumW(y[:, [1 - 1, 1 - 1]], fs, isSmoothing)

        # Now we are to compute the matrix whose eigenvalues are measures of coherence and to explicitly determine
        # the maximum coherence at each frequency. (The implementation is vectorised to speed up the calculations.)
        Sxx_aux = toMatrixArray(Sxx, N - 1, N - 1)
        Sxx_aux = np.asarray([linalg.inv(Sxxi) for Sxxi in Sxx_aux])

        Sxy_aux = toMatrixArray(Sxy, N - 1, 1)

        Syy_aux = toMatrixArray(Syy, 1, 1)
        Syy_aux = np.asarray([linalg.inv(Syyi) for Syyi in Syy_aux])

        Syx_aux = toMatrixArray(Syx, 1, N - 1)

        aux = np.asarray([np.matmul(Sxxi, Sxyi) for Sxxi, Sxyi in zip(Sxx_aux, Sxy_aux)])
        aux = np.asarray([np.matmul(auxi, Syyi) for auxi, Syyi in zip(aux, Syy_aux)])
        U = np.asarray([np.matmul(auxi, Syxi) for auxi, Syxi in zip(aux, Syx_aux)])

        aux = np.asarray([linalg.eig(Ui)[1 - 1] for Ui in U])
        aux = np.asarray([np.real(auxi) for auxi in aux])
        ev[:, :, ic - 1] = np.reshape(np.asarray([max(auxi) for auxi in aux]), (freqS_len, freqT_len), order = 'F')

        # Informational message to the standard output
        if (isInfo):
            print('CCWA: variate %d of %d processed' % (ic, N))

    # Finally we compute the total coherence frequency-wise
    evt = np.prod(ev, axis = 3 - 1) ** (1.0 / N)

    # We obtain the energy cones out of which the moving daughter wavelet does not affect the specific time moments (if any)
    if (energyThreshold):
        timeBorders = getBorders(freqS, energyThreshold, timesOfInterest, fs)
        return [evt, ev, freqS, coi, timeBorders]

    return [evt, ev, freqS, coi]

## Auxiliaries
# Reshaping a 4-D array into a sequence of 2-D matrices
def toMatrixArray(mtrx, nrows, ncols):
    aux = mtrx.reshape((mtrx.shape[1 - 1], mtrx.shape[4 - 1] * mtrx.shape[3 - 1] * mtrx.shape[2 - 1]), order = 'F')

    return aux.reshape((1, nrows, -1, ncols)).swapaxes(1, 2).reshape((-1, nrows, ncols))

# Computing the time borders (or energy cones) for specific time moments
def getBorders(freq, energyThreshold, timesOfInterest, fs):
    sigma = 6.0     # The Morlet wavelet parameter 'sigma' (see function MorletCWT() in file '../utils/getPowerSpectrumW.py')
    FourierFactor = 2.0 * np.pi / sigma

    scales = 1 * fs / freq / FourierFactor      # In samples
    t = sspec.erfinv(energyThreshold)           # In relative units
    t_deltas = t * scales / fs                  # In seconds

    if (np.isscalar(timesOfInterest)):
        timesOfInterest = np.array([timesOfInterest])
    N = len(timesOfInterest)

    timeBorders = np.zeros((len(freq), 2, N))
    for i in range(1, N + 1):
        timeBorders[:, 1 - 1, i - 1] = timesOfInterest[i - 1] - t_deltas
        timeBorders[:, 2 - 1, i - 1] = timesOfInterest[i - 1] + t_deltas

    return timeBorders