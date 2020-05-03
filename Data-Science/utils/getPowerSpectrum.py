#!/home/ubuntu/miniconda/bin/python
#
# This code implements the computation of the continuous Fourier auto- or cross- power spectrum using Welch's method.
# The code is based on MATLAB's functions pwelch() and cpsd().
#
# At the input, 'x' is a uni- or bivariate time series of the physical observables 'x_1(t)' and, possibly, 'x_2(t)',
# while 'fs' is the sampling rate. If 'x_2(t)' is omitted then the auto-spectrum of 'x_1(t)' is computed; otherwise
# the cross-spectrum between 'x_1(t)' and 'x_2(t)' is computed. The time series are to be provided column-wise
# for each variate.
#
# At the output, 'Ps' is the power spectrum, while 'freq' is the frequency range over which the spectrum is computed.
#
# The end user is granted perpetual permission to reproduce, adapt, and/or distribute this code, provided that
# an appropriate link is given to the original repository it was downloaded from.

import numpy as np

def getPowerSpectrum(x, fs = 1.0):
    # Initialisation
    x_len = x.shape[1 - 1]

    wndTotal = np.int(np.floor(x_len / 4.5))            # Window length
    wndOverlap = np.int(np.floor(0.5 * wndTotal))       # Length of the window overlap
    wnd = np.hamming(wndTotal)                          # Window

    nfft = max(256, 2 ** nextpow2(wndTotal))

    # We determine the number of segments into which we shall then split
    # the entire time interval of the data
    nSegments = np.int(np.floor((x_len - wndOverlap) / (wndTotal - wndOverlap)))

    wndNonOverlap = wndTotal - wndOverlap
    iFrom = np.arange(1, nSegments * wndNonOverlap + 1, wndNonOverlap)
    iTo = iFrom + wndTotal - 1

    # We compute periodograms for each segment of the time interval in a moving window ...
    Per = np.zeros(nfft)

    if (len(x.shape) == 1):
        for i in range(1, nSegments + 1):
            [Per_i, freq] = getPeriodogram(x[iFrom[i - 1] - 1 : iTo[i - 1]], [], wnd, nfft, fs)
            Per = Per + Per_i
    else:
        if (x.shape[2 - 1] == 1):
            for i in range(1, nSegments + 1):
                [Per_i, freq] = getPeriodogram(x[iFrom[i - 1] - 1 : iTo[i - 1], 1 - 1], [], wnd, nfft, fs)
                Per = Per + Per_i
        else:
            for i in range(1, nSegments + 1):
                [Per_i, freq] = getPeriodogram(x[iFrom[i - 1] - 1 : iTo[i - 1], 1 - 1], x[iFrom[i - 1] - 1 : iTo[i - 1], 2 - 1], wnd, nfft, fs)
                Per = Per + Per_i
    # ... and then determine the averaged periodogram all over the segments
    Per = Per / nSegments

    # Finally, we adjust the averaged periodogram to obtain the spectrum
    [Ps, freq] = getSpectrum(Per, freq, nfft)

    return [Ps, freq]

## Auxiliaries
# Computation of periodogram
def getPeriodogram(x, y, wnd, nfft, fs):
    # We compute the continuous window's energy
    U = np.sum(wnd ** 2.0) / fs

    # We compute the continuous Fourier transform(s) of the signal(s)
    Fx = np.fft.fft(np.multiply(x, wnd), nfft) / fs

    if (len(y) != 0):
        Fy = np.fft.fft(np.multiply(y, wnd), nfft) / fs
    else:
        Fy = Fx

    # Finally, we compute the auto- or the cross-periodogram ...
    Per = np.multiply(Fx, np.conj(Fy)) / U
    # ... and obtain the frequency range
    freq = getFrequencies(nfft, fs)

    return [Per, freq]

# Adjustment of frequency range
def getFrequencies(N, fs):
    freqResolution = fs / N
    freq = freqResolution * np.arange(0, N - 1 + 1)

    freqNyquist = fs / 2.0
    freqResolution2 = freqResolution / 2.0

    if (np.mod(N, 2) == 1):
        N2 = np.int((N + 1) / 2)

        freq[N2 - 1] = freqNyquist - freqResolution2
        freq[N2 + 1 - 1] = freqNyquist + freqResolution2
    else:
        N2 = np.int((N / 2) + 1)

        freq[N2 - 1] = freqNyquist

    freq[N - 1] = fs - freqResolution

    return freq

# Adjustment of periodogram to get spectrum
def getSpectrum(Per, freq, nfft):
    # We keep only the values corresponding to the frequency range [0, pi] or [0, pi)
    if (np.mod(nfft, 2) == 1):
        inds = np.arange(1, np.int((nfft + 1) / 2 + 1))

        Per = 2.0 * Per[inds - 1]
        Per[1 - 1] = Per[1 - 1] / 2.0
    else:
        inds = np.arange(1, np.int(nfft / 2 + 1 + 1))

        Per = 2.0 * Per[inds - 1]
        Per[1 - 1] = Per[1 - 1] / 2.0
        Per[-1] = Per[-1] / 2.0

    freq = freq[inds - 1]

    return [Per, freq]

# The closest power of two such that two to the power is greater than, or equal to, the given value
def nextpow2(val):
    p = 0
    while 2 ** p < val:
        p = p + 1

    return p