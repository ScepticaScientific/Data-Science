#!/home/ubuntu/miniconda/bin/python
#
# This code implements the computation of the auto- or cross- power spectrum using Welch's method.
# The code duplicates the corresponding code of MATLAB (see 'pwelch()' and 'cpsd()' there).
#
# At the input, 'x' and 'y' are time series of the physical observables 'x(t)' and 'y(t)'. If 'y' is omitted
# then the auto-spectrum of 'x' is computed; if 'y' is provided then the cross-spectrum is computed.
#
# At the output, 'Pxx' is the power spectrum, while 'freq' is the frequency range over which the spectrum is computed.
#
# The end user is granted perpetual permission to reproduce, adapt, and/or distribute this code, provided that
# an appropriate link is given to the original repository it was downloaded from.

import numpy as np

def getPowerSpectrum(x, y = []):
    ## Initialisation
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
    Sxx = np.zeros((nfft))
    if (len(y) > 0):
        for i in range(1, nSegments + 1):
            [Sxxk, freq] = getPeriodogram(x[iFrom[i - 1] - 1 : iTo[i - 1]], y[iFrom[i - 1] - 1 : iTo[i - 1]], wnd, nfft, 1.0)
            Sxx = Sxx + Sxxk
    else:
        for i in range(1, nSegments + 1):
            [Sxxk, freq] = getPeriodogram(x[iFrom[i - 1] - 1 : iTo[i - 1]], [], wnd, nfft, 1.0)
            Sxx = Sxx + Sxxk
    # ... and then determine the averaged periodogram all over the segments
    Sxx = Sxx / nSegments

    # Finally, we adjust the averaged periodogram to obtain the spectrum
    [Pxx, freq] = getSpectrum(Sxx, freq, nfft, 1.0)

    return [Pxx, freq]

## Auxiliaries
# Adjustment of spectrum
def getSpectrum(Sxx, freq, nfft, fs):
    # We keep only the values corresponding to the frequency  range[0, pi] or [0, pi)
    if (np.mod(nfft, 2) == 1):
        inds = np.arange(1, np.int((nfft + 1) / 2 + 1))

        Sxx = 2.0 * Sxx[inds - 1]
        Sxx[1 - 1] = Sxx[1 - 1] / 2.0
    else:
        inds = np.arange(1, np.int(nfft / 2 + 1 + 1))

        Sxx = 2.0 * Sxx[inds - 1]
        Sxx[1 - 1] = Sxx[1 - 1] / 2.0
        Sxx[-1] = Sxx[-1] / 2.0

    freq = freq[inds - 1]

    # Normalisation to the sampling frequency
    Pxx = Sxx / fs

    return [Pxx, freq]

# Computation of periodogram
def getPeriodogram(x, y, wnd, nfft, fs):
    # We compute the Fourier transform(s) of the signal(s)
    Fx = np.fft.fft(np.multiply(x, wnd), nfft)

    if (len(y) != 0):
        Fy = np.fft.fft(np.multiply(y, wnd), nfft)

    # We compute the window's energy
    U = np.sum(wnd ** 2)

    # Finally, we compute the auto- or the cross- periodogram ...
    if (len(y) != 0):
        Pxx = np.multiply(Fx, np.conj(Fy)) / U
    else:
        Pxx = np.multiply(Fx, np.conj(Fx)) / U
    # ... and obtain the frequency range
    freq = getFrequencies(nfft, fs)

    return [Pxx, freq]

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

# The closest power of two such that two to the power is greater than, or equal to, the given value
def nextpow2(val):
    p = 0
    while 2 ** p < val:
        p = p + 1

    return p