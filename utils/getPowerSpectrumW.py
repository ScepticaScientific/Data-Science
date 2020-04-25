#!/home/ubuntu/miniconda/bin/python
#
# This code implements the computation of the continuous wavelet power auto- or cross-spectrum. It is based on
# MATLAB's function cwtft() (see also wcoherence() there).
#
# At the input, 'x' is a uni- or bivariate time series of the physical observables 'x_1(t)' and, possibly, 'x_2(t)',
# while 'fs' is the sampling rate. If 'x_2(t)' is omitted then the auto-spectrum of 'x_1(t)' is computed; otherwise
# the cross-spectrum between 'x_1(t)' and 'x_2(t)' is computed. The time series are to be provided column-wise
# for each variate.
#
# At the output, 'Ps' is the power spectrum; 'freqT' is the time-related frequency range; 'freqS' is the scale-related
# frequency range.
#
# The end user is granted perpetual permission to reproduce, adapt, and/or distribute this code, provided that
# an appropriate link is given to the original repository it was downloaded from.

import numpy as np
import scipy.signal as ss

def getPowerSpectrumW(x, fs = 1.0):
    # We compute the continuous wavelet transform of the first time series ...
    if (len(x.shape) == 1):
        [cwtX, freqT, freqS, scales] = getCWT(x, fs)
    else:
        [cwtX, freqT, freqS, scales] = getCWT(x[:, 1 - 1], fs)

    if (len(x.shape) == 1):
        cwtY = cwtX
    else:
        if (x.shape[2 - 1] == 2):
            [cwtY, _, _, scales] = getCWT(x[:, 2 - 1], fs)
        else:
            cwtY = cwtX
    # ... and then multiply it with the corresponding conjugate
    Ps = np.multiply(cwtX, np.conj(cwtY))         # Inner product of 'cwtX' with its conjugate

    # Optional smoothing and filtering. This line can be remarked.
    #Ps = smoothSpectrum(Ps, scales, fs)

    return [Ps, freqT, freqS]

## Auxiliaries
# Computation of continuous wavelet transform
def getCWT(x, fs):
    # Initialisation
    x_len = len(x)

    fpo = 12    # Frequencies per octave
    a0 = 2.0 ** (1.0 / fpo)
    maxNumOctaves = np.int(np.floor(np.log2(x_len)) - 1)
    scales = 2.0 / fs * a0 ** np.arange(0, maxNumOctaves * fpo + 1)

    # We expand the original time series to a power of two and then symmetrise it
    extP2 = np.int(1 + np.floor(np.log2(x_len) + 0.4999))
    xExpanded = symmetrise(x - np.mean(x), 2 ** extP2 - x_len)

    # We determine the frequency range for the expanded data
    xExtended_len = len(xExpanded)
    freqT = 2.0 * np.pi * fs / xExtended_len * np.arange(1, np.int(np.floor(xExtended_len / 2) + 1))
    aux_inds = np.arange(np.int(np.floor((xExtended_len - 1) / 2)), 1 - 1, -1)
    freqT = np.concatenate([np.array([0.0]), freqT, -freqT[aux_inds - 1]])

    # We compute the continuous Fourier transform of the expanded data, ...
    xExFFT = np.fft.fft(xExpanded) / fs
    # ... prepare the Fourier transform of the Morlet wavelet, ...
    [wft, freqS] = MorletFourierTransform(freqT, scales)
    # ... compute the continuous wavelet transform of the expanded data in the frequency domain
    # and then restore the wavelet-transformed time series in the time domain
    cwt = np.fft.ifft(np.multiply(np.tile(xExFFT, (len(scales), 1)), wft), axis = 2 - 1)

    # We normalise the continuous wavelet transform scale-wise
    scaleFactors = np.tile(1.0 / np.sqrt(scales), (cwt.shape[2 - 1], 1)).transpose()
    cwt = np.multiply(scaleFactors, cwt)

    # Finally we keep only those values that correspond to the original data's frequencies
    cwt = cwt[:, 1 - 1 : x_len]
    freqT = freqT[1 - 1 : x_len]

    return [cwt, freqT, freqS, scales]

# Expansion and symmetrisation of time series
def symmetrise(x, x_len):
    xl = len(x)

    if (xl == 1):
        return np.tile(x, 1, x_len + 1)

    while (x_len + 1 > xl):
        x_len = x_len - (xl - 1)
        x = symmetrise(x, xl - 1)
        xl = len(x)

    inds = np.concatenate([np.arange(1, xl + 1), np.arange(xl - 1, xl - x_len - 1, -1)])
    return x[inds - 1]

# Fourier transform of the Morlet wavelet
def MorletFourierTransform(freqT, scales):
    nFreq = len(freqT)
    nScales = len(scales)

    wft = np.zeros((nScales, nFreq))

    amplitude = (1.0 / np.pi ** 0.25) * np.sqrt(freqT[2 - 1]) * np.sqrt(nFreq)
    for j in range(1, nScales + 1):
        expnt = np.multiply(-(scales[j - 1] * freqT - 6.0) ** 2.0 / 2.0, freqT > 0.0)
        wft[j - 1, :] = amplitude * np.sqrt(scales[j - 1]) * np.multiply(np.exp(expnt), freqT > 0.0)

    freqS = 3.0 / np.pi / scales

    return [wft, freqS]

# Smoothing the spectrum
def smoothSpectrum(mtrx, scales, fs):
    N = mtrx.shape[2 - 1]

    nPads = 2 ** nextpow2(N)    # We have to enlarge the number of samples to the closest power of two

    freqT = 2.0 * np.pi * fs / nPads * np.arange(1, np.int(np.floor(nPads / 2) + 1))
    aux_inds = np.arange(np.int(np.floor((nPads - 1) / 2)), 1 - 1, -1)
    freqT = np.concatenate([np.array([0.0]), freqT, -freqT[aux_inds - 1]])

    # We perform smoothing by multiplication in the Fourier domain
    for k in range(1, mtrx.shape[1 - 1] + 1):
        aux = np.exp(-0.25 * (scales[k - 1] * freqT) ** 2.0)
        smth = np.fft.ifft(np.multiply(aux, np.fft.fft(mtrx[k - 1, :], nPads)))
        mtrx[k - 1, :] = smth[1 - 1 : N]

    # We filter across several scales
    scalesToSmooth = 12
    flt = 1.0 / scalesToSmooth * np.ones((scalesToSmooth, 1))
    mtrx = ss.convolve2d(mtrx, flt, 'same')

    return mtrx

# The closest power of two such that two to the power is greater than, or equal to, the given value
def nextpow2(val):
    p = 0
    while 2 ** p < val:
        p = p + 1

    return p