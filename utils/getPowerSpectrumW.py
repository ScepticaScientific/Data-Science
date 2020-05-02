#!/home/ubuntu/miniconda/bin/python
#
# This code implements the computation of the continuous wavelet auto- or cross- power spectrum. It is based on
# MATLAB's function cwtft() (see also wcoherence() there). A Morlet wavelet is employed. For details, please refer
# to [1].
#
# At the input, 'x' is a uni- or bivariate time series of the physical observables 'x_1(t)' and, possibly, 'x_2(t)',
# 'fs' is the sampling rate, while 'isSmoothing' is an optional parameter prescribing to smooth the spectrum.
# If 'x_2(t)' is omitted then the auto-spectrum of 'x_1(t)' is computed; otherwise the cross-spectrum between 'x_1(t)'
# and 'x_2(t)' is computed. The time series are to be provided column-wise for each variate.
#
# At the output, 'Ps' is the power spectrum, while 'freqS' is the scale-related frequency range over which
# the spectrum is computed.
#
# REFERENCES:
# [1] C. Torrence and G.P. Compo, Bull. Am. Meteorol. Soc., 79 (1998) 61-78.
#
# The end user is granted perpetual permission to reproduce, adapt, and/or distribute this code, provided that
# an appropriate link is given to the original repository it was downloaded from.

import numpy as np
import scipy.signal as ss

def getPowerSpectrumW(x, fs = 1.0, isSmoothing = False):
    # We preliminarily compute the continuous Fourier transform of the Morlet wavelet
    # used in the subsequent calculations
    adjustedP2 = np.int(1 + np.floor(np.log2(x.shape[1 - 1]) + 0.4999))
    xAdjusted_len = 2 ** adjustedP2
    [wft, freqS, scales] = MorletCWT(xAdjusted_len, x.shape[1 - 1], fs)

    # We compute the continuous wavelet transform of the first time series ...
    if (len(x.shape) == 1):
        cwtX = getCWT(x, xAdjusted_len, wft)
    else:
        cwtX = getCWT(x[:, 1 - 1], xAdjusted_len, wft)

    if (len(x.shape) == 1):
        cwtY = cwtX
    else:
        if (x.shape[2 - 1] == 2):
            cwtY = getCWT(x[:, 2 - 1], xAdjusted_len, wft)
        else:
            cwtY = cwtX
    # ... and then multiply it with the corresponding conjugate
    Ps = np.multiply(cwtX, np.conj(cwtY))

    # Optional smoothing. This is recommended for removing noise appearing due to using inverse Fourier transform
    if (isSmoothing == True):
        Ps = smoothSpectrum(Ps, scales, fs)

    return [Ps, freqS]

## Auxiliaries
# Computation of the continuous wavelet transform of the data
def getCWT(x, xAdjusted_len, wft):
    # Initialisation
    x_len = len(x)

    # We centralise the original time series, expand it to a power of two and then symmetrise
    xAdjusted = symmetrise(x - np.mean(x), xAdjusted_len - x_len)

    # We compute the continuous Fourier transform of the adjusted time series ...
    xAdjustedFFT = np.fft.fft(xAdjusted)
    # ... and compute the continuous wavelet transform of the adjusted time series via inverse Fourier transform
    cwt = np.fft.ifft(np.multiply(np.tile(xAdjustedFFT, (wft.shape[1 - 1], 1)), wft), axis = 2 - 1)

    # Finally we keep only those values that correspond to the original time series' frequencies
    cwt = cwt[:, 1 - 1 : x_len]

    return cwt

# Expansion and symmetrisation of time series
def symmetrise(x, x_len):
    xl = len(x)

    if (xl == 1):
        return np.tile(x, (1, x_len + 1)).transpose()

    while (x_len + 1 > xl):
        x_len = x_len - (xl - 1)
        x = symmetrise(x, xl - 1)
        xl = len(x)

    inds = np.concatenate([np.arange(1, xl + 1), np.arange(xl - 1, xl - x_len - 1, -1)])

    return x[inds - 1]

# Computing the CWT of a Morlet wavelet
def MorletCWT(xAdjusted_len, x_len, fs):
    # We prepare the time-related frequencies ...
    freqT = 2.0 * np.pi * fs / xAdjusted_len * np.arange(1, np.int(np.floor(xAdjusted_len / 2) + 1))
    aux_inds = np.arange(np.int(np.floor((xAdjusted_len - 1) / 2)), 1 - 1, -1)
    freqT = np.concatenate([np.array([0.0]), freqT, -freqT[aux_inds - 1]])

    nFreq = len(freqT)

    # ... and scales
    fpo = 12    # Frequencies per octave. The greater this parameter, the more accurate the CWT
    a0 = 2.0 ** (1.0 / fpo)
    maxNumOctaves = np.int(np.floor(np.log2(x_len)) - 1)
    scales = 2.0 / fs * a0 ** np.arange(0, maxNumOctaves * fpo + 1)

    nScales = len(scales)

    # Now we compute the CWT of a Morlet wavelet
    wft = np.zeros((nScales, nFreq))
    amplitude = 1.0 / np.pi ** 0.25 * np.sqrt(fs)
    for j in range(1, nScales + 1):
        eexp = np.multiply(-(scales[j - 1] * freqT - 6.0) ** 2.0 / 2.0, freqT > 0.0)
        wft[j - 1, :] = amplitude * np.sqrt(scales[j - 1]) * np.multiply(np.exp(eexp), freqT > 0.0)

    freqS = 3.0 / np.pi / scales

    return [wft, freqS, scales]

# Smoothing the spectrum
def smoothSpectrum(mtrx, scales, fs):
    # We filter along the time axis ...
    x_len = mtrx.shape[2 - 1]

    nfft = 2 ** nextpow2(x_len)    # We enlarge the number of samples to the closest power of two

    freqT = 2.0 * np.pi * fs / nfft * np.arange(1, np.int(np.floor(nfft / 2) + 1))
    aux_inds = np.arange(np.int(np.floor((nfft - 1) / 2)), 1 - 1, -1)
    freqT = np.concatenate([np.array([0.0]), freqT, -freqT[aux_inds - 1]])

    for k in range(1, mtrx.shape[1 - 1] + 1):
        aux = np.exp(-0.25 * (scales[k - 1] * freqT) ** 2.0)
        smth = np.fft.ifft(np.multiply(aux, np.fft.fft(mtrx[k - 1, :], nfft)))
        mtrx[k - 1, :] = smth[1 - 1 : x_len]

    # ... and then across the scales
    scalesToSmooth = 12      # The greater this parameter, the less accurate the spectrum
    flt = 1.0 / scalesToSmooth * np.ones((scalesToSmooth, 1))
    mtrx = ss.convolve2d(mtrx, flt, 'same')

    return mtrx

# The closest power of two such that two to the power is greater than, or equal to, the given value
def nextpow2(val):
    p = 0
    while 2 ** p < val:
        p = p + 1

    return p
