#!/home/ubuntu/miniconda/bin/python
#
# This code implements the computation of the continuous wavelet auto- or cross- power spectrum. It is based on
# MATLAB's function cwtft() (see also wcoherence() there). For details, please refer to [1].
#
# At the input, 'x' is a uni- or bivariate time series of the physical observables 'x_1(t)' and, possibly, 'x_2(t)';
# 'fs' is the sampling rate; 'isSmoothing' is an optional parameter prescribing to smooth the spectrum time-wise
# over several scales. If 'x_2(t)' is omitted then the auto-spectrum of 'x_1(t)' is computed; otherwise the cross-
# spectrum between 'x_1(t)' and 'x_2(t)' is computed. The time series are to be provided column-wise for each variate.
#
# At the output, 'Ps' is the power spectrum, while 'freqS' is the scale-related frequency range.
#
# REFERENCES:
# [1] C. Torrence and G.P. Compo, Bull. Am. Meteorol. Soc., 79 (1998) 61-78.
#
# The end user is granted perpetual permission to reproduce, adapt, and/or distribute this code, provided that
# an appropriate link is given to the original repository it was downloaded from.

import numpy as np
import scipy.signal as ss

def getPowerSpectrumW(x, fs = 1.0, isSmoothing = False):
    # We compute the continuous wavelet transform of the first time series ...
    if (len(x.shape) == 1):
        [cwtX, freqS] = getCWT(x, fs)
    else:
        [cwtX, freqS] = getCWT(x[:, 1 - 1], fs)

    if (len(x.shape) == 1):
        cwtY = cwtX
    else:
        if (x.shape[2 - 1] == 2):
            [cwtY, _] = getCWT(x[:, 2 - 1], fs)
        else:
            cwtY = cwtX
    # ... and then multiply it with the corresponding conjugate
    Ps = np.multiply(cwtX, np.conj(cwtY))         # Inner product of 'cwtX' with its conjugate

    # Optional smoothing
    if (isSmoothing == True):
        Ps = smoothSpectrum(Ps)

    return [Ps, freqS]

## Auxiliaries
# Computation of continuous wavelet transform
def getCWT(x, fs):
    # Initialisation
    x_len = len(x)

    # We compute the continuous Fourier transform of the data, ...
    xFFT = np.fft.fft(x)
    # ... prepare the continuous Fourier transform of the Morlet wavelet ...
    [wft, freqS] = MorletCWT(x_len, fs)
    # ... and compute the continuous wavelet transform of the data via the inverse Fourier transform
    cwt = np.fft.ifft(np.multiply(np.tile(xFFT, (len(freqS), 1)), wft), axis = 2 - 1)

    return [cwt, freqS]

def MorletCWT(x_len, fs):
    freqT = 2.0 * np.pi * fs / x_len * np.arange(1, np.int(np.floor(x_len / 2) + 1))
    aux_inds = np.arange(np.int(np.floor((x_len - 1) / 2)), 1 - 1, -1)
    freqT = np.concatenate([np.array([0.0]), freqT, -freqT[aux_inds - 1]])

    nFreq = len(freqT)

    fpo = 12    # Frequencies per octave. The greater this parameter, the more accurate the CWT
    a0 = 2.0 ** (1.0 / fpo)
    maxNumOctaves = np.int(np.floor(np.log2(x_len)) - 1)
    scales = 2.0 / fs * a0 ** np.arange(0, maxNumOctaves * fpo + 1)

    nScales = len(scales)

    wft = np.zeros((nScales, nFreq))

    amplitude = 1.0 / np.pi ** 0.25 * np.sqrt(fs)
    for j in range(1, nScales + 1):
        eexp = np.multiply(-(scales[j - 1] * freqT - 6.0) ** 2.0 / 2.0, freqT > 0.0)
        wft[j - 1, :] = amplitude * np.sqrt(scales[j - 1]) * np.multiply(np.exp(eexp), freqT > 0.0)

    freqS = 3.0 / np.pi / scales

    return [wft, freqS]

# Smoothing the spectrum
def smoothSpectrum(mtrx):
    # We filter across several scales
    scalesToSmooth = 5      # The greater this parameter, the less accurate the spectrum
    flt = 1.0 / scalesToSmooth * np.ones((scalesToSmooth, 1))
    mtrx = ss.convolve2d(mtrx, flt, 'same')

    return mtrx