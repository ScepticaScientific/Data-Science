#!/home/ubuntu/miniconda/bin/python
#
# This code implements the computation of the continuous wavelet auto- or cross- power spectrum. It is based on
# MATLAB's function cwtft() (see also wcoherence() there). The Morlet wavelet of the time-frequency trade-off
# 'sigma = 6' is employed. For details, please refer to [1, 2].
#
# At the input, 'x' is a uni- or bivariate time series of the physical observables 'x_1(t)' and, possibly, 'x_2(t)',
# 'fs' is the sampling rate (optional), while 'isSmoothing' is a parameter prescribing to smooth the spectrum (optional,
# recommended). If 'x_2(t)' is omitted then the auto-spectrum of 'x_1(t)' is computed; otherwise the cross-spectrum
# between 'x_1(t)' and 'x_2(t)' is computed. The time series are to be provided column-wise for each variate.
#
# At the output, 'Ps' is the power spectrum, 'freqS' is the scale-related frequency range over which the spectrum is
# computed, while 'coi' is the cone of influence.
#
# NOTE:
# The function smoothSpectrum() should contain smoothings both in time and in scale-related frequency. However, because
# Python is inaccurate in the computation of two-dimensional convolution, the latter smoothing is not implemented. This
# may provoke overestimated coherence when getPowerSpectrumW() is called from ../coherence/getCanonicalCoherenceW().
#
# REFERENCES:
# [1] C. Torrence and G.P. Compo, Bull. Am. Meteorol. Soc., 79 (1998) 61-78.
# [2] Michael X. Cohen, Parameters of Morlet wavelet, https://www.youtube.com/watch?v=LMqTM7EYlqY
#
# The end user is granted perpetual permission to reproduce, adapt, and/or distribute this code, provided that
# an appropriate link is given to the original repository it was downloaded from.

import numpy as np
import scipy.signal as ss

def getPowerSpectrumW(x, fs = 1.0, isSmoothing = True):
    # We precompute the scale-related frequency limits to subsequently compute the CWT within ...
    [preFreqs, _] = ss.welch(x[:, 1 - 1], fs)
    # ... and prepare a padding to later add at both sides of the signals
    x_len_adjusted = adjustSignalLength(x.shape[1 - 1])

    # We compute the continuous Fourier transform of the Morlet wavelet
    [wft, freqS, scales, coi] = MorletCWT(x.shape[1 - 1], fs, np.array([preFreqs[1 - 1], preFreqs[-1]]), x_len_adjusted)

    # Now we compute the CWT of the first time series ...
    if (len(x.shape) == 1):
        cwtX = getCWT(x, x_len_adjusted, wft)
    else:
        cwtX = getCWT(x[:, 1 - 1], x_len_adjusted, wft)

    if (len(x.shape) == 1):
        cwtY = cwtX
    else:
        if (x.shape[2 - 1] == 2):
            cwtY = getCWT(x[:, 2 - 1], x_len_adjusted, wft)
        else:
            cwtY = cwtX
    # ... and then multiply it with the corresponding conjugate
    Ps = np.multiply(cwtX, np.conj(cwtY))

    # Optional smoothing. This is recommended for removing noise which appears after inverse Fourier transform.
    # If smoothing is not performed, singular spectral matrices may appear in function getCanonicalCoherenceW()
    # in file '../coherence/getCanonicalCoherenceW.py') due to nearly the same wavelet spectra
    if (isSmoothing == True):
        Ps = smoothSpectrum(Ps, scales)

    return [Ps, freqS, coi]

## Auxiliaries
# Computation of the continuous wavelet transform of the data
def getCWT(x, x_len_adjusted, wft):
    x_len = x.shape[1 - 1]

    # We adjust the signal
    auxIndLeft = np.arange(x_len_adjusted, 1 - 1, -1)
    auxIndRight = np.arange(x_len, x_len - x_len_adjusted + 1 - 1, -1)
    xAdjusted = np.concatenate((np.conj(x[auxIndLeft - 1]), x, np.conj(x[auxIndRight - 1])))

    # We compute the continuous Fourier transform of the adjusted time series ...
    xAdjustedFFT = np.fft.fft(xAdjusted)
    # ... and compute the continuous wavelet transform of the adjusted time series via inverse Fourier transform
    cwt = np.fft.ifft(np.multiply(np.tile(xAdjustedFFT, (wft.shape[1 - 1], 1)), wft), axis = 2 - 1)

    # Finally we keep only those values that correspond to the original time series' frequencies
    cwt = cwt[:, x_len_adjusted + 1 - 1 : x_len_adjusted + x_len]

    return cwt

# Computing the CWT of Morlet wavelet
def MorletCWT(x_len, fs, preFreqs, x_adjustment):
    # Before computing the wavelet transform, we make a correction for the minimum admissible scale-related frequency
    sigma = 6.0                                 # The Morlet wavelet's time-frequency trade-off, in radians
    FourierFactor = 2.0 * np.pi / sigma         # A constant, in relative units

    centralFrequency = 8.14597      # The global maximum of the Fourier transform of the Morlet wavelet under 'sigma = 6' and a slight shift -0.075 in the ordinate

    minScale = 1 * centralFrequency / np.pi         # Preliminary minimum admissible scale, in samples
    maxScale = x_len / (2.0 * np.sqrt(2.0))         # Preliminary maximum admissible scale, in samples

    fpo = 12                        # Number of scale-related frequencies per octave. The greater this parameter, the more accurate the CWT
    a0 = 2.0 ** (1.0 / fpo)         # Scale ratio (or geometric increment) between two successive voices

    if (maxScale < minScale * a0):
        maxScale = minScale * a0

    minFreq = 1.0 / (maxScale * FourierFactor) * fs # Minimum admissible scale-related frequency, in hertz

    # We make the correction, if required
    if (preFreqs[1 - 1] < minFreq):
        preFreqs[1 - 1] = minFreq                   # The corrected minimum admissible scale-related frequency, in hertz

    # Now we prepare the time-related frequencies (in radians) ...
    freqT = makeFourierTransformGrid(x_len, x_adjustment)

    nFreq = len(freqT)

    # ... and scales (in samples)
    preFreqs = 2.0 * np.pi * preFreqs / fs  # We convert the scale-related frequencies from hertz to radians per sample
    minScale = sigma / preFreqs[2 - 1]      # The minimum and ...
    maxScale = sigma / preFreqs[1 - 1]      # ... maximum admissible scales are updated but keep the same dimensions: samples

    maxNumOctaves = np.log2(maxScale / minScale)
    scales = minScale * a0 ** np.arange(0, np.int(np.floor(maxNumOctaves * fpo)) + 1)        # Scales, in samples

    nScales = len(scales)

    # Finally, we compute the CWT of the Morlet wavelet
    wft = np.zeros((nScales, nFreq))
    amplitude = 1.0 / np.pi ** 0.25
    for j in range(1, nScales + 1):
        eexp = np.multiply(-(scales[j - 1] * freqT - sigma) ** 2.0 / 2.0, freqT > 0.0)
        wft[j - 1, :] = amplitude * np.sqrt(scales[j - 1]) * np.multiply(np.exp(eexp), freqT > 0.0)

    freqS = 1.0 / FourierFactor / scales * fs       # Scale-related frequencies, in hertz

    # We compute the cone of influence as well
    xl = np.int(np.ceil(x_len / 2.0))
    if (np.mod(x_len, 2) == 1):
        samples = np.concatenate((np.arange(1, xl + 1), np.arange(xl - 1, 1 - 1, -1)))
    else:
        samples = np.concatenate((np.arange(1, xl + 1), np.arange(xl, 1 - 1, -1)))

    coi = FourierFactor / np.sqrt(2.0) * samples / fs       # In seconds
    coi = 1.0 / coi                                         # In hertz

    # We truncate the excessive values of the cone of influence
    coi[coi > freqS[1 - 1]] = freqS[1 - 1]
    coi[coi < freqS[-1]] = freqS[-1]

    return [wft, freqS, scales, coi]

# Adjusting signal length
def adjustSignalLength(x_len):
    if (x_len <= 100000):
        x_len_adjusted = np.int(np.floor(x_len / 2.0))
    else:
        x_len_adjusted = np.int(np.ceil(np.log2(x_len)))

    return x_len_adjusted

# Making time-related frequency grid for performing the Fourier transforms
def makeFourierTransformGrid(x_len, x_adjustment):
    # This is the size of the time-related frequency range
    N = x_len + 2 * x_adjustment

    # This is the first half of the frequency range ...
    freqT = np.arange(1, np.int(np.floor(N / 2.0) + 1))
    freqT = freqT * (2.0 * np.pi) / N
    # ... while these are the indices for covering the second half
    aux_inds = np.arange(np.int(np.floor((N - 1) / 2)), 1 - 1, -1)

    # So, we make the frequency range in radians
    freqT = np.concatenate([np.array([0.0]), freqT, -freqT[aux_inds - 1]])

    return freqT

# Smoothing the spectrum
def smoothSpectrum(mtrx, scales):
    # We are to filter along the time axis ...
    x_len = mtrx.shape[2 - 1]

    nfft = 2 ** nextpow2(x_len)    # We enlarge the number of samples to the closest power of two

    freqT = 2.0 * np.pi / nfft * np.arange(1, np.int(np.floor(nfft / 2) + 1))
    aux_inds = np.arange(np.int(np.floor((nfft - 1) / 2)), 1 - 1, -1)
    freqT = np.concatenate([np.array([0.0]), freqT, -freqT[aux_inds - 1]])          # In radians

    for k in range(1, mtrx.shape[1 - 1] + 1):
        aux = np.exp(-0.5 * (scales[k - 1] * freqT) ** 2.0)
        smth = np.fft.ifft(np.multiply(aux, np.fft.fft(mtrx[k - 1, :], nfft)))
        mtrx[k - 1, :] = smth[1 - 1 : x_len]

    # ... and then across the scales
    ## Python is inaccurate in the computation of two-dimensional convolution, so this smoothing is dismissed
    #scalesToSmooth = 12      # The greater this parameter, the less accurate the spectrum
    #flt = 1.0 / scalesToSmooth * np.ones((scalesToSmooth, 1))
    #mtrx = ss.convolve2d(np.real(mtrx), flt, 'same')

    return mtrx

# The closest power of two such that two to the power is greater than, or equal to, the given value
def nextpow2(val):
    p = 0
    while 2 ** p < val:
        p = p + 1

    return p
