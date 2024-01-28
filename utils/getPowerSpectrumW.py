#!/home/ubuntu/miniconda/bin/python
#
# This code implements the computation of the continuous wavelet auto- or cross- power spectrum. It is based on
# MATLAB's function 'cwtft()' (see also 'wcoherence()' there). The complex-valued Morlet wavelets of arbitrary time-
# frequency trade-off are employed. For details, please refer to [1-3].
#
# At the input:
#   - 'x' is a uni- or bivariate time series of the physical observables 'x_1(t)' and, optionally, 'x_2(t)'
#   - 'fs' is the sampling rate (optional)
#   - 'waveletSigma' is the Morlet wavelet's time-frequency resolution (optional)
#   - 'isSmoothing' is the parameter prescribing to smooth the spectrum (optional, recommended).
#
# If 'x_2(t)' is omitted then the auto-spectrum of 'x_1(t)' is computed; otherwise the cross-spectrum between 'x_1(t)'
# and 'x_2(t)' is computed. The time series are to be provided column-wise for each variate.
#
# At the output:
#   - 'Ps' is the power spectrum
#   - 'freqS' is the scale-related frequency range over which the spectrum is computed
#   - 'coi' is the cone of influence.
#
# NOTE:
# The implementation of the function 'smoothSpectrum()' should contain smoothings both in time and in scale-related
# frequency. However, because Python is inaccurate in the computation of two-dimensional convolution, the latter
# smoothing is not implemented. This may provoke overestimated coherence when 'getPowerSpectrumW()' is called from
# '../coherence/getCanonicalCoherenceW()'.
#
# REFERENCES:
# [1] J. Ashmead, Quanta, 1 (2012) 58-70.
# [2] C. Torrence and G.P. Compo, Bull. Am. Meteorol. Soc., 79 (1998) 61-78.
# [3] M.X. Cohen, Parameters of Morlet wavelet (time-frequency trade-off), https://www.youtube.com/watch?v=LMqTM7EYlqY
#
# The end user is granted perpetual permission to reproduce, adapt, and/or distribute this code, provided that
# an appropriate link is given to the original repository it was downloaded from.

import numpy as np
import scipy.signal as ss
from scipy.optimize import brentq

def getPowerSpectrumW(x, fs = 1.0, waveletSigma = 6.0, isSmoothing = True):
    # We precompute the scale-related frequency limits to subsequently compute the CWT within ...
    [preFreqs, _] = ss.welch(x[:, 1 - 1], fs)
    # ... and prepare a padding to later add at both sides of the signals
    x_len_adjusted = adjustSignalLength(x.shape[1 - 1])

    # We compute the continuous Fourier transform of the Morlet wavelet
    [wft, freqS, scales, coi] = MorletCWT(x.shape[1 - 1], fs, np.array([preFreqs[1 - 1], preFreqs[-1]]), x_len_adjusted, waveletSigma)

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
    # If smoothing is not performed, singular spectral matrices may appear in the function 'getCanonicalCoherenceW()'
    # in '../coherence/getCanonicalCoherenceW.py') due to nearly the same wavelet spectra
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
def MorletCWT(x_len, fs, preFreqs, x_adjustment, waveletSigma):
    # The list of dimensions of the quantities used:
    # [x_len] = [samples] = samples
    # [waveletSigma] = [omegaCutoff] = [freqT] = radians per second (angular frequencies)
    # [minScale] = [maxScale] = [scales] = relative units (dimensionless quantities)
    # [waveletFreq] = [preFreqs] = [minFreq] = [freqS] = [coi] = hertz, or cycles per second
    # [fs] = samples per second

    fpo = 12                                                    # Number of scale-related frequencies per octave. The greater this parameter, the more accurate the CWT
    a0 = 2.0 ** (1.0 / fpo)                                     # Scale ratio (or geometric increment) between two successive voices

    kappaSigma = np.exp(-0.5 * waveletSigma ** 2.0)
    cSigma = 1.0 / np.sqrt(1.0 + np.exp(-waveletSigma ** 2.0) - 2.0 * np.exp(-0.75 * waveletSigma ** 2.0))
    amplitude = cSigma / np.pi ** 0.25

    # We determine the cutoff angular frequency, in radians per second
    cutoffLevel = 0.1
    params = (cutoffLevel, amplitude, waveletSigma, kappaSigma)

    omegaMax = np.sqrt(2.0 * 750.0) + waveletSigma
    if (waveletFourierTransform(waveletSigma, *params) > 0.0):
        omegaCutoff = omegaMax
    else:
        omegaCutoff = brentq(waveletFourierTransform, waveletSigma, omegaMax, args = params)

    # We make a correction for the preset scale-related frequency (if relevant)
    minScale = omegaCutoff / np.pi / 1.0                        # Preliminary minimum admissible scale, in relative units (radians per second / radians per cycle / cycles per second)
    maxScale = x_len / (2.0 * np.sqrt(2.0))                     # Preliminary maximum admissible scale, in relative units (samples / samples)

    if (maxScale < minScale * a0):
        maxScale = minScale * a0

    waveletFreq = waveletSigma / (2.0 * np.pi)                  # The frequency of the sine harmonic of the Morlet wavelet, a constant, in hertz (or in cycles per second)
    minFreq = waveletFreq * 1.0 / maxScale * fs                 # Minimum admissible scale-related frequency, in hertz

    # So, we make the correction
    if (preFreqs[1 - 1] < minFreq):
        preFreqs[1 - 1] = minFreq                               # The corrected preset scale-related frequency, in hertz

    # Now we prepare the time-related frequencies (usually denoted in scientific literature as 'omega'), in radians per second ...
    freqT = makeFourierTransformGrid(x_len, x_adjustment)

    nFreq = len(freqT)

    # ... and scales, in relative units
    preFreqs = 2.0 * np.pi * preFreqs / fs                      # We convert the scale-related frequencies from hertz to radians per sample
    minScale = waveletSigma / preFreqs[2 - 1] / 1.0             # The minimum and ...
    maxScale = waveletSigma / preFreqs[1 - 1] / 1.0             # ... maximum admissible scales are updated but keep the same dimensions: relative units

    maxNumOctaves = np.log2(maxScale / minScale)
    scales = minScale * a0 ** np.arange(0, int(np.floor(maxNumOctaves * fpo)) + 1)       # Scales, in relative units

    nScales = len(scales)

    # Finally, we compute the CWT of the Morlet wavelet on the determined time-related frequencies and scales
    wft = np.zeros((nScales, nFreq))
    for j in range(1, nScales + 1):
        eexp1 = np.multiply(-(scales[j - 1] * freqT - waveletSigma) ** 2.0 / 2.0, freqT > 0.0)
        eexp2 = np.multiply(-(scales[j - 1] * freqT) ** 2.0 / 2.0, freqT > 0.0)
        wft[j - 1, :] = amplitude * np.sqrt(scales[j - 1]) * np.multiply(np.exp(eexp1) - kappaSigma * np.exp(eexp2), freqT > 0.0)

    # As the last action, we compute the actual scale-related frequencies ...
    freqS = waveletFreq * 1.0 / scales * fs                     # Actual scale-related frequencies, in hertz

    # ... and the cone of influence
    xl = int(np.ceil(x_len / 2.0))
    if (np.mod(x_len, 2) == 1):
        samples = np.concatenate((np.arange(1, xl + 1), np.arange(xl - 1, 1 - 1, -1)))
    else:
        samples = np.concatenate((np.arange(1, xl + 1), np.arange(xl, 1 - 1, -1)))

    coi = waveletFreq * np.sqrt(2.0) / samples * fs * 1.0       # In hertz

    # We truncate the excessive values of the cone of influence
    coi[coi > freqS[1 - 1]] = freqS[1 - 1]
    coi[coi < freqS[-1]] = freqS[-1]

    return [wft, freqS, scales, coi]

# Fourier transform of the Morlet wavelet minus the cutoff value
def waveletFourierTransform(omega, *params):
    [cutoffLevel, amplitude, waveletSigma, kappaSigma] = params

    root = cutoffLevel - amplitude * (np.exp(-0.5 * (omega - waveletSigma) ** 2.0) - kappaSigma * np.exp(-0.5 * omega ** 2.0))

    return root

# Adjusting signal length
def adjustSignalLength(x_len):
    if (x_len <= 100000):
        x_len_adjusted = int(np.floor(x_len / 2.0))
    else:
        x_len_adjusted = int(np.ceil(np.log2(x_len)))

    return x_len_adjusted

# Making time-related frequency grid for performing the Fourier transforms
def makeFourierTransformGrid(x_len, x_adjustment):
    # This is the size of the time-related frequency range
    N = x_len + 2 * x_adjustment

    # This is the first half of the frequency range ...
    freqT = np.arange(1, int(np.floor(N / 2.0) + 1))
    freqT = freqT * (2.0 * np.pi) / N
    # ... while these are the indices for covering the second half
    aux_inds = np.arange(int(np.floor((N - 1) / 2)), 1 - 1, -1)

    # So, we make the frequency range, in radians per second
    freqT = np.concatenate([np.array([0.0]), freqT, -freqT[aux_inds - 1]])

    return freqT

# Smoothing the spectrum
def smoothSpectrum(mtrx, scales):
    # We are to filter along the time axis ...
    x_len = mtrx.shape[2 - 1]

    nfft = 2 ** nextpow2(x_len)                                                     # We enlarge the number of samples to the closest power of two

    freqT = 2.0 * np.pi / nfft * np.arange(1, int(np.floor(nfft / 2) + 1))
    aux_inds = np.arange(int(np.floor((nfft - 1) / 2)), 1 - 1, -1)
    freqT = np.concatenate([np.array([0.0]), freqT, -freqT[aux_inds - 1]])          # In radians per second

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
