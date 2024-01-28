#!/home/ubuntu/miniconda/bin/python
#
# This file is the entry point for testing the canonical coherence analysis methods provided by the functions
# 'getCanonicalCoherence()' and 'getCanonicalCoherenceW()'. Set the testID parameter to run a certain test.

from getCanonicalCoherence import getCanonicalCoherence
from getCanonicalCoherenceW import getCanonicalCoherenceW
import numpy as np
import matplotlib.pyplot as plt

## Preparing data
testID = 4      # Any integer between 1 and 4

if (testID == 1):                           # No coherence
    N = 3                                   # Number of variates in the vector (multivariate) time series

    fs = 1000.0                             # Discretisation frequency
    t = np.arange(0.0, 1.0, 1.0 / fs)

    fcommon = 200.0                         # Base frequency

    c = np.zeros((len(t), N))
    c[:, 1 - 1] = 1.0 * np.cos(2.0 * np.pi * t * fcommon)          # Signal No. 1's trend
    c[:, 2 - 1] = 0.5 * np.cos(2.0 * np.pi * t * fcommon / 2.0)    # Signal No. 2's trend
    c[:, 3 - 1] = 1.0 * np.cos(2.0 * np.pi * t * fcommon / 4.0)    # Signal No. 3's trend

    ddx = np.ndarray((len(t), N))
    for i in range(1, N + 1):
        ddx[:, i - 1] = c[:, i - 1] + np.random.randn(len(t))      # Each of the scalar time series has its own frequency plus a random noise
elif (testID == 2):                                         # Coherence at three different frequencies
    N = 2

    fs = 1000.0
    t = np.arange(0.0, 1.0, 1.0 / fs)

    fcommon = 200.0

    c1 = 1.0 * np.cos(2.0 * np.pi * t * fcommon)            # Trend No. 1 is sampled at the frequency 'fcommon'
    c2 = 1.0 * np.cos(2.0 * np.pi * t * fcommon / 4.0)      # Trend No. 2 is sampled at the frequency 'fcommon / 4'
    c = c1 + c2                                             # The total trend has two frequencies

    ddx = np.ndarray((len(t), N))
    for i in range(1, N + 1):
        ddx[:, i - 1] = c + np.random.randn(len(t))         # Each of the scalar time series within the multivariate time series contains the common total trend plus a random noise
elif (testID == 3):                                         # Coherence at two different frequencies in time-shifted time ranges
    N = 2

    fs = 1000.0
    t = np.arange(0.0, 2.0, 1.0 / fs)

    fcommon1 = 10.0
    fcommon2 = 50.0

    ddx = np.ndarray((len(t), N))
    ddx[:, 1 - 1] = np.multiply(np.cos(2.0 * np.pi * fcommon1 * t), np.logical_and(t >= 0.5, t < 1.1)) + np.multiply(np.cos(2.0 * np.pi * fcommon2 * t), np.logical_and(t >= 0.2, t < 1.4)) + 0.25 * np.random.randn(len(t))
    ddx[:, 2 - 1] = np.multiply(np.sin(2.0 * np.pi * fcommon1 * t), np.logical_and(t >= 0.6, t < 1.2)) + np.multiply(np.sin(2.0 * np.pi * fcommon2 * t), np.logical_and(t >= 0.4, t < 1.6)) + 0.35 * np.random.randn(len(t))
elif (testID == 4):
    N = 3

    fs = 1000.0
    t = np.arange(0.0, 2.0, 1.0 / fs)

    fcommon1 = 10.0
    fcommon2 = 50.0

    ddx = np.ndarray((len(t), N))
    ddx[:, 1 - 1] = np.multiply(np.cos(2.0 * np.pi * fcommon1 * t), np.logical_and(t >= 0.5, t < 1.1)) + np.multiply(np.cos(2.0 * np.pi * fcommon2 * t), np.logical_and(t >= 0.2, t < 1.4)) + 0.25 * np.random.randn(len(t))
    ddx[:, 2 - 1] = np.multiply(np.sin(2.0 * np.pi * fcommon1 * t), np.logical_and(t >= 0.6, t < 1.2)) + np.multiply(np.sin(2.0 * np.pi * fcommon2 * t), np.logical_and(t >= 0.4, t < 1.6)) + 0.35 * np.random.randn(len(t))
    ddx[:, 3 - 1] = np.multiply(np.sin(2.0 * np.pi * fcommon1 * t), np.logical_and(t >= 1.5, t < 1.8)) + np.multiply(np.sin(2.0 * np.pi * fcommon2 * t), np.logical_and(t >= 1.3, t < 1.7)) + 0.15 * np.random.randn(len(t))

waveletSigma = 6.0      # Default value is 6.0
energyLimit = 0.95

## Computing
[evt, ev, freq] = getCanonicalCoherence(ddx, fs, True)
[wevt, wev, freqS, coi, timeBorders] = getCanonicalCoherenceW(ddx, fs, np.array([t[int(np.floor(len(t) / 4) - 1)], t[int(np.floor(4 * len(t) / 5) - 1)]]), energyLimit, waveletSigma, True)

## Output
# (Fourier) CCA
plt.figure()
plt.subplot(2, 1, 1)
plt.semilogx(freq, evt, 'k-')
plt.xlabel(r'$f, \mathrm{Hz}$')
plt.ylabel(r'$C(f)$')
plt.grid('on')
plt.title('Canonical Coherence Analysis')

plt.subplot(2, 1, 2)
leg_txt = []
for i in range(1, N + 1):
    plt.semilogx(freq, ev[:, i - 1])
    leg_txt.append('%d' % (i))
plt.legend(leg_txt)
plt.xlabel(r'$f, \mathrm{Hz}$')
plt.ylabel(r'$c_i(f)$')
plt.grid('on')

# Wavelet CCA
plt.figure()
plt.pcolormesh(t, freqS, wevt)
#plt.xlabel(r'$t$')
plt.xlabel(r'$t - b$')
plt.ylabel(r'$f, \mathrm{Hz}$')
plt.title('Wavelet CCA, Total Coherence')
plt.yscale('log')
plt.plot(t, coi, 'w--')
plt.colorbar()
ax = plt.gca()
ax.autoscale(enable = True, axis = 'y', tight = True)

if ('timeBorders' in locals()):
    nb = timeBorders.shape[3 - 1]
    for ib in range(1, nb + 1):
        # We cut off the time moments (if any) which are out of the period of observation
        timeBorders[timeBorders[:, 1 - 1, ib - 1] < t[1 - 1], 1 - 1, ib - 1] = np.nan
        timeBorders[timeBorders[:, 2 - 1, ib - 1] > t[-1], 2 - 1, ib - 1] = np.nan
        plt.plot(timeBorders[:, 1 - 1, ib - 1], freqS, 'w--', timeBorders[:, 2 - 1, ib - 1], freqS, 'w--')

plt.figure()
for i in range(1, N + 1):
    plt.subplot(N, 1, i)
    plt.pcolormesh(t, freqS, wev[:, :, i - 1])
    if (i == 1):
        plt.title('Wavelet CCA, Partial Coherences')
    if (i == N):
        #plt.xlabel(r'$t$')
        plt.xlabel(r'$t - b$')
    plt.ylabel(r'$f, \mathrm{Hz}$')
    plt.yscale('log')
    plt.plot(t, coi, 'w--')
    plt.colorbar()
    ax = plt.gca()
    ax.autoscale(enable = True, axis = 'y', tight = True)

    if ('timeBorders' in locals()):
        for ib in range(1, nb + 1):
            plt.plot(timeBorders[:, 1 - 1, ib - 1], freqS, 'w--', timeBorders[:, 2 - 1, ib - 1], freqS, 'w--')

plt.show()
