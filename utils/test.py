#!/home/ubuntu/miniconda/bin/python
#
# This file is the entry point for testing the computation of continuous Fourier and wavelet power spectra
# provided by the functions getPowerSpectrum() and getPowerSpectrumW().

from getPowerSpectrum import getPowerSpectrum
from getPowerSpectrumW import getPowerSpectrumW
import numpy as np
import matplotlib.pyplot as plt

## Preparing data
fs = 1000.0                             # Discretisation frequency
t = np.arange(0.0, 1.0, 1.0 / fs)
fcommon = 100.0                         # Base frequency
c = np.cos(2.0 * np.pi * t * fcommon)

N = 2                                   # Number of variates (1 or 2)
x = np.ndarray((len(t), N))
for i in range(1, N + 1):
    x[:, i - 1] = c + np.random.randn(len(t))

## Computing
# (Fourier) power spectra
[Ps, freq] = getPowerSpectrum(x, fs)

# Wavelet power spectra
[wPs, freqS, coi] = getPowerSpectrumW(x, fs, True)

## Output
# (Fourier) power spectra
plt.figure()
plt.semilogx(freq, np.abs(Ps), 'k-')        # We plot the absolute value of the spectrum, i.e. its magnitude, or energy, because, if x_2 != x_1, it is complex-valued
plt.xlabel(r'$\omega$')
plt.ylabel(r'$E(\omega)$')
plt.grid('on')
plt.title('Power Spectrum')

# Wavelet power spectra
plt.figure()
plt.pcolormesh(t, freqS, np.abs(wPs))       # Similarly to the Fourier case, we plot the absolute value of the wavelet spectrum
plt.xlabel(r'$t$')
plt.ylabel(r'$\omega$')
plt.title('Wavelet Power Spectrum')
plt.yscale('log')
plt.plot(t, coi, 'w--')
plt.colorbar()
ax = plt.gca()
ax.autoscale(enable = True, axis = 'y', tight = True)

plt.show()
