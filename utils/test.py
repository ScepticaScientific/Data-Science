#!/home/ubuntu/miniconda/bin/python
#
# This file is the entry point for testing the computation of power spectrum ('getPowerSpectrum.py').
# The time series are simulated below.

from getPowerSpectrum import getPowerSpectrum
import numpy as np
import matplotlib.pyplot as plt

## Preparing data
N = 10000       # Number of samples

x = np.random.randn(N)
y = np.random.randn(N)

## Computing
[Pxx, freq] = getPowerSpectrum(y)           # Auto-spectrum
#[Pxy, freq] = getPowerSpectrum(x, y)       # Cross-spectrum

## Output
plt.figure()
plt.semilogx(freq, Pxx, 'k-')
#plt.semilogx(freq, np.abs(Pxy), 'k-')
plt.xlabel(r'$\omega$')
plt.ylabel(r'$E(\omega)$')
plt.grid('on')
plt.title('Power Spectrum')

plt.show()
