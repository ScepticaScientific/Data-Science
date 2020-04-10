#!/home/ubuntu/miniconda/bin/python
#
# This file is the entry point for testing the canonical coherence analysis ('getCanonicalCoherence.py').
# The initial multivariate dataset is simulated below.

from getCanonicalCoherence import getCanonicalCoherence
import numpy as np
import matplotlib.pyplot as plt

## Preparing data
N = 3                                   # Number of variates in the vector (multivariate) time series
fs = 1000.0                             # Discretisation frequency
t = np.arange(0.0, 1.0, 1.0 / fs)

fcommon = 200.0                         # Base frequency

c1 = 1.0 * np.cos(2.0 * np.pi * t * fcommon)            # Trend No. 1 is sampled at the frequency 'fcommon'
c2 = 0.5 * np.cos(2.0 * np.pi * t * fcommon / 2.0)      # Trend No. 2 is sampled at the frequency 'fcommon / 2'
c3 = 1.0 * np.cos(2.0 * np.pi * t * fcommon / 4.0)      # Trend No. 3 is sampled at the frequency 'fcommon / 4'
c = c1 + c2 + c3                                        # The total trend has three frequencies

ddx = np.ndarray((len(t), N))
for i in range(1, N + 1):
    # Each of the scalar time series within the multivariate time series contains the total trend plus a random noise
    ddx[:, i - 1] = c + np.random.randn(len(t))

## Computing
[evt, ev, freq] = getCanonicalCoherence(ddx)

## Output
plt.figure()
plt.subplot(2, 1, 1)
plt.semilogx(freq * fs, evt, 'k-')
plt.xlabel(r'$\omega$')
plt.ylabel(r'$C(\omega)$')
plt.grid('on')
plt.title('Canonical Coherence Analysis')

plt.subplot(2, 1, 2)
leg_txt = []
for i in range(1, N + 1):
    plt.semilogx(freq * fs, ev[:, i - 1])
    leg_txt.append('%d' % (i))
plt.legend(leg_txt)
plt.xlabel(r'$\omega$')
plt.ylabel(r'$c_i(\omega)$')
plt.grid('on')

plt.show()
