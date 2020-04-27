# Data-Science/utils
This folder contains various auxiliary functions used in data science computations. To get familiar, run `test.py`.

### Power spectrum algorithm
The function `getPowerSpectrum()` implements the non-parametric Welch method for computing Fourier auto- and cross- power spectra. It was developed once the standard Python function `scipy.signal.csd()` had been found to provide inaccurate results on long time series containing 10,000+ samples. The function `getPowerSpectrumW()` implements the computation of wavelet auto- and cross- power spectra using the Morlet wavelet.