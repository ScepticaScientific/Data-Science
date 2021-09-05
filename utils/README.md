# Data-Science/utils
This folder contains various auxiliary functions used in data science computations. To get familiar, run the corresponding `test.*` files.

### Power spectrum algorithm
The Python function `getPowerSpectrum()` implements the non-parametric Welch method for computing Fourier auto- and cross- power spectra. It was developed once the standard Python function `scipy.signal.csd()` had been found to provide inaccurate results on long time series containing 10,000+ samples.

### Others
The functions `getPowerSpectrumW()` implement the computation of wavelet auto- and cross- power spectra using complex-valued Morlet wavelets of arbitrary time-frequency resolution.