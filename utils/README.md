# Data-Science/utils
This folder contains various auxiliary functions used in data science computations. To get familiar, run `test.py`.

### Power spectrum algorithm
The function `getPowerSpectrum()` implements the non-parametric Welch method for computing power auto- and cross-spectra. It was developed once the standard Python function `scipy.signal.csd()` had been found to provide inaccurate results on long time series containing 10,000+ samples.