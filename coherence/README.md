# Data-Science/coherence
This folder contains utilities implementing canonical coherence analysis (CCA) methods. To get familiar, run the corresponding `test.*` files.

### Multivariate CCAs
These functions realise multivariate CCAs when the number of scalar time series in a multivariate dataset can be any natural, not necessarily two. For more details on multivariate CCAs, please refer to [1-3].

### Technical issues
The Python CCA codes involve the function `getPowerSpectrum()` from `../utils`. The code of `getPowerSpectrum()` ensures an accurate computation of cross-spectra of multivariate datasets regardless of the length of time series. The function `getPowerSpectrum()` was developed once the standard Python function `scipy.signal.cpsd()` had been found to provide inaccurate results on long time series containing 10,000+ samples each.

### References
1. A.A. Lyubushin, *Data Analysis of Geophysical and Environmental Monitoring Systems*, Nauka, Moscow, 2008, available at https://search.rsl.ru/en/record/01003114864, also http://alexeylyubushin.narod.ru/Index.htm (in Russian).
2. A.A. Lyubushin, Synchronization of Geophysical Fields Fluctuations, in: T. Chelidze, L. Telesca and F. Vallianatos (eds.), *Complexity of Seismic Time Series: Measurement and Applications*, Elsevier, Amsterdam, 2018, pp. 161-197.
3. D.M. Filatov and A.A. Lyubushin, Precursory Analysis of GPS Time Series for Seismic Hazard Assessment, *Pure Appl. Geophys.*, 177 (2020) 509-530.
