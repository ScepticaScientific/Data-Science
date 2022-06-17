# Data-Science/coherence
This folder contains utilities implementing canonical coherence analysis (CCA) methods. To get familiar, run the corresponding `test.*` files.

### Multivariate CCAs
These functions realise multivariate Fourier- and wavelet-based CCAs when the number of scalar time series in a multivariate dataset is any natural, not necessarily two. For wavelet-based CCA complex-valued Morlet wavelets of arbitrary, user-specified time-frequency resolution are used. For more details on multivariate CCAs, please refer to [1-3].

### Technical issues
The Python Fourier CCA code involves calling the function `getPowerSpectrum()` from `../utils`. The code of `getPowerSpectrum()` ensures an accurate computation of the power spectra of multivariate datasets regardless of the length of the time series. The function `getPowerSpectrum()` was developed once the standard Python function `scipy.signal.csd()` had been found to provide inaccurate results on long time series containing 10,000+ samples each.

### References
1. A.A. Lyubushin, *Data Analysis of Systems of Geophysical and Ecological Monitoring*, Nauka, Moscow, 2007, available at https://search.rsl.ru/en/record/01003114864, also http://alexeylyubushin.narod.ru/Index.htm (in Russian).
2. A.A. Lyubushin, Synchronization of Geophysical Fields Fluctuations, in: T. Chelidze, L. Telesca and F. Vallianatos (eds.), *Complexity of Seismic Time Series: Measurement and Applications*, Elsevier, Amsterdam, 2018, pp. 161-197. <p><a href = "https://doi.org/10.1016/B978-0-12-813138-1.00006-7" rel = "nofollow"><img src = "https://zenodo.org/badge/DOI/10.1016/B978-0-12-813138-1.00006-7.svg" alt = "DOI:10.1016/B978-0-12-813138-1.00006-7" style = "vertical-align: top; max-width: 100%;"></a></p>
3. D.M. Filatov and A.A. Lyubushin, Stochastic Dynamical Systems Always Undergo Trending Mechanisms of Transition to Criticality, *Physica A: Stat. Mech. Appl.*, 527 (2019) 121309. <p><a href = "https://doi.org/10.1016/j.physa.2019.121309" rel = "nofollow"><img src = "https://zenodo.org/badge/DOI/10.1016/j.physa.2019.121309.svg" alt = "DOI:10.1016/j.physa.2019.121309" style = "vertical-align: top; max-width: 100%;"></a></p>
4. D.M. Filatov, On Spatial Synchronisation as a Manifestation of Irregular Energy Cascades in Continuous Media under the Transition to Criticality, *Nonlinear Dyn.*, (2022), in press. <p><a href = "http://dx.doi.org/10.1007/s11071-022-07580-7" rel = "nofollow"><img src = "https://zenodo.org/badge/DOI/10.1007/s11071-022-07580-7.svg" alt = "DOI:10.1007/s11071-022-07580-7" style = "vertical-align: top; max-width: 100%;"></a></p>
