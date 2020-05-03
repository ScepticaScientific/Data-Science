# Data-Science/fractals
This folder contains utilities implementing various fractal analysis methods. To get familiar, run `test.py`.

### Fractal analysis algorithms
The functions `getHurstByUpscaling()` and `getMSSByUpscaling()` realise modifications of the standard first-order uni- and multifractal analysis algorithms, respectively. They also implement modifications of the classical Hurst (a.k.a. R/S) analysis methods. For details, please refer to [1, 2]. The function `getScalingExponents()` implements the determination of a crossover and the related fractal dimensions, or scaling exponents at major and minor scales, in the data measure returned by the unifractal analysis algorithm.

### References
1. D.M. Filatov, A Method for Identification of Critical States of Open Stochastic Dynamical Systems Based on the Analysis of Acceleration, *J. Stat. Phys.*, 165 (2016) 681-692. <p><a href = "https://doi.org/10.1007/s10955-016-1641-6" rel = "nofollow"><img src = "https://zenodo.org/badge/DOI/10.1007/s10955-016-1641-6.svg" alt = "DOI:10.1007/s10955-016-1641-6" style = "vertical-align: top; max-width: 100%;"></a></p>
2. J.W. Kantelhardt, Fractal and Multifractal Time Series, available at https://arxiv.org/abs/0804.0747, 2008.