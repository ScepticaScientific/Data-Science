# Data-Science/fractals
This folder contains utilities implementing various fractal analysis methods. To get familiar, run `test.py`.

### Fractal analysis algorithms
Both the uni- (`getHurstByUpscaling()`) and multifractal (`getMSSByUpscaling()`) analysis algorithms realise modifications of the standard first-order detrended fluctuation analysis (DFA/MF-DFA) methods. They also implement the classical Hurst (a.k.a. R/S) analysis algorithms. For details, please refer to [1, 2].

### References
1. D.M. Filatov, A Method for Identification of Critical States of Open Stochastic Dynamical Systems Based on the Analysis of Acceleration, *J. Stat. Phys.*, 165 (2016) 681-692. <p><a href="https://doi.org/10.1007/s10955-016-1641-6" rel = "nofollow"><img src = "https://zenodo.org/badge/DOI/10.1007/s10955-016-1641-6.svg" alt = "DOI:10.1007/s10955-016-1641-6" style = "vertical-align: top; max-width: 100%;"></a></p>
2. J.W. Kantelhardt, Fractal and Multifractal Time Series, available at https://arxiv.org/abs/0804.0747, 2008.