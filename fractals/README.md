# Data-Science/fractals
This folder contains utilities implementing various fractal analysis methods. To get familiar, run `test.py`.

### Fractal analysis algorithms
Both the uni- (`getHurstByUpscaling()`) and multifractal (`getMSSByUpscaling()`) analysis algorithms realise modifications of the standard first-order detrended fluctuation analysis (DFA/MF-DFA) methods. They also implement the classical Hurst (a.k.a. R/S) analysis algorithms. For details, please refer to [1, 2].

### References
1. D.M. Filatov, A Method for Identification of Critical States of Open Stochastic Dynamical Systems Based on the Analysis of Acceleration, *J. Stat. Phys.*, 165 (2016) 681-692, <a href="https://doi.org/10.1007/s10955-016-1641-6" rel="nofollow"><img src="https://camo.githubusercontent.com/c40d7831d685ef53656a8296050bcf378e3fadf8/68747470733a2f2f7a656e6f646f2e6f72672f62616467652f444f492f31302e313030372f3937382d332d3331392d37363230372d345f31352e737667" alt="DOI:10.1007/s10955-016-1641-6" data-canonical-src="https://zenodo.org/badge/DOI/10.1007/s10955-016-1641-6.svg" style="max-width:100%;"></a>
2. J.W. Kantelhardt, Fractal and Multifractal Time Series, available at https://arxiv.org/abs/0804.0747, 2008.