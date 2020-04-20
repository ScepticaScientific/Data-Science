# Data-Science/fractals
This folder contains utilities implementing various fractal analysis methods. To get familiar, run `test.py`.

### Fractal analysis algorithms
Both the uni- (`getHurstByUpscaling()`) and multifractal (`getMSSByUpscaling()`) analysis algorithms realise modifications of the standard first-order detrended fluctuation analysis (DFA/MF-DFA) methods. They also implement the classical Hurst (a.k.a. R/S) analysis algorithms. For details, please refer to [1, 2].

### References
1. D.M. Filatov, A Method for Identification of Critical States of Open Stochastic Dynamical Systems Based on the Analysis of Acceleration, *J. Stat. Phys.*, 165 (2016) 681-692.
2. J.W. Kantelhardt, Fractal and Multifractal Time Series, available at https://arxiv.org/abs/0804.0747, 2008.