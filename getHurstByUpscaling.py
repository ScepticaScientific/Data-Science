#!/home/ubuntu/miniconda/bin/python
#
# This code implements a modification of the first-order detrended fluctuation analysis (DFA) method
# originally described in [1]. For more details on the DFA and Hurst analysis techniques, please refer to [2, 3].
#
# INPUT:
# dx              - a time series of increments of the quantity 'x(t)', of the length equal to
#                   an integer power of two greater than two (i.e. 4, 8, 16, 32, etc.)
# normType_p      - any real greater than or equal to one, specifying the p-norm
# isDFA           - a boolean value prescribing to use either the DFA-based algorithm (1) or
#                   the standard Hurst (R/S) analysis (0)
# normType_q      - any real greater than or equal to one, specifying the q-norm
#
# OUTPUT:
# timeMeasure     - the time measure of the data's support at different scales
# meanDataMeasure - the data measure at different scales
# scales          - the scales at which the data measure is computed
#
# The conventional way of using the output values is to plot the data measure vs the scales; the time measure,
# being the inverse quantity to the scales, is computed for an alternative representation and may be ignored.
#
# REFERENCES:
# [1] D.M. Filatov, J. Stat. Phys., 165 (2016) 681-692. DOI: 10.1007/s10955-016-1641-6
# [2] J.W. Kantelhardt, Fractal and Multifractal Time Series, available at http://arxiv.org/abs/0804.0747, 2008
# [3] J. Feder, Fractals, Plenum Press, New York, 1988
#
# The end user is granted perpetual permission to reproduce, adapt, and/or distribute this code, provided that
# an appropriate link is given to the original code at the repository it was downloaded from.

import numpy as np

def getHurstByUpscaling(dx, normType_p = np.inf, isDFA = 1, normType_q = 1.0):
    ## Some initialiation
    dx_len = np.int(len(dx))

    # We have to reserve the most major scale for shifts, so we divide the data
    # length by two. (As a result, the time measure starts from 2.0, not from
    # 1.0, see below.)
    dx_len = np.int(dx_len / 2)

    dx_shift = np.int(dx_len / 2)

    nScales = np.int(np.round(np.log(dx_len * 1.0) / np.log(2.0)))    # Number of scales involved. P.S. We use 'round()' to prevent possible malcomputing of the logarithms
    j = 2 ** (np.arange(1, nScales + 1) - 1) - 1

    meanDataMeasure = np.zeros(nScales)

    ## Computing the data measure
    for ji in range(1, nScales + 1):
        # At the scale 'j(ji)' we deal with '2 * (j(ji) + 1)' elements of the data 'dx'
        dx_k_len = 2 * (j[ji - 1] + 1)
        n = np.int(dx_len / dx_k_len)

        dx_leftShift = np.int(dx_k_len / 2)
        dx_rightShift = np.int(dx_k_len / 2)

        for k in range(1, n + 1):
            # We get a portion of the data of the length '2 * (j(ji) + 1)' plus the data from the left and right boundaries
            dx_k_withShifts = dx[(k - 1) * dx_k_len + 1 + dx_shift - dx_leftShift - 1 : k * dx_k_len + dx_shift + dx_rightShift]

            # Then we perform free upscaling and, using the above-selected data (provided at the scale j = 0),
            # compute the velocities at the scale 'j(ji)'
            j_dx = np.convolve(dx_k_withShifts, np.ones(dx_rightShift), 'valid')

            # Then we compute the accelerations at the scale 'j(ji) + 1'
            r = (j_dx[1 + dx_rightShift - 1 : ] - j_dx[1 - 1 : -dx_rightShift]) / 2.0

            # Finally, we compute the range ...
            if (normType_p == 0):
                R = np.max(r[2 - 1 : ]) - np.min(r[2 - 1 : ])
            elif (np.isinf(normType_p)):
                R = np.max(np.abs(r[2 - 1 : ]))
            else:
                R = (np.sum(r[2 - 1 : ] ** normType_p) / len(r[2 - 1 : ])) ** (1.0 / normType_p)
            # ... and the normalisation factor ("standard deviation")
            S = np.sqrt(np.sum(np.abs(np.diff(r)) ** 2.0) / (len(r) - 1))
            if (isDFA == 1):
                S = 1.0

            meanDataMeasure[ji - 1] += (R / S) ** normType_q
        meanDataMeasure[ji - 1] = (meanDataMeasure[ji - 1] / n) ** (1.0 / normType_q)

    # We pass from the scales ('j') to the time measure; the time measure at the scale j(nScales) (the most major one)
    # is assumed to be 2.0, while it is growing when the scale is tending to j(1) (the most minor one).
    # (The scale j(nScales)'s time measure is NOT equal to 1.0, because we reserved the highest scale for shifts
    # in the very beginning of the function.)
    timeMeasure = 2.0 * dx_len / (2 * (j + 1))

    scales = j + 1

    return [timeMeasure, meanDataMeasure, scales]
