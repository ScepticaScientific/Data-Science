#!/home/ubuntu/miniconda/bin/python
#
# This code implements a modification of the first-order unifractal analysis algorithm originally described in [1].
# It covers both the detrended fluctuation analysis (DFA) and the Hurst (a.k.a. R/S) analysis methods. For more details
# on the DFA and Hurst analysis methods, please refer to [2, 3].
#
# At the input, 'dx' is a time series of increments of the physical observable 'x(t)', of an arbitrary length greater 
# than 3, 'normType_p' is any real greater than or equal to one specifying the p-norm, 'isDFA' is a boolean value 
# prescribing to use either the DFA-based algorithm or the standard Hurst (a.k.a. R/S) analysis, 'normType_q' is 
# any real greater than or equal to one specifying the q-norm.
#
# At the output, 'timeMeasure' is the time measure of the data's support at different scales, 'meanDataMeasure' is
# the data measure at different scales, while 'scales' is the scales at which the data measure is computed.
#
# The conventional way of using the output values is to plot the data measure vs the scales; the time measure,
# being the inverse quantity to the scales, is computed for an alternative representation and may be ignored.
#
# REFERENCES:
# [1] D.M. Filatov, J. Stat. Phys., 165 (2016) 681-692. DOI: 10.1007/s10955-016-1641-6.
# [2] J.W. Kantelhardt, Fractal and Multifractal Time Series, available at http://arxiv.org/abs/0804.0747, 2008.
# [3] J. Feder, Fractals, Plenum Press, New York, 1988.
#
# The end user is granted perpetual permission to reproduce, adapt, and/or distribute this code, provided that
# an appropriate link is given to the original repository it was downloaded from.

import numpy as np

def getHurstByUpscaling(dx, nMajorScalesToDo = np.inf, normType_p = np.inf, isDFA = 1, normType_q = 1.0):
    ## Some initialiation
    dx_backward = np.flipud(dx)
    
    dx_len = len(dx)

    # We have to reserve the most major scale for shifts, so we divide the data
    # length by two. (As a result, the time measure starts from 2.0, not from
    # 1.0, see below.)
    dx_len = int(dx_len / 2)

    dx_shift = int(dx_len / 2)

    j = np.arange(0, int(dx_len / 2) - 1 + 1)
    nScales = len(j)
    
    meanDataMeasure = np.zeros(nScales)

    ## Computing the data measure
    for ji in range(1, nScales + 1):
        # At the scale 'j(ji)' we deal with '2 * (j(ji) + 1)' elements of the data 'dx'
        dx_k_len = 2 * (j[ji - 1] + 1)
        n = int(dx_len / dx_k_len)

        dx_leftShift = int(dx_k_len / 2)
        dx_rightShift = int(dx_k_len / 2)
        
        for ip in range(1, 2 + 1):
            if (ip == 1):
                dx_aux = dx                 # Forward processing
            else:
                dx_aux = dx_backward        # Backward processing

            for k in range(1, n + 1):
                # We get a portion of the data of the length '2 * (j(ji) + 1)' plus the data from the left and right boundaries
                dx_k_withShifts = dx_aux[(k - 1) * dx_k_len + 1 + dx_shift - dx_leftShift - 1 : k * dx_k_len + dx_shift + dx_rightShift]
                if (ip == 2):
                    dx_k_withShifts = np.flipud(dx_k_withShifts)

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
        meanDataMeasure[ji - 1] = (meanDataMeasure[ji - 1] / n / 2.0) ** (1.0 / normType_q)

    # We pass from the scales ('j') to the time measure; the time measure at the scale j(nScales) (the most major one)
    # is assumed to be 2.0, while it is growing when the scale is tending to j(1) (the most minor one).
    # (The scale j(nScales)'s time measure is NOT equal to 1.0, because we reserved the highest scale for shifts
    # in the very beginning of the function.)
    timeMeasure = 2.0 * dx_len / (2 * (j + 1))

    scales = j + 1
    
    # This code is the determination of the Hurst parameter 'H' over all the scales (the classical approach)
    log10tm = np.log10(timeMeasure)
    log10dm = np.log10(meanDataMeasure)

    D = np.sum(np.multiply(log10tm, (log10dm - np.mean(log10dm)))) / np.sum(np.multiply(log10tm, (log10tm - np.mean(log10tm))))
    H = -D

    # This code is the determination of the Hurst parameter 'H' over the major 'nMajorScalesToDo' scales ...
    nMajorScalesDone = np.min(np.array([nScales, nMajorScalesToDo]))
    log10tm = np.log10(timeMeasure[-nMajorScalesDone + 1 - 1 : ])
    log10dm = np.log10(meanDataMeasure[-nMajorScalesDone + 1 - 1 : ])

    DMajor = np.sum(np.multiply(log10tm, (log10dm - np.mean(log10dm)))) / np.sum(np.multiply(log10tm, (log10tm - np.mean(log10tm))))
    HMajor = -DMajor

    # ... and the Hurst parameter 'H' over the rest minor scales
    log10tm = np.log10(timeMeasure[1 - 1 : -nMajorScalesDone + 1])
    log10dm = np.log10(meanDataMeasure[1 - 1 : -nMajorScalesDone + 1])

    DMinor = np.sum(np.multiply(log10tm, (log10dm - np.mean(log10dm)))) / np.sum(np.multiply(log10tm, (log10tm - np.mean(log10tm))))
    HMinor = -DMinor

    return [timeMeasure, meanDataMeasure, scales, H, HMajor, HMinor, nMajorScalesDone]
    
