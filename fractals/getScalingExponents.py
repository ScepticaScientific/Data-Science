#!/home/ubuntu/miniconda/bin/python
#
# This function determines the optimal linear approximations of the data measure using two segments and returns
# the index of the corresponding boundary scale, the boundary scale itself (a.k.a. crossover), as well as the
# unifractal characteristics at the major and minor scales. For examples of using crossovers, see [1, 2].
#
# At the input, 'timeMeasure' is a time measure at different scales, while 'dataMeasure' is a data measure at the same
# scales.
#
# At the output, 'bScale' is the boundary scale, or crossover, separating the major and minor scales, 'bDM' is the
# data measure at the boundary scale, 'bsIndex' is the crossover's index with respect to the time measure, 'HMajor' is
# the unifractal dimension at the major scales, 'HMinor' is the unifractal dimension at the minor scales.
#
# REFERENCES:
# [1] D.M. Filatov, J. Stat. Phys., 165 (2016) 681-692. DOI: 10.1007/s10955-016-1641-6.
# [2] C.-K. Peng, S. Havlin, H.E. Stanley and A.L. Goldberger, Chaos, 5 (1995) 82–87. DOI: 10.1063/1.166141.
#
# The end user is granted perpetual permission to reproduce, adapt, and/or distribute this code, provided that
# an appropriate link is given to the original repository it was downloaded from.

import numpy as np

def getScalingExponents(timeMeasure, dataMeasure):
    ## Initialisation
    nScales = len(timeMeasure)

    log10tm = np.log10(timeMeasure)
    log10dm = np.log10(dataMeasure)

    res = 1.0e+07
    bsIndex = nScales

    ## Computing
    # We find linear approximations for major and minor subsets of the data measure and determine the index of the
    # boundary scale at which the approximations are optimal in the sense of best fitting to the data measure
    for i in range(3, nScales - 2 + 1):
        # Major 'i' scales are approximated by the function 'k * x + b' ...
        curr_log10tm = log10tm[nScales - i + 1 - 1 : nScales]
        curr_log10dm = log10dm[nScales - i + 1 - 1 : nScales]
        detA = i * np.sum(curr_log10tm ** 2.0) - np.sum(curr_log10tm) ** 2.0
        detK = i * np.sum(np.multiply(curr_log10tm, curr_log10dm)) - np.sum(curr_log10tm) * np.sum(curr_log10dm)
        detB = np.sum(curr_log10dm) * np.sum(curr_log10tm ** 2.0) - np.sum(curr_log10tm) * np.sum(np.multiply(curr_log10tm, curr_log10dm))
        k = detK / detA
        b = detB / detA
        # ... and the maximum residual is computed
        resMajor = max(np.abs(k * curr_log10tm + b - curr_log10dm))

        # Minor 'nScales - i + 1' scales are approximated by the function 'k * x + b' ...
        curr_log10tm = log10tm[1 - 1 : nScales - i + 1]
        curr_log10dm = log10dm[1 - 1 : nScales - i + 1]
        detA = (nScales - i + 1) * np.sum(curr_log10tm ** 2.0) - np.sum(curr_log10tm) ** 2.0
        detK = (nScales - i + 1) * np.sum(np.multiply(curr_log10tm, curr_log10dm)) - np.sum(curr_log10tm) * np.sum(curr_log10dm)
        detB = np.sum(curr_log10dm) * np.sum(curr_log10tm ** 2.0) - np.sum(curr_log10tm) * np.sum(np.multiply(curr_log10tm, curr_log10dm))
        k = detK / detA
        b = detB / detA
        # ... and the maximum residual is computed
        resMinor = max(np.abs(k * curr_log10tm + b - curr_log10dm))

        if (resMajor ** 2.0 + resMinor ** 2.0 < res):
            res = resMajor ** 2.0 + resMinor ** 2.0
            bsIndex = i

    # Now we determine the boundary scale and the boundary scale's data measure, ...
    bScale = 2.0 * timeMeasure[1 - 1] / timeMeasure[nScales - bsIndex + 1 - 1] / 2.0
    bDM = dataMeasure[nScales - bsIndex + 1 - 1]
    # ... as well as compute the unifractal dimensions using the boundary scale's index:
    # at the major 'bsIndex' scales ...
    curr_log10tm = log10tm[nScales - bsIndex + 1 - 1 : nScales]
    curr_log10dm = log10dm[nScales - bsIndex + 1 - 1 : nScales]
    detA = bsIndex * np.sum(curr_log10tm ** 2.0) - np.sum(curr_log10tm) ** 2.0
    detK = bsIndex * np.sum(np.multiply(curr_log10tm, curr_log10dm)) - np.sum(curr_log10tm) * np.sum(curr_log10dm)
    DMajor = detK / detA
    HMajor = -DMajor
    # ... and at the minor 'nScales - bsIndex + 1' scales
    curr_log10tm = log10tm[1 - 1 : nScales - bsIndex + 1]
    curr_log10dm = log10dm[1 - 1 : nScales - bsIndex + 1]
    detA = (nScales - bsIndex + 1) * np.sum(curr_log10tm ** 2.0) - np.sum(curr_log10tm) ** 2.0
    detK = (nScales - bsIndex + 1) * np.sum(np.multiply(curr_log10tm, curr_log10dm)) - np.sum(curr_log10tm) * np.sum(curr_log10dm)
    DMinor = detK / detA
    HMinor = -DMinor

    return [bScale, bDM, bsIndex, HMajor, HMinor]
