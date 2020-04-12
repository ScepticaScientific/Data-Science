#!/home/ubuntu/miniconda/bin/python
#
# This file is the entry point for testing the modified first-order DFA ('getHurstByUpscaling.py') and
# the modified first-order multifractal DFA ('getMSSByUpscaling.py'). The initial dataset is a binomial cascade [1].
#
# References:
# [1] J. Feder, Fractals, Plenum Press, New York, 1988.

from getHurstByUpscaling import getHurstByUpscaling
from getMSSByUpscaling import getMSSByUpscaling
import numpy as np
from numpy import loadtxt
import matplotlib.pyplot as plt

## Loading data
dx = loadtxt('binomial.txt')
dx = dx[1 - 1 : 8192]               # We take the first 8192 samples

## Computing
# Modified first-order DFA
[timeMeasure, meanDataMeasure, scales] = getHurstByUpscaling(dx)                    # Set of parameters No. 1
#[timeMeasure, meanDataMeasure, scales] = getHurstByUpscaling(dx, 3.0, 0, 2.0)       # Set of parameters No. 2

# Modified first-order MF-DFA
[_, dataMeasure, _, out, q] = getMSSByUpscaling(dx, isNormalised = 1)

## Output
# Modified first-order DFA
plt.figure()
plt.subplot(2, 1, 1)
plt.loglog(timeMeasure, meanDataMeasure, 'ko-')
plt.xlabel(r'$\mu(t)$')
plt.ylabel(r'$\mu(\Delta x)$')
plt.grid('on', which = 'minor')
plt.title('Modified First-Order DFA of a Binomial Cascade')
plt.subplot(2, 1, 2)
plt.loglog(scales, meanDataMeasure, 'ko-')
plt.xlabel(r'$j$')
plt.ylabel(r'$\mu(\Delta x)$')
plt.grid('on', which = 'minor')

# Modified first-order MF-DFA
print('alpha_min = %g, alpha_max = %g, dalpha = %g' % (out['LH_min'], out['LH_max'], out['LH_max'] - out['LH_min']))
print('h_min = %g, h_max = %g, dh = %g\n' % (out['h_min'], out['h_max'], out['h_max'] - out['h_min']))

plt.figure()
plt.subplot(2, 1, 1)
plt.plot(out['LH'], out['f'], 'ko-')
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$f(\alpha)$')
plt.grid('on', which = 'major')
plt.title('Statistics of Modified First-Order MF-DFA of a Binomial Cascade')

plt.subplot(2, 1, 2)
plt.plot(q, out['tau'], 'ko-')
plt.xlabel(r'$q$')
plt.ylabel(r'$\tau(q)$')
plt.grid('on', which = 'major')

plt.figure()
nq = np.int(len(q))
leg_txt = []
for qi in range(1, nq + 1):
    llh = plt.loglog(scales, dataMeasure[qi - 1, :], 'o-')
    leg_txt.append('tau = %g (q = %g)' % (out['tau'][qi - 1], q[qi - 1]))
plt.xlabel(r'$j$')
plt.ylabel(r'$\mu(\Delta x, q)$')
plt.grid('on', which = 'minor')
plt.title('Modified First-Order MF-DFA of a Binomial Cascade')
plt.legend(leg_txt)

plt.show()
