#
#  run_OUFit.py 
#  OUFit
#
#  Created by nicain on 5/18/10.
#  Copyright (c) 2010 __MyCompanyName__. All rights reserved.
#

################################################################################
# Preamble:
################################################################################

# Import necessary packages:
from subprocess import call as call

# Compile cython extension cythonOU.pyx:
call('python setup.py build_ext --inplace', shell=True)

# Import cython extension:
from autocorrChop import autocorr
import pbsTools as pt

################################################################################
# Call the main function:
################################################################################

# Settings:
FFTN = 2**18	#23 max
dt = .01		#.001 min
saveFileName = 'autocorr.dat'
saveResults = False
displayResults = True

# Parameters:
mu = .1
sigma = 2
tau = .1
R = 0
maxS = .3

# Initializations:
T = FFTN*dt
		
# Do the computation:
A, lags =  autocorr(mu, sigma, tau, R, T, dt, maxS)

################################################################################
# Generating Output:
################################################################################

# Save Results:
if saveResults is True:
	pt.pickle({
		'mu':mu,
		'sigma':sigma,
		'tau':tau,
		'R':R,
		'A':h0,
		'lags':Ez,
		'dt':dt,
		'FFTN':FFTN
		}, saveFileName = saveFileName)

# Display results:
if displayResults is True:
	print A, lags





