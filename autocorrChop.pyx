#  autocorr.pyx 
#  Created by nicain on 5/17/09.
#  Copyright (c) 2009 __MyCompanyName__. All rights reserved.

################################################################################
# Preamble: 
################################################################################

# Import floating point division from python 3.0
from __future__ import division

# Import necessary python packages:
import random, uuid
import numpy as np
from math import ceil
from math import floor
import time
cimport numpy as np
cimport cython


# Compile-time type initialization for numpy arrays:
ctypedef np.float_t DTYPE_t	

################################################################################
# Useful C++ functions:
################################################################################

# Wrapper for the RNG:
cdef extern from "MersenneTwister.h":
	ctypedef struct c_MTRand "MTRand":
		double randNorm( double mean, double stddev)
		void seed( unsigned long bigSeed[])

# External math wrapper functions that might be needed:
cdef extern from "math.h":
	float sqrt(float sqrtMe)
	float fabs(float absMe)		# Unused in this file, but can't hurt!
	float abs(float absMe)		# Unused in this file, but can't hurt!
	float exp(float expMe)
	
def trapezoidalRule(x,y):
	return (x[1]-x[0])*(y[0] + 2*np.sum(y[1:-2]) + y[-1])/2

################################################################################
# Create autocorrelation function:
################################################################################
@cython.boundscheck(False)
def autocorr(
			float mu,			# Steady-state mean
			float sigma,		# Steady-state std dev
			float tau,			# Time-constant
			float R,			# Chop
			float T,			# Max Time
			float dt,			# Time-step
			float maxS):
#			float maxY,			# max plotting S
#			int plotsOn):		# to plot, or not to plot...

	############################################################################
	# Initializations:
	############################################################################

	# Simple C initializations
	cdef int i										# Array indicies
	cdef float tau_recip = 1/tau					# Reciprocal of tau,
	cdef int maxTInd = ceil(T/dt)					# Compute max time index
	
	# RNG initializations:
	cdef unsigned long mySeed[624]		# Seed array for the RNG, length 624
	cdef c_MTRand myTwister				# RNG object construction
	
	# numpy array initializations:
	DTYPE = np.float					# Initialize a data-type for the array
	cdef np.ndarray[DTYPE_t, ndim=1] X = np.zeros(maxTInd, dtype=DTYPE) # 1D array of floats
	
	# Initialization of random number generator:
	myUUID = uuid.uuid4()
	random.seed(myUUID.int)
	for i in range(624): mySeed[i] = random.randint(0,int(exp(21)))
	myTwister.seed(mySeed)
	


	############################################################################
	# Simulation:
	############################################################################

	# Loop across times:
	tBegin = time.mktime(time.localtime())

	# Set up initial conditions:
	X[0] = mu + sigma*myTwister.randNorm(0,1)
	XOld = X[0]
	for i in range(1,maxTInd):
		X[i] = mu + exp(-dt*tau_recip)*(XOld - mu) + sigma*sqrt(1-exp(-2*dt*tau_recip))*myTwister.randNorm(0,1)
		XOld = X[i]
		if fabs(X[i]) < R*sigma:
			X[i] = 0

	# Create xcorr limits:
	XNew = X-X.mean()
	fx= np.fft.fft(XNew)
	fy = np.conj(fx)
	fxfy = fx*fy
	result = np.fft.fftshift(np.real(np.fft.ifft(fxfy)))/maxTInd
	v,ind = result.max(0),result.argmax(0)
	lInd, rInd = int(ind-maxS/dt),int(ind+maxS/dt)
	autoCorr = result[lInd:rInd]
	lags = (np.array(range(lInd,rInd))-ind)*dt

	# End time loop:
	tEnd = time.mktime(time.localtime())
	print 'Total Computation Time: ', time.strftime("H:%H M:%M S:%S",time.gmtime(tEnd - tBegin))

#	if plotsOn:
#		import pylab as pl
#		pl.plot(lags, autoCorr)
#		pl.xlim([lags.min(0),lags.max(0)])
#		pl.ylim([0,v])
#		pl.show()
	
	return autoCorr, lags
	
	
	
	
	
	
	
	
	