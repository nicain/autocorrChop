def autocorr(
			mu,			# Steady-state mean
			sigma,		# Steady-state std dev
			tau,			# Time-constant
			R,			# Chop
			T,			# Max Time
			dt,			# Time-step
			maxS,
			maxY,			# max plotting S
			plotsOn):		# to plot, or not to plot...

	############################################################################
	# Initializations:
	############################################################################

	# Packages:
	import pylab as pl
	import numpy as np
	import random
	import time

	# Stuff I need:
	tau_recip = 1/tau					# Reciprocal of tau,
	maxTInd = int(np.ceil(T/dt))				# Compute max time index
	X = np.zeros(maxTInd, type=float32) # 1D array of floats
	
	############################################################################
	# Simulation:
	############################################################################

	# Loop across times:
	tBegin = time.mktime(time.localtime())

	# Set up initial conditions:
	X[0] = mu + sigma*random.normalvariate(0,1)
	XOld = X[0]
	for i in range(1,maxTInd):
		X[i] = mu + np.exp(-dt*tau_recip)*(XOld - mu) + sigma*np.sqrt(1-np.exp(-2*dt*tau_recip))*random.normalvariate(0,1)
		XOld = X[i]
		if np.abs(X[i]) < R*sigma:
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

	if plotsOn:
		pl.plot(lags, autoCorr)
		pl.xlim([lags.min(0),lags.max(0)])
		pl.ylim([0,maxY])
		pl.show()
	
	return autoCorr, lags