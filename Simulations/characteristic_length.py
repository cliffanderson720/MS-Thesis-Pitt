def decel_rate(RPM,
				RPM_max,
				thresholds = (2000,700,420,350,300,0),
				decel_rates = (200,50,21.48,12.1209,8.398,1.953)):
				
	'''Gives the expected deceleration rate at a given RPM'''
	if RPM>RPM_max:
		print('Input RPM is greater than RPM_max')
		return
		
	# get indices of thresholds where RPM_max is below any of the threshold values
	indices = [i for i,x in enumerate(thresholds) if RPM_max>x]
	
	# subset thresholds and decel_rates
	thresholds = ((RPM_max,)+thresholds[min(indices):max(indices)+1])
	decel_rates = decel_rates[min(indices):max(indices)+1]

	# get position between thresholds
	nthresholds=len(thresholds)
	
	for i in range(nthresholds-1):

		is_below_left  = RPM < thresholds[i]
		is_above_right = RPM > thresholds[i+1]
		
		# if RPM is between threshold values, return the rate at that index
		if (is_below_left and is_above_right):
			return decel_rates[i]
	return 

def profile(t,t_hold,RPM_max,accel_rate=500,**kwargs):
	'''generates position in a spin profile'''
	return
	
	
	
	
	
	