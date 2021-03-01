def error_propagation(func,variables,values,uncerts):
	r'''
	Function to compute propagated errors of an aribitrary function.
	sigma_f(x1,x2...) = sqrt(sigma_x1^2*(df/fx1)^2+sigma_x2^2*(df/dx2)^2+...)

	Parameters
	----------
	func : function
		function for which error propagation will be calculated
	vars : tuple
		tuple of strings definining the independent variables in func
	values : tuple
		tuple of arrays of values for each var
	uncerts : tuple
		tuple of arrays of uncertainties associated with each value
		
	Returns
	-------
	func_eval : array
		array of value of function evalulated at variables=values
	uncerts_prop : array
		array of propagated uncertaintities for func given values and unverts
	
	Notes
	-----
	Required modules: numpy, sympy
	Author: R. C. Levy (rlevy.astro@gmail.com)
	Last updated: 2021-03-01
	Change log:
		2021-03-01 : file created, RCL
		
	Examples
	--------
	>>> from ErrorPropagation import error_propagation
	>>> import numpy as np
	>>> def func(x,y):
	>>>	# z = 2xy^2-5
	>>> 	return 2.0*x*y**2-5.0
	>>> variables = ('x','y')
	>>> values = (np.linspace(1,10,10),np.linspace(-10,-1,10))
	>>> uncerts = (np.linspace(1,10,10)/10,np.linspace(-10,-1,10)/5)
	>>> z,ez = error_propagation(func,variables,values,uncerts)
	'''

	import numpy as np
	from sympy import (symbols, diff)

	#make the given variables symbolic
	vv = symbols(variables)

	#now take the derivative of func wtr each variable
	derivs = np.zeros((len(vv),len(values[0])))
	for i in range(len(vv)):
		d = diff(func(*vv),vv[i])
		for j in range(len(values[0])):
			#make a dictionary of the variable, value pairs
			subs = {}
			for k in range(len(variables)):
				subs[variables[k]] = values[k][j]
			#evaluate the derivative at that point
			derivs[i,j] = d.evalf(subs=subs)

	#evaluate the function at variables=values
	func_eval = func(*values)	

	#compute the final uncertainty
	uncerts_prop = np.sqrt(np.sum(np.asarray(uncerts)**2*derivs**2,axis=0))

	# return func_eval +/- uncerts_prop
	return func_eval,uncerts_prop






