from collections import namedtuple

from pysb.integrate import odesolve


SensitivityResult = namedtuple('SensitivityResult', ['param', 'left', 'right', 'delta'])

def time_dep_sensitivity(model, t, results = None, delta=0.1, parameters=[], normalize=True): 
	if not parameters:
		return
	all_param_results = {}

	if results is None:
		results = odesolve(model, t)
	print 'Solving time-dependent sensitivity analysis for %d parameters' % len(parameters)

	for param in parameters:
		if param.name.startswith('__'):
			continue
		print 'Solving for', param.name
		
		oldval = param.value

		newval_left = param.value * (1. + 0.5 * delta)
		model.parameters[param.name].value = newval_left
		new_results_left = odesolve(model, t)

		newval_right = param.value * (1. - 0.5 * delta)
		model.parameters[param.name].value = newval_right
		new_results_right = odesolve(model, t)

		model.parameters[param.name].value = oldval

		all_param_results[param.name] = SensitivityResult(param.name,
												 new_results_left,
												 new_results_right,
												 delta * param.value)

	sens_coeffs = {}
	for pname, sens_result in all_param_results.iteritems():
		coeff = {}
		for obs in model.observables:
			coeff[obs.name] = (sens_result.right[obs.name] - sens_result.left[obs.name]) / sens_result.delta
			if normalize:
				coeff[obs.name] *= model.parameters[pname].value / resuts[obs.name]
		sens_coeffs[pname] = coeff

	
	return sens_coeffs

