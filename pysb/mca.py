from collections import namedtuple

from pysb.integrate import odesolve


sens_res = namedtuple('SensitivityResult', ['param', 'left', 'right', 'delta'])

def time_dep_sensitivity(model, t, results = None, delta=0.1, parameters=[]):
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

		all_param_results[param.name] = sens_res(param.name,
												 new_results_left,
												 new_results_right,
												 delta * param.value)
	
	return all_param_results
