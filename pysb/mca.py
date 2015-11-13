from random import random, shuffle
from collections import namedtuple, defaultdict

import numpy as np
from pysb.integrate import odesolve


SensitivityResult = namedtuple('SensitivityResult', ['param', 'left', 'right', 'delta'])

def time_dep_sensitivity(model, t, results = None, delta=0.001, parameters=[], normalize=True): 
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
			coeff[obs.name] = (sens_result.left[obs.name] - sens_result.right[obs.name]) / sens_result.delta
			if normalize:
				coeff[obs.name] *= model.parameters[pname].value / results[obs.name]
				if results[obs.name][0] == 0.:
					coeff[obs.name][0] = 0.
		sens_coeffs[pname] = coeff

	
	return sens_coeffs


def average(model, sensitivities):
	avg_sens = {}
	for pname, coeff in sensitivities.iteritems():
		averages = {}
		for obs in model.observables:
			averages[obs.name] = np.average(coeff[obs.name])
		avg_sens[pname] = averages
	return avg_sens

def latin_hypercube_sampling(parameters, M):
	# parameters: list of tuples (name, values)
	random_values = []
	for i in range(M):
		random_values.append([])
	for pname, value in parameters:
		l = 0
		r = 2 * value
		_len = r - l
		intervals = [(l + _len / M * i, l + _len / M * (i + 1)) for i in range(M)]
		print pname
		print len(intervals)
		print intervals
		vals = [li + ri * random() for li, ri in intervals]
		order = range(M)
		shuffle(order)
		for idx, val in enumerate(order):
			random_values[idx].append(vals[val])
	return np.array(random_values)

def global_sensitivity(model, t, results = None, delta=0.001, parameters=[], normalize=True):
	if not parameters:
		return
	params = [(p.name, p.value) for p in parameters]
	param_lhs = latin_hypercube_sampling(params, 10)

	if results is None:
		results = odesolve(model, t)

	sens_coeffs = {}
	for idx, (pname, pvalue) in enumerate(params):
		sens_coeffs[pname] = {}

		lhs_values = param_lhs[:,idx]
		print 'Solving for', parameters[idx]

		for param_val in lhs_values:
			coeff = {}
			model.parameters[pname].value = param_val
			perturbed_result = odesolve(model, t)
			res = (SensitivityResult(pname,
									results,
									perturbed_result,
									param_val - pvalue))
			for obs in model.observables:
				coeff[obs.name] = (res.right[obs.name] / res.left[obs.name])
			sens_coeffs[pname][param_val] = coeff




		model.parameters[pname].value = pvalue
	return sens_coeffs