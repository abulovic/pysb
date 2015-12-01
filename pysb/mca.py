from random import random, shuffle
from collections import namedtuple, defaultdict

from sympy import diff
from sympy.core.numbers import Integer
from sympy.core.symbol import Symbol
from sympy.core.add import Add

import numpy as np
from pysb.integrate import odesolve

def get_param_bounds(value, delta):
	return value * (1. - 0.5 * delta), value * (1. + 0.5 * delta)


def estimate_point_sensitivity(model, t, results, pname, pvalue, delta):
	param = model.parameters[pname]
	p_oldval = param.value

	# move slightly to the left and to the right
	# and thus derive an approximation of the derivative
	# at that point
	newval_left, newval_right = get_param_bounds(param.value, delta)

	model.parameters[pname].value = newval_left
	nr_left = odesolve(model, t)

	model.parameters[pname].value = newval_right
	nr_right = odesolve(model, t)

	model.parameters[pname].value = p_oldval

	result = PointSens(model, pname, value, delta, results, nr_left, nr_right)


class SensData(object):
	def __init__(self, dtc, scaled_dtc, compute_stats=True):
		self.dtc = dtc
		self.scaled_dtc = scaled_dtc
		self.scaled_dtc[np.isnan(self.scaled_dtc)] = 0.
		if compute_stats:
			self.mean = np.mean(self.dtc)
			self.scaled_mean = np.mean(self.scaled_dtc)
			self.std = np.std(self.dtc)
			self.scaled_std = np.std(self.scaled_dtc)


class PointSens(object):
	def __init__(self, model, pname, value, delta, results, results_left, results_right):
		self.pname = pname
		self.pvalue = value
		self.sensitivities = {}
		for obs in model.observables:
			cname = obs.name
			obs_der = (results_right[obs.name] - results_left[obs.name]) / delta
			scaled_obs_der = obs_der / value * results[obs.name]
			self.sensitivities[obs.name] = SensData(obs_der, scaled_obs_der, True)



SensitivityResult = namedtuple('SensitivityResult', ['param', 'left', 'right', 'delta'])

def point_derivative(model, t, results, pname, value, delta):
	'''
	Estimate local sensitivity for a particular parameter
	around a particular point value.
	'''
	param = model.parameters[pname]
	p_oldval = param.value

	# move slightly to the left and to the right
	# and thus derive an approximation of the derivative
	# at that point
	newval_left = param.value * (1. - 0.5 * delta)
	newval_right = param.value * (1. + 0.5 * delta)

	model.parameters[pname].value = newval_left
	nr_left = odesolve(model, t)

	model.parameters[pname].value = newval_right
	nr_right = odesolve(model, t)

	model.parameters[pname].value = p_oldval

	ps = PointSens(model, pname, p_oldval, delta, results, nr_left, nr_right)
	return ps


def LSA(model, t, results=None, delta=0.001, parameters=[]):
	if not parameters:
		return
	if results is None:
		results = odesolve(model, t)

	param_sens = {}
	for param in filter(lambda p: not p.name.startswith('__'), parameters):
		presult = point_derivative(model, t, results, param.name, param.value, delta)
		param_sens[param.name] = presult

	return param_sens



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
	param_lhs = latin_hypercube_sampling(params, 2)

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

def analyze_eq(eq):
	'''Lists species and rate constants mentioned in the equation'''
	eq_atoms = list(eq.atoms())
	rate_constants = []
	species = []
	for atom in eq_atoms:
		if isinstance(atom, Integer):
			continue
		elif isinstance(atom, Symbol):
			try:
				float(atom.name)
			except ValueError, e:
				if atom.name.startswith('__'):
					species.append(atom)
				else:
					rate_constants.append(atom)
	return species, rate_constants

def force_eval(eq):
	atoms = list(filter(lambda a: isinstance(a, Symbol), eq.atoms()))
	subs_vals = [(a, float(a.name)) for a in atoms]
	return eq.subs(subs_vals)


def symLSA(model, t, results):
	ALL = sum(model.odes)
	species, rate_constants = analyze_eq(ALL)
	all_species = {s.name: s for s in species}
	all_rate_constants = {rc.name: rc for rc in rate_constants}

	spec_values = [()]
	spec_values = [(sp, results[sp.name][-1]) for sp in all_species.values()]
	rc_values = [(rc, model.parameters[rc.name].value) for rc in all_rate_constants.values()]

	obs_eqs = {}
	derivatives = {}
	der_values = {}
	scaled_values = {}
	total = 0
	errors = 0
	for obs in model.observables:
		master_eq = sum([model.odes[idx] for idx in obs.species])
		if 'pp' in obs.name:
			print master_eq
		if not master_eq:
			continue
		species, rate_constants = analyze_eq(master_eq)
		derivatives[obs.name] = {}
		der_values[obs.name] = {}
		scaled_values[obs.name] = {}
		for rc in rate_constants:
			D = diff(master_eq, rc)
			_species, _rates = analyze_eq(D)
			spec_values = [(sp, results[sp.name][-1]) for sp in _species]

			derivatives[obs.name][rc.name] = D
			subs_eq = derivatives[obs.name][rc.name].subs(spec_values).subs(rc_values)
			if isinstance(subs_eq, Add):
				try:
					subs_eq = force_eval(subs_eq)
				except Exception:
					import pdb
					pdb.set_trace()
			der_values[obs.name][rc.name] = subs_eq
			scaled_values[obs.name][rc.name] = subs_eq * model.parameters[rc.name].value / results[obs.name][-1]
			print obs.name, rc.name
			print D
			print der_values[obs.name][rc.name]
			print scaled_values[obs.name][rc.name]
			print
	print errors, '/', total
				
		
