from __future__ import print_function, division, absolute_import

import sys

import qminospy
from qminospy.me2 import ME_NLP

# python imports
from copy import copy
import re
from os.path import join
from collections import defaultdict
import pickle

# third party imports
import pandas
import cobra
from tqdm import tqdm
import numpy as np
import scipy

# COBRAme
import cobrame
from cobrame.util import building, mu, me_model_interface
from cobrame.io.json import save_json_me_model, save_reduced_json_me_model

# ECOLIme
import ecolime
from ecolime import (transcription, translation, flat_files, generics, formulas, compartments)
from ecolime.util.helper_functions import *


# Sensitivity functions

def add_dummy_demand(me,met_id,flux=0,fraction=0.01):
	'''
	This function adds a dummy demand reactions to calculate sensitivity and cost
	'''
	all_rxns_in_model = [rxn.id for rxn in me.reactions]
	met = me.metabolites.get_by_id(met_id)
	rxn_id = 'dummy_demand'
	weight = met.formula_weight/1000

	if rxn_id not in all_rxns_in_model:
	    rxn = cobrame.MEReaction(rxn_id)
	    me.add_reaction(rxn)
	    if isinstance(met,cobrame.Metabolite):
	        rxn.reaction = met_id + ' ->'
	        mass_correction = 0
	    elif isinstance(met,cobrame.Complex):
	        mass_correction = 'protein_biomass'
	        rxn.reaction = met_id + ' + ' + str(weight) + ' ' + mass_correction + ' ->'
	    else:
	        raise ValueError('Metabolite type not supported')
	else:
	    rxn = me.reactions.get_by_id(rxn_id)
	    
	if not flux and mass_correction:
	    # fraction 1% of total group production. This is done if production of a metabolite causes the
	    # production of biomass, e.g. protein_biomass, lipid_biomass, etc.
	    group_production = me.solution.x_dict[mass_correction+'_to_biomass']
	    flux = fraction * group_production / weight
	elif not flux:
	    flux = 0.001
	    
	rxn.lower_bound = flux
	rxn.upper_bound = flux
    
def single_flux_response(me,met_id,mu_fix=False,precision=1e-6):
	'''
	This function calculates flux response of a single metabolite. Response of growth and energy are
	sensitivity and cost, respectively.
	'''
	add_dummy_demand(me,met_id)
	solve_me_model(me, 0.2, min_mu = .05, using_soplex=False,\
		 precision = precision,verbosity=0,mu_fix=mu_fix)
	return met_id, me.solution.x_dict

def all_flux_responses(me,met_ids,mu_fix=False,solution=0,NP=1,precision=1e-6):
	'''
	This function calculates flux responses of several metabolites. Response of growth and energy are
	sensitivity and cost, respectively.
	'''
	# Correct number of threads if necessary
	NP = min([NP,len(met_ids)])
	if not solution:
	    solve_me_model(me, 0.2, min_mu = .05, using_soplex=False,\
	    	 precision = precision,mu_fix=mu_fix,verbosity=0)
	flux_dict = dict()
	flux_dict['base'] = me.solution.x_dict
	
	obj = [rxn.id for rxn in me.objective][0]

	if NP == 1:
	    for met_id in met_ids:
	        single_flux_response(me,met_id)
	        flux_dict[met_id] = me.solution.x_dict
	else:
	    import multiprocessing as mp
	    pool = mp.Pool(NP)
	    pbar = tqdm(total=len(met_ids))
	    pbar.set_description('{} response ({} threads)'.format(obj,NP))
	    def collect_result(result):
	    	pbar.update(1)
	        flux_dict[result[0]] = result[1]
	        
	    for met_id in met_ids:
	        pool.apply_async(single_flux_response, args=(me,met_id,mu_fix,\
	        		precision), callback=collect_result)
	    pool.close()
	    pool.join()   

	flux_results_df = pd.DataFrame.from_dict(flux_dict)
	return flux_results_df

def process_flux_responses(me,flux_results_df,base_rxn):
	base_flux = flux_results_df.loc[base_rxn]['base']
	met_ids = flux_results_df.columns.values
	processed_results_df = pd.DataFrame(index=met_ids,\
			columns=['molecular_weight','sensitivity'])

	for met_id in met_ids:
		if met_id == 'base':
			continue
		new_flux = flux_results_df.loc[base_rxn][met_id]
		met = me.metabolites.get_by_id(met_id)
		met_flux = flux_results_df.loc['dummy_demand'][met_id]
		response = (base_flux-new_flux)/met_flux
		MW = met.formula_weight/1000

		processed_results_df.loc[met_id]['molecular_weight'] = MW
		processed_results_df.loc[met_id][base_rxn+'_response'] = response

	return processed_results_df.drop('base')

def sensitivity(me,met_ids,NP=1,solution=0,biomass_dilution='biomass_dilution',\
		precision=1e-6):
	'''
	This function calculates sensitivity of all metabolites. 
	Parallel processing is activated when NP>1.
	'''

	flux_results_df = all_flux_responses(me,met_ids,mu_fix=False,\
			solution=solution,NP=NP,precision=1e-6)

	# Process results to get sensitivity
	sensitivity_df = process_flux_responses(me,flux_results_df,biomass_dilution)

	return sensitivity_df,flux_results_df

def biosynthetic_cost(me,met_ids,cost_rxn='ATPM',\
		biomass_dilution='biomass_dilution', NP=1,solution=0,precision=1e-6):
	'''
	This function calculates sensitivity of all metabolites. 
	Parallel processing is activated when NP>1.
	'''
	if not solution:
	    solve_me_model(me, max_mu = 0.2, min_mu = .05, \
	    		using_soplex=False, precision = 1e-6,verbosity=0)
	base_flux = me.solution.x_dict[cost_rxn]
	base_mu = me.solution.f

	# Growth rate-fixed solution.
	me.objective = cost_rxn
	mu_fix = 0.9*base_mu

	flux_results_df = all_flux_responses(me,met_ids,NP=NP,\
			mu_fix=mu_fix,solution=0,precision=1e-6)

	# Process results to get cost
	cost_df = process_flux_responses(me,flux_results_df,cost_rxn)

	return cost_df,flux_results_df