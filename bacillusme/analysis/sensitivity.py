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
import bacillusme
from bacillusme import (transcription, translation, flat_files, generics, formulas, compartments)
from bacillusme.util.helper_functions import *


# Sensitivity functions

def add_dummy_demand(me,met_id,flux=0,fraction=0.01):
    '''
    This function adds a dummy demand reactions to calculate sensitivity and cost
    '''
    met = me.metabolites.get_by_id(met_id)
    rxn_id = 'dummy_demand'
    weight = met.formula_weight/1000

    # Replace dummy_demand if existent
    if hasattr(me.reactions,rxn_id):
        me.reactions.get_by_id(rxn_id).remove_from_model()
        
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
    
    if not flux and mass_correction:
        # fraction 1% of total group production. This is done if production of a metabolite causes the
        # production of biomass, e.g. protein_biomass, lipid_biomass, etc.
        group_production = me.solution.x_dict[mass_correction+'_to_biomass']
        flux = fraction * group_production / weight
    elif not flux:
        flux = 0.001
        
    rxn.lower_bound = flux
    rxn.upper_bound = flux
    
def single_flux_response(me,met_id,mu_fix=False,precision=1e-6,single_change_function=False,growth_key='mu'):
    '''
    This function calculates flux response of a single metabolite. Response of growth and energy are
    sensitixvity and cost, respectively.
    '''
    
    if isinstance(met_id,list):
        for met in met_id:
            if single_change_function=='transporter': # close_transporter function
                close_transporter(me,met)
            elif single_change_function == 'overexpress':
                overexpress_transporter(me,met)
            elif single_change_function == 'group_knockout':
                group_knockout(me,met)
            else: # Just normal sensitivity/cost calculation
                add_dummy_demand(me,met)
        met_id = met
    else:
        if single_change_function=='transporter': # close_transporter function
            close_transporter(me,met_id)
        elif single_change_function == 'overexpress':
            overexpress_transporter(me,met_id)
        elif single_change_function == 'group_knockout':
            group_knockout(me,met_id)
        else: # Just normal sensitivity/cost calculation
            add_dummy_demand(me,met_id)

    solve_me_model(me, 0.2, min_mu = .05, using_soplex=False,\
        precision = precision,verbosity=0,mu_fix=mu_fix,growth_key=growth_key)
    
    if me.solution:
        flux_dict = me.solution.x_dict
    else:
        flux_dict = {r.id:0. for r in me.reactions}
    return met_id, flux_dict

def all_flux_responses(me,met_ids,mu_fix=False,solution=0,NP=1,precision=1e-6,
                       single_change_function=False,growth_key = 'mu',sequential=False):
    '''
    This function calculates flux responses of several metabolites. Response of growth and energy are
    sensitivity and cost, respectively.
    '''
    import copy
    # Correct number of threads if necessary
    NP = min([NP,len(met_ids)])
    if not solution:
        solve_me_model(me, 0.2, min_mu = .05, using_soplex=False,precision = precision,mu_fix=mu_fix,verbosity=0)
    flux_dict = dict()
    flux_dict['base'] = me.solution.x_dict

    obj = [rxn.id for rxn in me.objective][0]
    
    
    if NP == 1:
        for met_id in tqdm(met_ids):
            og_me = copy.deepcopy(me)
            single_flux_response(og_me,met_id,single_change_function=single_change_function,\
                                 mu_fix=mu_fix,precision=precision,growth_key = growth_key)
            flux_dict[met_id] = og_me.solution.x_dict
    else:
        import multiprocessing as mp
        
        pool = mp.Pool(NP)
        pbar = tqdm(total=len(met_ids))
        pbar.set_description('{} response ({} threads)'.format(obj,NP))
        def collect_result(result):
            pbar.update(1)
            flux_dict[result[0]] = result[1]

        for met_id in met_ids:
            if sequential:
                met_idx = met_ids.index(met_id)
                met_arg = met_ids[:met_idx+1] # All mets until met
            else:
                met_arg = met_id
            args = (me,met_arg)
            kwds = {'single_change_function':single_change_function,\
                'mu_fix':mu_fix,'precision':precision,'growth_key':growth_key}
            pool.apply_async(single_flux_response,args,kwds,
                callback=collect_result)
        pool.close()
        pool.join()
    flux_results_df = pd.DataFrame.from_dict(flux_dict)
    return flux_results_df

def process_flux_responses(me,flux_results_df,base_rxn):
	base_flux = flux_results_df.loc[base_rxn]['base']
	met_ids = flux_results_df.columns.values
	processed = {}
	for met_id in met_ids:
		if met_id == 'base':
			continue
		processed[met_id] = {}
		new_flux = flux_results_df.loc[base_rxn][met_id]
		met = me.metabolites.get_by_id(met_id)
		met_flux = flux_results_df.loc['dummy_demand'][met_id]
		response = (base_flux-new_flux)/met_flux
		MW = met.formula_weight/1000
		processed[met_id]['molecular_weight'] = MW
		processed[met_id]['sensitivity'] = response

	processed_results_df = pd.DataFrame.from_dict(processed).T

	return processed_results_df

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


def close_transporter(me,transport_id):
        r = me.reactions.get_by_id(transport_id)
        r.lower_bound = 0
        r.upper_bound = 0
        
def overexpress_transporter(me,transport_id):
        r = me.reactions.get_by_id(transport_id)
        base_flux = me.solution.x_dict[r.id]
        
        r.lower_bound = base_flux*2
        r.upper_bound = base_flux*2
        
def group_knockout(me,met_id):
    transport_reactions = get_transport_reactions(me,met_id,comps=['c','s']) \
                + get_transport_reactions(me,met_id,comps=['s','c'])
    
    for r in transport_reactions:
        if 'SPONT' not in r.id:
            r.lower_bound = 0
            r.upper_bound = 0    

def transporter_knockout(me,transport_ids,NP=1,solution=0,precision=1e-6,growth_key='mu',\
                        biomass_dilution='biomass_dilution',single_change_function='transporter',sequential=False):
    '''
    This function calculates the response of shutting down
    transporters.
    ''' 
    
    print('Chosen change function: {}      Sequential = {}'.format(single_change_function,str(sequential)))
        
    flux_results_df = all_flux_responses(me,transport_ids,mu_fix=False,\
            solution=solution,NP=NP,precision=1e-6,\
            single_change_function=single_change_function,growth_key=growth_key,sequential=sequential)

    return flux_results_df