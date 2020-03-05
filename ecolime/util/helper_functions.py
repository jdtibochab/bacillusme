import cobrame
from tqdm import tqdm
import pandas as pd

def get_base_complex_data(model, complex_id):
    """If a complex is modified in a metabolic reaction it will not
    have a formation reaction associated with it. This function returns
    the complex data of the "base" complex, which will have the subunit
    stoichiometry of that complex"""

    # First try unmodified complex id
    try_1 = complex_id.split('_')[0]
    if try_1 in model.process_data:
        return model.process_data.get_by_id(try_1)

    try_2 = complex_id.split('_')[0] + '_'
    count = 0
    for i in model.process_data.query(try_2):
        if isinstance(i, cobrame.ComplexData):
            count += 1
            data = i
    if count == 0:
        raise UserWarning('No base complex found for %s' % complex_id)
    if count > 1:
        raise UserWarning('More than one possible base complex found for %s' %
                          complex_id)
    return data

def get_identical_reactions(ref,rxn):
    candidate_rxns = []
    metabolites = rxn.metabolites
    met_ids = []
    # Get reaction metabolite IDs
    for metabolite in metabolites:
        met_ids.append(metabolite.id)
    # Look for identical reactions in reference model
    for ref_rxn in ref.reactions:
        ref_metabolites = ref_rxn.metabolites
        ref_met_ids = []
        if len(metabolites) == len(ref_metabolites):
            for ref_metabolite in ref_metabolites:
                ref_met_ids.append(ref_metabolite.id)
            if len(list(set(ref_met_ids) & set(met_ids))) == len(metabolites):
                candidate_rxns.append(ref_rxn)
        
    return candidate_rxns

def get_gene_info(gb_file,info,ID,element_types):
    output = None
    for feature in gb_file.features:
    # Skip if not a gene used in ME construction
        if feature.type not in element_types or 'pseudo' in feature.qualifiers:
            continue
        if feature.qualifiers["locus_tag"][0] == ID:
            output = feature.qualifiers[info]
    return output

    ############################### NEW ###############################

def test_metabolite_production(me,metabolites,muf = 0.):
	from qminospy.me2 import ME_NLP
	gap_mets = []

	if not muf and me.global_info['k_deg'] != 0:
		print ('Updating model with kdeg = 0 for mu = 0')
		me.global_info['k_deg'] = 0.
		me.update()

	for met_id in metabolites:

		r_id = 'DM_' + met_id
		r = cobrame.MEReaction(r_id)
		try:
			me.add_reaction(r)
			r.reaction = met_id + '->'
		except:
			print(me.reactions.get_by_id(r_id).id,' already in model')
        #print_reactions_of_met(me,met_id)
		me.objective = r_id
        
		me_nlp = ME_NLP(me, growth_key='mu')
		x,status,hs = me_nlp.solvelp(muf)
        
		f = me.solution.x_dict[r_id]

		if not status == 'optimal' or f < 0.01:
			gap_mets.append(met_id)
		print(met_id, status, f)
	return gap_mets

def identify_precursors(me,metabolite_id,only_direct_precursors = False,ignore_classes = None, force_classes = None):
	
	### NOT WORKING YET

	import copy
	precursors = []
	formation_reactions = []
	metabolite = me.metabolites.get_by_id(metabolite_id)
	for rxn in me.reactions:
		if metabolite in rxn.products:
			formation_reactions.append(rxn)
            
	for rxn in formation_reactions:
		for reactant in rxn.reactants:
			if reactant.id not in precursors:
				precursors.append(reactant.id)
	direct_precursors = copy.copy(precursors)
	#print(precursors)
	if not only_direct_precursors:
		for precursor in precursors:
			for rxn in me.metabolites.get_by_id(precursor).reactions:
				products_of_rxn = [product.id for product in rxn.products]
				if precursor in products_of_rxn:
					reactants = rxn.reactants
					reactant_ids = [met.id for met in reactants]
					if metabolite_id not in reactant_ids:            
						for reactant in reactants:
							if reactant.id not in precursors:
								precursors.append(reactant.id)
								#print(reactant.id)
	test_precursors = copy.copy(precursors)
	if ignore_classes:
	    for precursor_id in test_precursors:
	    	precursor = me.metabolites.get_by_id(precursor_id)
	    	for ignore_class in ignore_classes:
		    	if isinstance(precursor,ignore_class):
		    		precursors.remove(precursor_id)
		    		break
	print(len(precursors))
	test_precursors = copy.copy(precursors)
	if force_classes:
		for precursor_id in test_precursors:
			precursor = me.metabolites.get_by_id(precursor_id)
			e = 1
			for force_class in force_classes:
				if isinstance(precursor,force_class):
					e = 0
			if e:
				precursors.remove(precursor_id)	
	print(len(precursors))
	return precursors,direct_precursors

def get_reactions_of_met(me,met,s = 0, ignore_ids = [], verbose = True):
	import copy
	met_stoich = 0

	reactions = []

	all_mets = [met.id for met in me.metabolites]
	if met not in all_mets:
		return reactions

	for rxn in me.metabolites.get_by_id(met).reactions:
		if ignore_ids:
			i = 0
			for identifier in ignore_ids:
				if identifier in rxn.id:
					i = 1
					break
			if i:
				continue

		reactants = [met.id for met in rxn.reactants]
		products = [met.id for met in rxn.products]

		pos = 1 if met in products else -1
		rev = 1 if rxn.lower_bound < 0 else 0
		fwd = 1 if rxn.upper_bound > 0 else 0

		try:
			if not s:
				reactions.append(rxn)
				if verbose:
					print('(',rxn.id,rxn.lower_bound,rxn.upper_bound,')', '\t',rxn.reaction)

			elif s == pos*fwd or s == pos*rev:
				reactions.append(rxn)
				if verbose:
					print('(',rxn.id,rxn.lower_bound,rxn.upper_bound,')', '\t',rxn.reaction)
			
		except:
			if verbose:
				print(rxn.id, 'no reaction')
			else:
				pass
	return reactions

def add_exchange_reactions(me,metabolites):
	for met in metabolites:
		rxn_id = "EX_" + met
		try:
			r = cobrame.MEReaction(rxn_id)
			me.add_reaction(r)
			r.reaction = met + " <=> "
		except:
			r = me.reactions.get_by_id(rxn_id)
		r.lower_bound = -1000
		#print(r.id,r.lower_bound,r.upper_bound,r.reaction)
	return me

def brute_force_check(me,metabolites_to_add,objective_function = 'biomass_dilution',muf = 0.01, min_f = 0.01):

	me.objective = objective_function

	from qminospy.me2 import ME_NLP
	print('Added exchange reactions ')
	me = add_exchange_reactions(me,metabolites_to_add)
	print('Objective: ', objective_function, me.reactions.get_by_id(objective_function).reaction)

	me_nlp = ME_NLP(me, growth_key='mu')
	x,status,hs = me_nlp.solvelp(muf)
	initial_f = me.solution.x_dict[objective_function]
	print('Initial objective function value of ', initial_f, status)

	if not status =='optimal':
		return

	if initial_f < min_f:
		print('No production capacity of objective')

	print(me.solution.x_dict['formation_ribosome'])

	eliminate_mets = []
	for met_id in metabolites_to_add:
		ex_rxn_id = "EX_" + met_id
		ex_rxn_flux = me.solution.x_dict[ex_rxn_id]
		ex_rxn = me.reactions.get_by_id(ex_rxn_id)
		if ex_rxn_flux > 0:
			me.reactions.get_by_id(ex_rxn_id).lower_bound = 0
			me.reactions.get_by_id(ex_rxn_id).upper_bound = 1000
			print(ex_rxn_id, ex_rxn_flux, ex_rxn.reaction)
		elif ex_rxn_flux < 0:
			me.reactions.get_by_id(ex_rxn_id).lower_bound = -1000
			me.reactions.get_by_id(ex_rxn_id).upper_bound = 0
			print(ex_rxn_id, ex_rxn_flux, ex_rxn.reaction)
		elif ex_rxn_flux == 0:
			me.reactions.get_by_id(ex_rxn_id).lower_bound = 0
			me.reactions.get_by_id(ex_rxn_id).upper_bound = 0
			print(ex_rxn_id, ' carrying no flux ... eliminated')

			eliminate_mets.append(met_id)


	for el_met_id in eliminate_mets:
		el_rxn_id = 'EX_' + el_met_id
		metabolites_to_add.remove(el_met_id)

	print('Processing ', len(metabolites_to_add), ' metabolites')

	gap_mets = []
	for met_id in metabolites_to_add:
		ex_rxn_id =  "EX_" + met_id

		lb = me.reactions.get_by_id(ex_rxn_id).lower_bound
		ub = me.reactions.get_by_id(ex_rxn_id).lower_bound

		me.reactions.get_by_id(ex_rxn_id).lower_bound = 0
		me.reactions.get_by_id(ex_rxn_id).upper_bound = 0

		me_nlp = ME_NLP(me, growth_key='mu')
		x,status,hs = me_nlp.solvelp(muf)
		f = me.solution.x_dict[objective_function]

		el_bool = ''
		if not status == 'optimal' or f < min_f:
			me.reactions.get_by_id(ex_rxn_id).lower_bound = lb
			me.reactions.get_by_id(ex_rxn_id).lower_bound = ub
			gap_mets.append(met_id)
			el_bool = ' gap'
            
		print(met_id, status, f, el_bool, '... Gaps: ', len(gap_mets))
	return gap_mets

def solve_me_model(me, max_mu=1., precision=1e-6, min_mu=0, using_soplex=True,
                  compiled_expressions=None, verbosity = 2, mu_fix = False):
	from qminospy.me1 import ME_NLP1

	## If fixed growth rate, solve as LP
	if mu_fix:
		me_nlp = ME_NLP1(me)
		me_nlp.solvelp(mu_fix)
	else:
		## Full solve. USE ME_NLP.solvenlp???
		if using_soplex:
			from cobrame.solve.algorithms import binary_search
			binary_search(me, min_mu=min_mu, max_mu=max_mu, debug=True, mu_accuracy=precision,
				compiled_expressions=compiled_expressions)
		else:
			# The object containing solveME methods--composite that uses a ME model object 
			me_nlp = ME_NLP1(me, growth_key='mu')
			# Use bisection for now (until the NLP formulation is worked out)
			muopt, hs, xopt, cache = me_nlp.bisectmu(precision=precision, mumax=max_mu, verbosity=verbosity)
			try:
				me.solution.f = me.solution.x_dict['biomass_dilution']
			except:
				pass

def show_escher_map(me, solution=None):
    import escher
    view = escher.Builder("iJO1366.Central metabolism")
    view.reaction_data = me.get_metabolic_flux(solution=solution)
    return view

def open_all_exchange(me):
	for rxn in me.reactions:
		rxn_id = rxn.id
		if 'EX_' in rxn_id:
			rxn.upper_bound = 1000
			rxn.lower_bound = -1000
	return me

def is_same_reaction(rxn,ref_rxn):
	met_ids = [met.id for met in rxn.metabolites]
	ref_met_ids = [met.id for met in ref_rxn.metabolites]

	i = 0
	if (len(met_ids) == len(ref_met_ids)) and (len(set(met_ids) & set(ref_met_ids)) == len(met_ids)):
		i = 1

	return i

def homogenize_reactions(model,ref_model):
	all_ref_rxns = [rxn.id for rxn in ref_model.reactions] 
	rxn_dict = dict()
	rxn_id_dict = {}
	for rxn in tqdm(model.reactions):
		if rxn.id in all_ref_rxns:
			ref_rxn = ref_model.reactions.get_by_id(rxn.id)
			# Check if rxn has same ID in ref_model to avoid iterating
			if is_same_reaction(rxn,ref_rxn):
				rxn_dict[rxn] = ref_rxn
				rxn_id_dict[rxn.id] = ref_rxn.id
			# If rxn is not in ref_model, iterate
			else:
				for ref_rxn in ref_model.reactions:
					if is_same_reaction(rxn,ref_rxn):
						rxn_dict[rxn] = ref_rxn
						rxn_id_dict[rxn.id] = ref_rxn.id
	return rxn_dict,rxn_id_dict

def exchange_single_model(me, solution = 0):
	import pandas as pd

	complete_dict = {'id':[],'name':[],'reaction':[],'lower_bound':[],'upper_bound':[],'flux':[]}
	

	if not solution:
		solution = me.solution

	for rxn in me.reactions:
		if 'EX_' in rxn.id:
			flux = solution.x_dict[rxn.id]

			if not flux:
				continue
			rxn_name = rxn.name
			reaction = rxn.reaction
			lb = rxn.lower_bound
			ub = rxn.upper_bound
			
			
			complete_dict['id'].append(rxn.id)
			complete_dict['name'].append(rxn_name)
			complete_dict['reaction'].append(reaction)
			complete_dict['lower_bound'].append(lb)
			complete_dict['upper_bound'].append(ub)
			complete_dict['flux'].append(flux)


	df = pd.DataFrame(complete_dict).set_index('id')
	return df

def get_metabolites_from_pattern(model,pattern):
	met_list = []
	for met in model.metabolites:
		if pattern in met.id:
			met_list.append(met.id)
	return met_list

def flux_based_reactions(model,met_id,s,ignore_ids=[],threshold = 0.2,solution=0):
	if not solution:
		solution = model.solution

	flux_dict = solution.x_dict
	reactions = get_reactions_of_met(model,met_id,s,ignore_ids,False)

	met_stoich = {}
	for rxn in reactions:
		for rxn_met,stoich in rxn.metabolites.items():
			if rxn_met.id == met_id:
				if hasattr(stoich, 'subs'):
					met_stoich[rxn.id] = float(stoich.subs('mu',solution.f))
				else:
					met_stoich[rxn.id] = stoich
				break

	fluxes = [flux_dict[rxn.id]*met_stoich[rxn.id] for rxn in reactions]
	max_precursor_flux = max(fluxes)

	reactions_with_flux = []
	for rxn in reactions:
		flux = flux_dict[rxn.id]*met_stoich[rxn.id]
		if abs(flux) > abs(threshold*max_precursor_flux):
			print(rxn.id,'(',flux,')',rxn.reaction)
			reactions_with_flux.append(rxn)
	return reactions_with_flux

def generate_gene_field(me,identifier):
	import cobra
	current_gene_ids = [gene.id for gene in me.genes]
	for met in me.metabolites:
	    met_id = met.id
	    if identifier in met_id:
	        gene_id = met_id.split('_')[1]
	        if gene_id and gene_id not in current_gene_ids:
	            try:
	                gene = cobra.Gene(gene_id)
	                me.genes.append(gene)
	                print(gene_id)
	            except:
	                pass

	return me

def solution_summary(me):
	reactions = [rxn.id for rxn in me.reactions]
	summary_df = pd.DataFrame(columns=['lb','ub','flux','formula'],index=reactions)

	for rxn_id in tqdm(reactions):
		rxn = me.reactions.get_by_id(rxn_id)
		summary_df.loc[rxn_id]['lb'] = rxn.lower_bound
		summary_df.loc[rxn_id]['ub'] = rxn.upper_bound
		summary_df.loc[rxn_id]['flux'] = me.solution.x_dict[rxn_id]
		summary_df.loc[rxn_id]['formula'] = rxn.reaction
		
	return summary_df

def get_flux_for_escher(model,type='m'):

	if type == 'm':
		flux_dict = model.solution.x_dict
	elif type == 'me':
		flux_dict = model.get_metabolic_flux(solution=me.solution)

	return pd.DataFrame.from_dict({'flux':flux_dict})

