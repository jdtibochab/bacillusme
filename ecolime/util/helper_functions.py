import cobrame

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

def print_reactions_of_met(me,met,s = 0):
	import copy
	met_stoich = 0
	for rxn in me.metabolites.get_by_id(met).reactions:
		for key,value in rxn.metabolites.items():
			if str(key) == str(met):
				met_stoich = copy.copy(value)
		if s == 1 and met_stoich > 0:
			print('(',rxn.id,rxn.lower_bound,rxn.upper_bound,')', '\t',rxn.reaction)
		elif s == -1 and met_stoich < 0:
			print('(',rxn.id,rxn.lower_bound,rxn.upper_bound,')', '\t',rxn.reaction)
		elif s == 0:
			print('(',rxn.id,rxn.lower_bound,rxn.upper_bound,')', '\t',rxn.reaction)
        
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

def solve_me_model(me, max_mu, precision=1e-6, min_mu=0, using_soplex=True,
                  compiled_expressions=None):
    if using_soplex:
        from cobrame.solve.algorithms import binary_search
        binary_search(me, min_mu=min_mu, max_mu=max_mu, debug=True, mu_accuracy=precision,
                      compiled_expressions=compiled_expressions)
    else:
        from qminospy.me1 import ME_NLP1
        # The object containing solveME methods--composite that uses a ME model object 
        me_nlp = ME_NLP1(me, growth_key='mu')
        # Use bisection for now (until the NLP formulation is worked out)
        muopt, hs, xopt, cache = me_nlp.bisectmu(precision=precision, mumax=max_mu)
        me.solution.f = me.solution.x_dict['biomass_dilution']

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

def homogenize_reactions(model,ref_model):
	rxn_dict = dict()
	for rxn in model.reactions:
		for ref_rxn in ref_model.reactions:
			mets = rxn.metabolites
			ref_mets = ref_rxn.metabolites
			if (len(mets) == len(ref_mets)): # and (len(set(mets) & set(ref_mets)) == len(mets)):
				rxn_dict[rxn] = ref_rxn
	return rxn_dict
