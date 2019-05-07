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