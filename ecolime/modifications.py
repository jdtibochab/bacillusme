from __future__ import print_function, absolute_import, division


import cobrame
from ecolime import corrections

# these guys can transfer assembled iron sulfur clusters to the various enzymes
fes_transfer = {}

# Add known specific chaperone
fes_chaperones = {}  # FE-S modification

# complexes that can transfer an iron sulfur cluster to an enzyme target
generic_fes_transfer_complexes = {
}

lipoate_modifications = {
    "mod_lipo_c": {"enzyme": 'BSU07100-MONOMER',
                   "stoich": {"lipoamp_c": -1,
                              "amp_c": 1}},

    "mod_lipo_c_alt": {"enzyme": 'BSU24530-MONOMER',
                       "stoich": {'BSU15920-MONOMER_mod_pan4p_mod_lipo': -1,
                                  'BSU15920-MONOMER_mod_pan4p': 1,
                                  'h_c': -1}}
    }

bmocogdp_chaperones = {}


def add_iron_sulfur_modifications(me_model):

    for name, complexes in generic_fes_transfer_complexes.items():
        generic_fes_transfer = cobrame.GenericData(name, me_model, complexes)
        generic_fes_transfer.create_reactions()

    for fes in ['2fe2s_c', '4fe4s_c']:
        me_model.add_metabolites([cobrame.Metabolite(fes)])
        for name in fes_transfer.values():
            rxn = cobrame.MEReaction('_'.join([name, fes, 'unloading']))
            me_model.add_reactions([rxn])
            rxn.add_metabolites({name + '_mod_1:' + fes.replace('_c', ''): -1,
                                 fes: 1,
                                 name: 1})

    # add fes transfer enzymes to proper modification data
    mod_2fe2s = me_model.process_data.mod_2fe2s_c
    mod_2fe2s.enzyme = 'generic_2fe2s_transfer_complex'
    mod_2fe2s.stoichiometry = {'2fe2s_c': -1.}

    mod_4fe4s = me_model.process_data.mod_4fe4s_c
    mod_4fe4s.enzyme = 'generic_4fe4s_transfer_complex'
    mod_4fe4s.stoichiometry = {'4fe4s_c': -1.}

    for chaperone in set(fes_chaperones.values()):
        new_mod = cobrame.SubreactionData('mod_2fe2s_c_' + chaperone, me_model)
        new_mod.enzyme = [chaperone, 'generic_2fe2s_transfer_complex']
        new_mod.stoichiometry = {'2fe2s_c': -1.}

    for cplx_data in me_model.process_data.get_by_id(
            'mod_2fe2s_c').get_complex_data():
        cplx_id = cplx_data.id.split('_mod')[0]
        if cplx_id in fes_chaperones:
            cplx_data.subreactions['mod_2fe2s_c_' + fes_chaperones[
                cplx_id]] = \
                cplx_data.subreactions.pop('mod_2fe2s_c')


def add_lipoate_modifications(me_model):
    # two different reactions can add a lipoate modification.
    # We create a separate SubreactionData for each one

    for mod, info in lipoate_modifications.items():
        if mod in me_model.process_data:
            mod_data = me_model.process_data.get_by_id(mod)
        else:
            mod_data = cobrame.SubreactionData(mod, me_model)

        mod_data.stoichiometry = info["stoich"]
        mod_data.enzyme = info["enzyme"]

        # element count for lipoate modifications
        mod_data._element_contribution = {'C': 8, 'H': 11, 'O': 1, 'S': 2}

    lipo = me_model.process_data.get_by_id('mod_lipo_c')
    alt_lipo = me_model.process_data.get_by_id('mod_lipo_c_alt')
    for cplx_data in lipo.get_complex_data():
        alt_cplx_data = cobrame.ComplexData(cplx_data.id + "alt", me_model)
        alt_cplx_data.complex_id = cplx_data.complex_id
        alt_cplx_data.stoichiometry = cplx_data.stoichiometry
        alt_cplx_data.subreactions = cplx_data.subreactions.copy()
        alt_cplx_data.subreactions[alt_lipo.id] = \
            alt_cplx_data.subreactions.pop(lipo.id)
        alt_cplx_data.create_complex_formation()


def add_bmocogdp_modifications(me_model):
    for chaperone in set(bmocogdp_chaperones.values()):
        new_mod = \
            cobrame.SubreactionData('mod_bmocogdp_c_' + chaperone, me_model)
        new_mod.enzyme = chaperone
        new_mod.stoichiometry = {'bmocogdp_c': -1}

    for cplx_data in me_model.process_data.get_by_id(
            'mod_bmocogdp_c').get_complex_data():
        cplx_id = cplx_data.id.split('_mod')[0]
        if cplx_id in bmocogdp_chaperones:
            cplx_data.subreactions[
                'mod_bmocogdp_c_' + bmocogdp_chaperones[cplx_id]] = \
                cplx_data.subreactions.pop('mod_bmocogdp_c')


def add_modification_procedures(me_model):
    # add SubreactionData for iron sulfur clusters
    add_iron_sulfur_modifications(me_model)

    # lipoate modifications can be accomplished using two different mechanisms
    add_lipoate_modifications(me_model)

    # bmocogdp modifications have multiple selective chaperones that transfer
    # the metabolite to the target complexes
    add_bmocogdp_modifications(me_model)

    corrections.correct_complex_modifications(me_model)
