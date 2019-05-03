from __future__ import division, absolute_import, print_function

from cobrame import SubreactionData, Complex
from cobrame.util import dogma

# 1 machine + 1 atp + 1 aa + 1 h2o --> 1 machine-amp + 1 h + 1 ppi
# 1 machine-amp + 1 free tRNA --> 1 machine + 1 amp + 1 charged tRNA
special_trna_subreactions = {
    'sec_addition_at_UGA': {
        'enzymes': ['SelA_deca_mod_10:pydx5p',
                    'SelB_mono'],  # Selenocysteine loaders
        'stoich': {'h_c': 1, 'selnp_c': -1,
                   'pi_c': 1,
                   'generic_tRNA_UGA_cys__L_c': -1},
        'element_contribution': {'O': -1, 'Se': 1}}}

initiation_subreactions = {
    'Translation_initiation_factor_InfA':
        {'enzymes': 'InfA_mono',
         'stoich': {}},

    'Translation_initiation_factor_InfC':
        {'enzymes': 'InfC_mono',
         'stoich': {}},

    'Translation_gtp_initiation_factor_InfB':
        {'enzymes': 'InfB_mono',
         'stoich': {'gtp_c': -1,
                    'h2o_c': -1,
                    'h_c': 1,
                    'pi_c': 1,
                    'gdp_c': 1}},

    'fmet_addition_at_START':
        {'enzymes': ['InfB_mono',
                     'Fmt_mono_mod_mg2_mod_k'],
         # iOL had h_c:1 for fmet addition but this is not mass balanced
         'stoich': {'10fthf_c': -1, 'thf_c': 1,
                    # 'h_c': 1,
                    'generic_tRNA_START_met__L_c': -1},
         'element_contribution': {'C': 1, 'O': 1}}
   }

elongation_subreactions = {'FusA_mono_elongation': {'enzymes': ['FusA_mono'],
                                                    'stoich': {'gtp_c': -1,
                                                               'h2o_c': -1,
                                                               'h_c': 1,
                                                               'pi_c': 1,
                                                               'gdp_c': 1}},

                           'Tuf_gtp_regeneration': {'enzymes': ['Tsf_mono'],
                                                    'stoich': {}}}

# TODO Go through and double check elongation/termination ATP usage etc.
termination_subreactions = {'PrfA_mono_mediated_termination':
                            {'enzymes': ['PrfA_mono'],
                             'stoich': {}},

                            'PrfB_mono_mediated_termination':
                            {'enzymes': ['PrfB_mono'],
                             'stoich': {}},

                            'generic_RF_mediated_termination':
                            {'enzymes': ['generic_RF'],
                             'stoich': {}},

                            'N_terminal_methionine_cleavage':
                            {'enzymes': ['Map_mono_mod_2:fe2'],
                             'stoich': {'h2o_c': -1,
                                        'met__L_c': 1, 'h_c': 1},
                             'element_contribution': {'H': -10, 'O': -1,
                                                      'C': -5, 'N': -1,
                                                      'S': -1}},

                            'peptide_deformylase_processing':
                            {'enzymes': ['Def_mono_mod_1:fe2'],
                             'stoich': {'h2o_c': -1,
                                        'for_c': 1},
                             'element_contribution':
                                 {'H': 1, 'O': -1, 'C': -1}},

                            # This is a GTPS
                            'peptide_chain_release':
                            {'enzymes': ['PrfC_mono'],
                             'stoich': {'gtp_c': -1,
                                        'h2o_c': -1,
                                        'h_c': 1,
                                        'pi_c': 1,
                                        'gdp_c': 1}},

                            'ribosome_recycler':
                            {'enzymes': ['Rrf_mono'],
                             'stoich': {}},

                            'GroEL_dependent_folding':
                            {'enzymes': ['GroL_14', 'cisGroES_hepta',
                                         'transGroES_hepta'],
                             'stoich': {'atp_c': -7,
                                        'h2o_c': -7,
                                        'h_c': 7,
                                        'adp_c': 7,
                                        'pi_c': 7}},

                            # DnaK is correct
                            'DnaK_dependent_folding':
                            {'enzymes': ['DnaK_mono', 'DnaJ_dim_mod_4:zn2',
                                         'GrpE_dim'],
                             'stoich': {'atp_c': -1,
                                        'h2o_c': -1,
                                        'h_c': 1,
                                        'adp_c': 1,
                                        'pi_c': 1}}
                            }

# Subreaction for translation termination
translation_stop_dict = {'UAG': 'PrfA_mono',
                         'UGA': 'PrfB_mono',
                         'UAA': 'generic_RF'}

translation_start_codons = {"AUG", "GUG", "UUG", "AUU", "CUG"}

# Dictionary of frame shift mutations
frameshift_dict = {'b2891': '3033206:3034228,3034230:3034304'}

peptide_processing_subreactions = {"peptide_deformylase_processing",
                                   "peptide_chain_release",
                                   "ribosome_recycler"}


def add_translation_subreactions_to_model(me_model):
    # add general subreactions
    for rxn, info in elongation_subreactions.items():
        data = SubreactionData(rxn, me_model)
        data.enzyme = info['enzymes']
        data.stoichiometry = info['stoich']

    # add subreactions associated with termination and postprocessing
    for rxn, info in termination_subreactions.items():
        data = SubreactionData(rxn, me_model)
        data.enzyme = info['enzymes']
        data.stoichiometry = info['stoich']
        data._element_contribution = info.get('element_contribution', {})

    # add subreactions associated with translation initiation
    for rxn, info in initiation_subreactions.items():
        data = SubreactionData(rxn, me_model)
        data.enzyme = info['enzymes']
        data.stoichiometry = info['stoich']
        data._element_contribution = info.get('element_contribution', {})


def add_charged_trna_subreactions(me_model):
    # create subreaction for each codon. this will be used to model
    # the addition of charged tRNAs to the elongating peptide
    for codon in dogma.codon_table:
        if dogma.codon_table[codon] == '*':
            stop_codon = codon.replace('T', 'U')
            stop_enzyme = translation_stop_dict.get(stop_codon)
            me_model.add_metabolites([Complex(stop_enzyme)])

            subreaction_data = SubreactionData(
                stop_codon + '_' + stop_enzyme + '_mediated_termination',
                me_model)
            subreaction_data.enzyme = stop_enzyme
            subreaction_data.stoichiometry = {}
        else:
            full_aa = dogma.amino_acids[dogma.codon_table[codon]]
            amino_acid = full_aa.split('_')[0]
            subreaction_data = SubreactionData(
                amino_acid + '_addition_at_' + codon.replace('T', 'U'),
                me_model)
            trna = 'generic_tRNA_' + codon.replace('T', 'U') + '_' + full_aa
            subreaction_data.enzyme = 'generic_Tuf'  # Default AA loader enzyme

            # Accounts for GTP hydrolyzed by EF-TU and the ATP hydrolysis to
            # AMP required to add the amino acid to the tRNA
            subreaction_data.stoichiometry = {'gtp_c': -1, 'h2o_c': -2,
                                              'gdp_c': 1, 'h_c': 2, 'pi_c': 1,
                                              'ppi_c': 1, 'amp_c': 1,
                                              'atp_c': -1,
                                              trna: -1}

    # Add subreactions for start codon and selenocysteine
    for rxn, info in special_trna_subreactions.items():
        data = SubreactionData(rxn, me_model)
        data.enzyme = info['enzymes']
        data.stoichiometry = info['stoich']
        data._element_contribution = info.get('element_contribution', {})


# N terminal methionine cleaved
methionine_cleaved = ['BSU17410',
 'BSU10790',
 'BSU32150',
 'BSU17460',
 'BSU06480',
 'BSU39920',
 'BSU06180',
 'BSU16100',
 'BSU01390',
 'BSU28230',
 'BSU16150',
 'BSU31390',
 'BSU15990',
 'BSU01440',
 'BSU30540',
 'BSU01290',
 'BSU18000',
 'BSU06850',
 'BSU33540',
 'BSU13900',
 'BSU03520',
 'BSU28310',
 'BSU02890',
 'BSU25020',
 'BSU06150',
 'BSU01050',
 'BSU00730',
 'BSU19530',
 'BSU01410',
 'BSU23860',
 'BSU34790',
 'BSU33940',
 'BSU01100',
 'BSU04730',
 'BSU30190',
 'BSU23040',
 'BSU03130',
 'BSU25410',
 'BSU32710',
 'BSU32890',
 'BSU28440',
 'BSU29120',
 'BSU28500',
 'BSU01340',
 'BSU01020',
 'BSU16250',
 'BSU23470',
 'BSU30650',
 'BSU38550',
 'BSU33910',
 'BSU29660',
 'BSU28870',
 'BSU07000',
 'BSU06030',
 'BSU16170',
 'BSU16690',
 'BSU01310',
 'BSU01150',
 'BSU01700',
 'BSU00510',
 'BSU01200',
 'BSU01040',
 'BSU30120',
 'BSU14580',
 'BSU33400',
 'BSU04190',
 'BSU19550',
 'BSU12290',
 'BSU19240',
 'BSU01120',
 'BSU38140',
 'BSU21870',
 'BSU37660',
 'BSU00110',
 'BSU18030',
 'BSU25480',
 'BSU40030',
 'BSU13180',
 'BSU01780',
 'BSU35000',
 'BSU28430',
 'BSU08820',
 'BSU31350',
 'BSU01250',
 'BSU16500',
 'BSU04400',
 'BSU07830',
 'BSU36830',
 'BSU27320',
 'BSU13470',
 'BSU08070',
 'BSU00480',
 'BSU13160',
 'BSU33900',
 'BSU23820',
 'BSU01220',
 'BSU14600',
 'BSU18239',
 'BSU15080',
 'BSU29490',
 'BSU40890',
 'BSU28860',
 'BSU29080',
 'BSU16680',
 'BSU39340',
 'BSU16520',
 'BSU01300',
 'BSU01190',
 'BSU01030',
 'BSU18360',
 'BSU02900',
 'BSU25550',
 'BSU16490',
 'BSU16330',
 'BSU01110',
 'BSU16180',
 'BSU01500',
 'BSU07850',
 'BSU08550',
 'BSU00520',
 'BSU26660',
 'BSU37540',
 'BSU14590',
 'BSU25470']

folding_dict = {
    'GroEL_dependent_folding': ['b0014', 'b0015', 'b0061', 'b0062', 'b0064',
                                'b0114', 'b0115', 'b0116', 'b0130', 'b0134',
                                'b0143', 'b0144', 'b0167', 'b0170', 'b0172',
                                'b0185', 'b0209', 'b0369', 'b0404', 'b0439',
                                'b0593', 'b0607', 'b0628', 'b0660', 'b0726',
                                'b0727', 'b0755', 'b0782', 'b0783', 'b0797',
                                'b0870', 'b0902', 'b0929', 'b0930', 'b0957',
                                'b1062', 'b1095', 'b1107', 'b1109', 'b1189',
                                'b1190', 'b1215', 'b1241', 'b1243', 'b1398',
                                'b1718', 'b1719', 'b1748', 'b1779', 'b1831',
                                'b2096', 'b2097', 'b2140', 'b2149', 'b2155',
                                'b2284', 'b2285', 'b2286', 'b2296', 'b2435',
                                'b2441', 'b2478', 'b2533', 'b2551', 'b2557',
                                'b2607', 'b2608', 'b2614', 'b2620', 'b2830',
                                'b2834', 'b2916', 'b2925', 'b2926', 'b2942',
                                'b3162', 'b3256', 'b3260', 'b3295', 'b3340',
                                'b3357', 'b3390', 'b3405', 'b3433', 'b3607',
                                'b3650', 'b3651', 'b3708', 'b3725', 'b3741',
                                'b3775', 'b3780', 'b3845', 'b3847', 'b3850',
                                'b3865', 'b3941', 'b3957', 'b3962', 'b3987',
                                'b3988', 'b3990', 'b3991', 'b3993', 'b3997',
                                'b4039', 'b4154', 'b4177', 'b4381', 'b4382'],
    'DnaK_dependent_folding': ['b2507', 'b2508', 'b2557', 'b2697', 'b2699',
                               'b2764', 'b2780', 'b2913', 'b2925', 'b2926',
                               'b2935', 'b3067', 'b3189', 'b3212', 'b3295',
                               'b3339', 'b3340', 'b3384', 'b3686', 'b3687',
                               'b3708', 'b3744', 'b3783', 'b3829', 'b3831',
                               'b3870', 'b3893', 'b3931', 'b3942', 'b3980',
                               'b3987', 'b3988', 'b4019', 'b4129', 'b4131',
                               'b4147', 'b4177', 'b4232', 'b4239', 'b4258',
                               'b4260', 'b4375', 'b4382', 'b4383', 'b4384',
                               'b4391', 'b0008', 'b0032', 'b0033', 'b0059',
                               'b0095', 'b0114', 'b0115', 'b0118', 'b0194',
                               'b0438', 'b0439', 'b0642', 'b0680', 'b0726',
                               'b0728', 'b0755', 'b0893', 'b0903', 'b0930',
                               'b0932', 'b1014', 'b1095', 'b1114', 'b1136',
                               'b1175', 'b1224', 'b1241', 'b1275', 'b1479',
                               'b1612', 'b1614', 'b1676', 'b1713', 'b1719',
                               'b1779', 'b1945', 'b2029', 'b2036', 'b2114',
                               'b2231', 'b2234', 'b2284', 'b2287', 'b2297',
                               'b2463']}

# Codons are not unique to a tRNA
trna_to_codon = {'BSU_tRNA_1': ['TTT', 'TTC'],
                 'BSU_tRNA_10': ['ATG'],
                 'BSU_tRNA_11': ['GAG', 'GAA'],
                 'BSU_tRNA_12': ['GTA', 'GTC', 'GTG', 'GTT'],
                 'BSU_tRNA_13': ['ACA', 'ACG', 'ACT', 'ACC'],
                 'BSU_tRNA_14': ['AAG', 'AAA'],
                 'BSU_tRNA_15': ['CTT', 'CTG', 'CTA', 'CTC', 'TTA', 'TTG'],
                 'BSU_tRNA_16': ['GGT', 'GGG', 'GGA', 'GGC'],
                 'BSU_tRNA_17': ['CTT', 'CTG', 'CTA', 'CTC', 'TTA', 'TTG'],
                 'BSU_tRNA_18': ['AGG', 'AGA', 'CGA', 'CGG', 'CGT', 'CGC'],
                 'BSU_tRNA_19': ['CCT', 'CCG', 'CCA', 'CCC'],
                 'BSU_tRNA_2': ['GAT', 'GAC'],
                 'BSU_tRNA_20': ['GCA', 'GCC', 'GCG', 'GCT'],
                 'BSU_tRNA_21': ['ATG'],
                 'BSU_tRNA_22': ['GAT', 'GAC'],
                 'BSU_tRNA_23': ['AAC', 'AAT'],
                 'BSU_tRNA_24': ['ACA', 'ACG', 'ACT', 'ACC'],
                 'BSU_tRNA_25': ['GGT', 'GGG', 'GGA', 'GGC'],
                 'BSU_tRNA_26': ['AGG', 'AGA', 'CGA', 'CGG', 'CGT', 'CGC'],
                 'BSU_tRNA_27': ['CCT', 'CCG', 'CCA', 'CCC'],
                 'BSU_tRNA_28': ['GCA', 'GCC', 'GCG', 'GCT'],
                 'BSU_tRNA_29': ['AAC', 'AAT'],
                 'BSU_tRNA_3': ['GAG', 'GAA'],
                 'BSU_tRNA_30': ['AGC', 'AGT', 'TCT', 'TCG', 'TCC', 'TCA'],
                 'BSU_tRNA_31': ['GAG', 'GAA'],
                 'BSU_tRNA_32': ['GTA', 'GTC', 'GTG', 'GTT'],
                 'BSU_tRNA_33': ['ATG'],
                 'BSU_tRNA_34': ['GAT', 'GAC'],
                 'BSU_tRNA_35': ['TTT', 'TTC'],
                 'BSU_tRNA_36': ['ACA', 'ACG', 'ACT', 'ACC'],
                 'BSU_tRNA_37': ['TAT', 'TAC'],
                 'BSU_tRNA_38': ['TGG'],
                 'BSU_tRNA_39': ['CAT', 'CAC'],
                 'BSU_tRNA_4': ['AAG', 'AAA'],
                 'BSU_tRNA_40': ['CAA', 'CAG'],
                 'BSU_tRNA_41': ['GGT', 'GGG', 'GGA', 'GGC'],
                 'BSU_tRNA_42': ['TGT', 'TGC'],
                 'BSU_tRNA_43': ['CTT', 'CTG', 'CTA', 'CTC', 'TTA', 'TTG'],
                 'BSU_tRNA_44': ['CTT', 'CTG', 'CTA', 'CTC', 'TTA', 'TTG'],
                 'BSU_tRNA_45': ['AAC', 'AAT'],
                 'BSU_tRNA_46': ['AGC', 'AGT', 'TCT', 'TCG', 'TCC', 'TCA'],
                 'BSU_tRNA_47': ['GAG', 'GAA'],
                 'BSU_tRNA_48': ['CAA', 'CAG'],
                 'BSU_tRNA_49': ['AAG', 'AAA'],
                 'BSU_tRNA_5': ['ATC', 'ATA', 'ATT'],
                 'BSU_tRNA_50': ['CTT', 'CTG', 'CTA', 'CTC', 'TTA', 'TTG'],
                 'BSU_tRNA_51': ['CTT', 'CTG', 'CTA', 'CTC', 'TTA', 'TTG'],
                 'BSU_tRNA_52': ['GTA', 'GTC', 'GTG', 'GTT'],
                 'BSU_tRNA_53': ['ACA', 'ACG', 'ACT', 'ACC'],
                 'BSU_tRNA_54': ['AAG', 'AAA'],
                 'BSU_tRNA_55': ['CTT', 'CTG', 'CTA', 'CTC', 'TTA', 'TTG'],
                 'BSU_tRNA_56': ['GGT', 'GGG', 'GGA', 'GGC'],
                 'BSU_tRNA_57': ['CTT', 'CTG', 'CTA', 'CTC', 'TTA', 'TTG'],
                 'BSU_tRNA_58': ['AGG', 'AGA', 'CGA', 'CGG', 'CGT', 'CGC'],
                 'BSU_tRNA_59': ['CCT', 'CCG', 'CCA', 'CCC'],
                 'BSU_tRNA_6': ['GCA', 'GCC', 'GCG', 'GCT'],
                 'BSU_tRNA_60': ['GCA', 'GCC', 'GCG', 'GCT'],
                 'BSU_tRNA_61': ['ATG'],
                 'BSU_tRNA_62': ['ATG'],
                 'BSU_tRNA_63': ['AGC', 'AGT', 'TCT', 'TCG', 'TCC', 'TCA'],
                 'BSU_tRNA_64': ['ATG'],
                 'BSU_tRNA_65': ['GAT', 'GAC'],
                 'BSU_tRNA_66': ['TTT', 'TTC'],
                 'BSU_tRNA_67': ['CAT', 'CAC'],
                 'BSU_tRNA_68': ['GGT', 'GGG', 'GGA', 'GGC'],
                 'BSU_tRNA_69': ['ATC', 'ATA', 'ATT'],
                 'BSU_tRNA_7': ['AGC', 'AGT', 'TCT', 'TCG', 'TCC', 'TCA'],
                 'BSU_tRNA_70': ['AAC', 'AAT'],
                 'BSU_tRNA_71': ['AGC', 'AGT', 'TCT', 'TCG', 'TCC', 'TCA'],
                 'BSU_tRNA_72': ['GAG', 'GAA'],
                 'BSU_tRNA_73': ['CAA', 'CAG'],
                 'BSU_tRNA_74': ['CAA', 'CAG'],
                 'BSU_tRNA_75': ['GAG', 'GAA'],
                 'BSU_tRNA_76': ['ACA', 'ACG', 'ACT', 'ACC'],
                 'BSU_tRNA_77': ['TAT', 'TAC'],
                 'BSU_tRNA_78': ['GTA', 'GTC', 'GTG', 'GTT'],
                 'BSU_tRNA_79': ['AGG', 'AGA', 'CGA', 'CGG', 'CGT', 'CGC'],
                 'BSU_tRNA_8': ['ATC', 'ATA', 'ATT'],
                 'BSU_tRNA_80': ['GGT', 'GGG', 'GGA', 'GGC'],
                 'BSU_tRNA_81': ['AGG', 'AGA', 'CGA', 'CGG', 'CGT', 'CGC'],
                 'BSU_tRNA_82': ['GGT', 'GGG', 'GGA', 'GGC'],
                 'BSU_tRNA_83': ['GTA', 'GTC', 'GTG', 'GTT'],
                 'BSU_tRNA_84': ['AGG', 'AGA', 'CGA', 'CGG', 'CGT', 'CGC'],
                 'BSU_tRNA_85': ['AGG', 'AGA', 'CGA', 'CGG', 'CGT', 'CGC'],
                 'BSU_tRNA_86': ['GCA', 'GCC', 'GCG', 'GCT'],
                 'BSU_tRNA_9': ['GCA', 'GCC', 'GCG', 'GCT']}
