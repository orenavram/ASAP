from global_params import MMU, add_mass_spec_seq


IGK_id = 'RTVAA'
IGL_id = 'VTVLL'

isotype_A_id = 'ASPTSP'
isotype_A1_id = 'ASPTSPKVFPLSLCSTQP'
isotype_A2_id = 'ASPTSPKVFPLSLDSTPQ'
isotype_D_id = 'APTKAP'
isotype_E_id = 'ASTQSP'
isotype_G_id = 'ASTKGPS'
isotype_M_id = 'GSASAPT'

end_j_seq = 'VTVSS' #framework 4

mass_spec_seq = 'BLABLABLA' #'ASTK'

if MMU:
    IGK_id = 'RADAA'
    isotype_A_id = 'ESARNP'
    isotype_D_id = 'GDKKEP'
    isotype_E_id = 'ASIRNP'
    isotype_G_id = r'AKTT[AP]P'
    isotype_M_id = 'ESQSFP'
    #todo: handle throught the code
    #end_j_seq = ('VTVSS', 'LTVSS', 'VTVSA')  # framework 4
    mass_spec_seq = 'AK'

if not add_mass_spec_seq:
    # when the user does not want to add any mass_spec_seq to the final fasta file
    mass_spec_seq = ''
