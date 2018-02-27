from global_params import MMU

IGK_id = 'RTVAA'
IGL_id = 'VTVLL'
if MMU:
    IGK_id = 'RADAA'

iso_A_id = 'ASPTSP'
iso_A1_id = 'VVIACL'
iso_A2_id = 'VVVACL'
iso_G_id = 'ASTKGPS'
iso_M_id = 'SASAPT'

end_J_seq = 'VTVSS'
# mass_spec_seq = 'ASTK'
mass_spec_seq = ''

if MMU:
    iso_A_id = 'ESARNP'
    iso_G_id = 'AKTTPP'
    iso_M_id = 'ESQSFP'
    mass_spec_seq = 'AK'

adapter_seq = 'PPLIP'