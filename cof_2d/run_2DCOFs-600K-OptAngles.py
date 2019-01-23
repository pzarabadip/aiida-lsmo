from aiida.common.example_helpers import test_and_get_code  # noqa
from aiida.orm.data.base import Float
from aiida.work.run import submit
from ase.io import read
from glob import glob
from os import path
from aiida_lsmo_workflows.geoopt_charges import Cp2kGeoOptDdecWorkChain
ParameterData = DataFactory('parameter')
StructureData = DataFactory('structure')

# Test the codes and specify the nodes and walltime
cp2k_code = test_and_get_code('cp2k-5.1@fidis', expected_code_type='cp2k')
ddec_code = test_and_get_code('ddec@fidis', expected_code_type='ddec')
cp2k_options = {
    "resources": {
        "num_machines": 2,
    },
    "max_wallclock_seconds": 24 * 60 * 60,
    }
ddec_options = {
    "resources": {
        "num_machines": 1,
    },
    "max_wallclock_seconds": 5 * 60 * 60,
    "withmpi": False,
    }

# Set the settings for CP2K

''' First CELL_OPT (Stage1), hardcoded
params_dict = {
        'MOTION':{
            'CELL_OPT': {
                'MAX_ITER': 20,
                'KEEP_ANGLES' : True,
            },
        },
}
'''
params_dict = {
        'MOTION':{
            'MD':{
                'STEPS': 50,
                'TEMPERATURE': '[K] 600',
                },
            'GEO_OPT': {
                'MAX_ITER': 20,
            },
            'CELL_OPT': {
                'MAX_ITER': 140,
                'KEEP_ANGLES' : False,
            },
        },
}
cp2k_parameters = ParameterData(dict=params_dict)

# Using lists to specify the IDs
with open('2dN.list') as f:
    ids=f.read().splitlines()
all_structures = [ "/home/daniele/Documents/CoRE-COFs/cifs/{}.cif".format(x) for x in ids]
# Submit the calculations
for s in all_structures:
    s_ase = read(s)
    structure = StructureData(ase=s_ase)
    structure.label = s.split('/')[-1].split('.')[0]
    structure.store()
    submit(Cp2kGeoOptDdecWorkChain,
        structure=structure,
        cp2k_code=cp2k_code,
        cp2k_parameters=cp2k_parameters,
        _cp2k_options=cp2k_options,
        ddec_code=ddec_code,
        _ddec_options=ddec_options,
        _label='2DCOFs-600K-OptAngles',
        _guess_multiplicity=True,
        min_cell_size=Float(5.0)
        )
