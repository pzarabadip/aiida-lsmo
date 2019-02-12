from aiida.common.example_helpers import test_and_get_code  # noqa
from aiida.orm.data.base import Float
from aiida.work.run import submit
from ase.io import read
from glob import glob
from os import path
from aiida_lsmo_workflows.directcellopt import Cp2kDirectCellOptWorkChain
ParameterData = DataFactory('parameter')
StructureData = DataFactory('structure')

# Test the codes and specify the nodes and walltime
cp2k_code = test_and_get_code('cp2k-5.1@fidis', expected_code_type='cp2k')
ddec_code = test_and_get_code('ddec@fidis', expected_code_type='ddec')
cp2k_options = {
    "resources": {
        "num_machines": 4,
    },
    "max_wallclock_seconds": 24 * 60 * 60,
    }

# Set the settings for CP2K
params_dict = {
    'MOTION': {
        'CELL_OPT': {
            'MAX_ITER': 1000,
            'KEEP_ANGLES' : False,
        },
    },
}

cp2k_parameters = ParameterData(dict=params_dict)

# Using lists to specify the IDs
with open('list-directcellopt.list') as f:
    ids=f.read().splitlines()
all_structures = [ "./cifs/{}.cif".format(x) for x in ids]
# Submit the calculations
for s in all_structures:
    print('SUBMITTING: %s '%s)
    s_ase = read(s)
    structure = StructureData(ase=s_ase)
    structure.label = s.split('/')[-1].split('.')[0]
    structure.store()
    submit(Cp2kDirectCellOptWorkChain,
        structure=structure,
        cp2k_code=cp2k_code,
        cp2k_parameters=cp2k_parameters,
        _cp2k_options=cp2k_options,
        _label='directcellopt_randomized',
        _guess_multiplicity=True,
        min_cell_size=Float(10.0)
        )
