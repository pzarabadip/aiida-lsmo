from aiida.common.example_helpers import test_and_get_code  # noqa
from aiida.orm.data.structure import StructureData  # noqa
from aiida.orm.data.parameter import ParameterData  # noqa
from aiida.orm.data.base import Float
from aiida.work.run import submit
from ase.io import read

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
params_dict = {
        'MOTION':{
            'MD':{
                'STEPS': 20,
                },
            'GEO_OPT': {
                'MAX_ITER': 20,
            },
            'CELL_OPT': {
                'MAX_ITER': 140,
            },
        },
}
cp2k_parameters = ParameterData(dict=params_dict)

# Specify the path of the cif files to import, label and submit them
from glob import glob
from os import path
from aiida_lsmo_workflows.geoopt_charges import Cp2kGeoOptDdecWorkChain
ArrayData = DataFactory('array')
ParameterData = DataFactory('parameter')
StructureData = DataFactory('structure')
all_structures = glob(path.abspath("/home/daniele/Documents/CoRE-COFs/cifs/13*N.cif")) #18 COFs from 2005 to 2009
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
        _label='Test18CofsNoFixAngles',
        _guess_multiplicity=True,
        min_cell_size=Float(5.0)
        )
