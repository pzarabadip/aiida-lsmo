from aiida.common.example_helpers import test_and_get_code  # noqa
from aiida.orm.data.base import Float, Bool
from aiida.orm.data.cif import _get_aiida_structure_ase_inline
from aiida.orm.calculation.work import WorkCalculation
from aiida.work.run import submit
from ase.io import read
from glob import glob
from os import path
from aiida_lsmo_workflows.directcellopt_charges import Cp2kDirectCellOptDdecWorkChain
ParameterData = DataFactory('parameter')
StructureData = DataFactory('structure')
CifData = DataFactory('cif')

# Test the codes and specify the nodes and walltime
cp2k_code = test_and_get_code('cp2k_6.1_18464@daint-s888', expected_code_type='cp2k')
ddec_code = test_and_get_code('ddec@daint-s888', expected_code_type='ddec')

cp2k_options = {
    "resources": {
        "num_machines": 4,
    },
    "max_wallclock_seconds": 24 * 60 * 60,
    }
ddec_options = {
    "resources": {
        "num_machines": 1,
    },
    "max_wallclock_seconds": 10 * 60 * 60,
    "withmpi": False,
    }

# Set the settings for CP2K
params_dict = {
    'FORCE_EVAL': {
        'DFT': {
            'SCF': {
                'EPS_SCF': 1.0E-7,
                'OUTER_SCF': {
                    'EPS_SCF': 1.0E-7,
                    },
                },
            },
        },
    'MOTION': {
        'CELL_OPT': {
            'OPTIMIZER': 'BFGS',                       #CP2K default
            'MAX_ITER': 1000,
            'KEEP_ANGLES' : False,
            'PRESSURE_TOLERANCE' : '[bar] 100',        #CP2K default
            'MAX_DR':    '[bohr] 0.030',               #CP2K default
            'RMS_DR':    '[bohr] 0.015',               #CP2K default
            'MAX_FORCE': '[bohr^-1*hartree] 0.00045',  #CP2K default
            'RMS_FORCE': '[bohr^-1*hartree] 0.00003',  #CP2K default
        },
    },
}
cp2k_parameters = ParameterData(dict=params_dict)

# Set the data folder for DDEC (depending on the computer!)
params_dict = {
    "net charge"                               : 0.0,
    "charge type"                              : "DDEC6",
    "periodicity along A, B, and C vectors"    : [True, True, True,],
    "compute BOs"                              : False,
    "input filename"                           : "valence_density",
    "atomic densities directory complete path" : "/users/ongari//aiida-database/data/chargemol_09_26_2017/atomic_densities/",
}

ddec_parameters = ParameterData(dict=params_dict)


# Using lists to specify the IDs and get the output cif from test2
with open('list-directcellopt.list') as f:
    ids=f.read().splitlines()
all_structures = [ "/home/daniele/Documents/CoRE-COFs/cifs/{}.cif".format(x) for x in ids]
# Submit the calculations
for s in all_structures:
    print('SUBMITTING: %s '%s)
    s_ase = read(s)
    structure = StructureData(ase=s_ase)
    structure.label = s.split('/')[-1].split('.')[0]
    structure.store()
    submit(Cp2kDirectCellOptDdecWorkChain,
            structure=structure,
            cp2k_code=cp2k_code,
            cp2k_parameters=cp2k_parameters,
            _cp2k_options=cp2k_options,
            ddec_compute_charges=Bool(False),
            ddec_code=ddec_code,
            ddec_parameters=ddec_parameters,
            _ddec_options=ddec_options,
            _label='test4-OT',
            _guess_multiplicity=True,
            min_cell_size=Float(5.0)
            )
