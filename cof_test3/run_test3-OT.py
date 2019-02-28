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
        "num_machines": 6,
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
            'OPTIMIZER': 'LBFGS',                    #default: BFGS
            'LBFGS' : {
                'TRUST_RADIUS': '[angstrom] 0.5',     #default: None
            },
            'MAX_ITER': 1000,
            'KEEP_ANGLES' : False,
            'MAX_DR':    '[bohr] 0.030',               #default: [bohr] 0.0030
            'RMS_DR':    '[bohr] 0.015',               #default: [bohr] 0.0015
            'MAX_FORCE': '[bohr^-1*hartree] 0.00045',
            'RMS_FORCE': '[bohr^-1*hartree] 0.00003',
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
with open('../cof_test2/list-OT.list') as f:
    ids=f.read().splitlines()
ids=ids[:3]
prevWorkflow = 'test2-0'
for id in ids:
    print('SUBMITTING: %s '%id)
    q = QueryBuilder()
    q.append(StructureData, filters={'label': id}, tag='inp_struct')
    q.append(WorkCalculation, filters={'label':prevWorkflow},
                              output_of='inp_struct', tag='wf')
    q.append(CifData, edge_filters={'label': 'output_structure'},
                      output_of='wf')
    q.order_by({WorkCalculation:{'ctime':'desc'}})
    cif = q.all()[0][0] #take the last
    structure = _get_aiida_structure_ase_inline(cif)['structure']
    structure.store()
    structure.label = id
    # Submit the calculation
    submit(Cp2kDirectCellOptDdecWorkChain,
        structure=structure,
        cp2k_code=cp2k_code,
        cp2k_parameters=cp2k_parameters,
        _cp2k_options=cp2k_options,
        ddec_compute_charges=Bool(True),
        ddec_code=ddec_code,
        ddec_parameters=ddec_parameters,
        _ddec_options=ddec_options,
        _label='test3-OT',
        _guess_multiplicity=True,
        min_cell_size=Float(5.0)
        )
