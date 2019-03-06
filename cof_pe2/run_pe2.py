import os
from aiida.common.example_helpers import test_and_get_code
from aiida.orm import DataFactory
from aiida.orm.data.base import Float
from aiida.orm.calculation.work import WorkCalculation
from aiida.work.run import submit
from aiida_lsmo_workflows.volpo_Kh_isotherm import VolpoKhIsothermWorkChain
ParameterData = DataFactory('parameter')
StructureData = DataFactory('structure')
SinglefileData = DataFactory('singlefile')
CifData = DataFactory('cif')

def dict_merge_ez(dict1, dict2):
    sumdicts = dict1.copy()
    sumdicts.update(dict2)
    return sumdicts

# Test the codes and specify the nodes and walltime
zeopp_code = test_and_get_code('zeopp@deneb', expected_code_type='zeopp.network')
raspa_code = test_and_get_code('raspa@deneb', expected_code_type='raspa')

zeopp_options = {
    "resources": {
        "num_machines": 1,
        "tot_num_mpiprocs": 1,
    },
    "max_wallclock_seconds": 10 * 60 * 60,
    "withmpi": False,
    }
raspa_options = {
    "resources": {
        "num_machines": 1,
        "tot_num_mpiprocs": 1,
    },
    "max_wallclock_seconds": 72 * 60 * 60,
    "withmpi": False,
    }

# Settings for Zeopp and Raspa (Widom calculation)
zeopp_probe_radius_co2_trappe = Float(1.525)
zeopp_probe_radius_n2_trappe = Float(1.655)
zeopp_atomic_radii_file = SinglefileData(file=os.path.abspath("./UFF.rad"))
raspa_params_general = {
        "GeneralSettings":
        {
        "NumberOfInitializationCycles"     : 1000, # Widom will use 0
        "NumberOfCycles"                   : 10000, # Widom will use 10x
        "Forcefield"                       : "LSMO_UFF-TraPPE",
        "CutOff"                           : 12.0,
        "ExternalTemperature"              : 300.0,
        },
}
raspa_params_co2 = {
        "Component":
        [{
        "MoleculeName"                     : "CO2",
        "MoleculeDefinition"               : "TraPPE",
        }],
}
raspa_params_n2 = {
        "Component":
        [{
        "MoleculeName"                     : "N2",
        "MoleculeDefinition"               : "TraPPE",
        }],
}
raspa_minKh_co2 = Float(0) #(mol/kg/Pa) Use GCMC if Kh>minKh
raspa_molsatdens_co2 = Float(21.2) #(mol/l) Density of the molecule @ saturation
raspa_minKh_n2 = Float(0) #(mol/kg/Pa) Use GCMC if Kh>minKh
raspa_molsatdens_n2 = Float(28.3) #(mol/l) Density of the molecule @ saturation

# Take the structures from a RobustGeoOptDdec calculation and submit
with open('../cof_test2/list-OT.list') as f:
    ids=f.read().splitlines()
ids=ids[:50]
ids=['05001N2','07000N2','07013N3','11001N2','13011N2','13040N3','13072N2'] #6Mar, only N2, resubmitted after they crashed because of the "None" problem
prevwf_label = 'test2-0'
for id in ids:
    q = QueryBuilder()
    q.append(StructureData, filters={'label': id}, tag='inp_struct')
    q.append(WorkCalculation, filters={'label':prevwf_label},
                              output_of='inp_struct', tag='wf')
    q.append(CifData, edge_filters={'label': 'output_structure'},
                      output_of='wf')
    q.order_by({WorkCalculation:{'ctime':'desc'}})
    if len(q.all())==0:
        print('WARNING: {} not found from prevous workflow ({})'.format(id,prevwf_label))
    else:
        print('SUBMITTING: {} from prevous workflow ({})'.format(id,prevwf_label))
        structure = q.all()[0][0] #take the last
        # Run for CO2, using UFF-TraPPE force field
        submit(VolpoKhIsothermWorkChain,
            structure=structure,
            zeopp_code=zeopp_code,
            _zeopp_options=zeopp_options,
            zeopp_probe_radius=zeopp_probe_radius_co2_trappe,
            zeopp_atomic_radii=zeopp_atomic_radii_file,
            raspa_code=raspa_code,
            raspa_parameters=ParameterData(dict=dict_merge_ez(raspa_params_general,raspa_params_co2)),
            _raspa_options=raspa_options,
            _raspa_usecharges=True,
            raspa_minKh=raspa_minKh_co2,
            raspa_molsatdens=raspa_molsatdens_co2,
            _label='pe2-co2',
            )
        # Run for N2, using UFF-TraPPE force field
        submit(VolpoKhIsothermWorkChain,
            structure=structure,
            zeopp_code=zeopp_code,
            _zeopp_options=zeopp_options,
            zeopp_probe_radius=zeopp_probe_radius_n2_trappe,
            zeopp_atomic_radii=zeopp_atomic_radii_file,
            raspa_code=raspa_code,
            raspa_parameters=ParameterData(dict=dict_merge_ez(raspa_params_general,raspa_params_n2)),
            _raspa_options=raspa_options,
            _raspa_usecharges=True,
            raspa_minKh=raspa_minKh_n2,
            raspa_molsatdens=raspa_molsatdens_n2,
            _label='pe2-n2',
            )
