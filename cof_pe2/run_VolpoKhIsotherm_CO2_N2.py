import os
import pandas as pd
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
zeopp_code = test_and_get_code('zeopp@fidis', expected_code_type='zeopp.network')
raspa_code = test_and_get_code('raspa2@fidis', expected_code_type='raspa')

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
prevwf_label_dict={'OT':'test2-0','smear':'test2-smearing'}
df=pd.read_csv("../cof_test2/pk_final.csv")
df = df[df.extension_ok != 1].reset_index()

for i in range(len(df)):
    structure_label = df.at[i,'structure']
    prevwf_label = prevwf_label_dict[df.at[i,'dft']]
    q = QueryBuilder()
    q.append(StructureData, filters={'label':structure_label}, tag='inp_struct')
    q.append(WorkCalculation, filters={'label':prevwf_label},
                              output_of='inp_struct', tag='wf')
    q.append(CifData, edge_filters={'label': 'output_structure'},
                      output_of='wf')
    q.order_by({WorkCalculation:{'ctime':'desc'}})
    if len(q.all())==0:
        print('WARNING: label={} not found from prevous workflow (label={})'.format(
              structure_label,prevwf_label))
    else:
        structure = q.all()[0][0] #take the last
        print('SUBMITTING: label={} (pk={}) from prevous workflow (label={})'.format(
              structure_label,structure.pk,prevwf_label))

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
            _label='pe3-co2-fidis',
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
            _label='pe3-n2-fidis',
            )
