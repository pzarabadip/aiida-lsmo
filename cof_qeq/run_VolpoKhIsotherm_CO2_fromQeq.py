import os
import pandas as pd
from aiida_qeq.calculations.qeq import QeqCalculation
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

# setup codes and parameters
zeopp_code = test_and_get_code('zeopp@deneb', expected_code_type='zeopp.network')
raspa_code = test_and_get_code('raspa2@deneb', expected_code_type='raspa')

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


df = pd.read_csv("../cof_test2/pk_final.csv")
wc_label = "aiida_qeq Qeq with GMP parameters"

for i in df.index:
 if i>1:
    cif_label = df.at[i,'structure']
    q = QueryBuilder()
    q.append(CifData, filters={'label': cif_label}, tag='cif-inp')
    q.append(QeqCalculation, filters={'label': wc_label}, output_of='cif-inp', tag='qeq')
    if len(q.all())>0:
        pk_qeqcalc=q.all()[0][0].pk
    else:
        pk_cifqeq='missing'
    q.append(CifData, output_of='qeq', tag='cif-qeq')
    if len(q.all())>0:
        cifqeq = q.all()[0][0]
        pk_cifqeq = cifqeq.pk
        print("CifData.label= {} with Qeq >>> CifData.pk= {} ".format(
              cif_label,pk_cifqeq))
    else:
        print("CifData.label= {} failed QeqCalculation.pk= {}".format(
              cif_label,pk_qeqcalc))
        break

    # Submit the kH calculation if cif found: I only need blockpore + kH
    submit(VolpoKhIsothermWorkChain,
        structure=cifqeq,
        zeopp_code=zeopp_code,
        _zeopp_options=zeopp_options,
        zeopp_probe_radius=Float(1.525), #CO2
        zeopp_atomic_radii=SinglefileData(file=os.path.abspath("../cof_pe2/UFF.rad")),
        raspa_code=raspa_code,
        raspa_parameters=ParameterData(dict=dict_merge_ez(raspa_params_general,raspa_params_co2)),
        _raspa_options=raspa_options,
        _raspa_usecharges=True,
        raspa_minKh=Float(1e10), # *** HIGH to disable GCMC ***
        raspa_molsatdens=Float(21.2), #(mol/l) for CO2,
        _label='qeq-kh-co2',
        )
