import os
from aiida.plugins import DataFactory
from aiida.orm import Code, Dict, Float, Int
from aiida.engine import run, submit
from aiida_lsmo_workflows.separation_2comp import SeparationWorkChain

ParameterData = DataFactory('dict')
SinglefileData = DataFactory('singlefile')
CifData = DataFactory('cif')

# Import the structure
structure_zeopp = CifData(file=os.path.abspath("./HKUST1.cif"))
structure = CifData(file=os.path.abspath("./HKUST1.cif"))
structure_zeopp.label = "hkust1"
structure.label = "hkust1"

# Zeopp settings
zeopp_code = Code.get_from_string('zeo++@ovan')
zeopp_probe_radius_co2_trappe = Float(2.0) #(Angs) It will create 8 pore blocks for test purpose
zeopp_atomic_radii_file = SinglefileData(file=os.path.abspath("./UFF.rad")) # Radius file for the framework

# RASPA settings
raspa_code = Code.get_from_string('raspa2.0.37@ovan')

raspa_parameters = Dict(
    dict={
        "GeneralSettings": {
            "SimulationType": "MonteCarlo",
            "NumberOfCycles": 100,
            "NumberOfInitializationCycles": 200,
            "PrintEvery": 1000,
            "Forcefield": "GenericMOFs",
            "RemoveAtomNumberCodeFromLabel": True,
            "CutOff": 12.0,
        },
        "System": {
            "hkust1": {
                "type": "Framework",
                "UnitCells": "1 1 1",
                "HeliumVoidFraction": 0.149,
                "ExternalTemperature": 300.0,
                "ExternalPressure": 1e5,
            },
        },
        "Component": {
            "comp1": {
                "MoleculeDefinition": "TraPPE",
                "CreateNumberOfMolecules": 0,
            },
            "comp2": {
                "MoleculeDefinition": "TraPPE",
                "CreateNumberOfMolecules": 0,
            },
        },
    })

# raspa_options = {
#     "resources": {
#         "num_machines": 1,
#         "tot_num_cpus": 1,
#         "num_mpiprocs_per_machine": 1,
#         },
# }


# raspa_options = {
#     "resources": {
#         "num_machines": 1,
#         "tot_num_cpus": 1,
#         "num_mpiprocs_per_machine": 1,
#     },
    # "max_wallclock_seconds": 1 * 30 * 60,# 30 min
    # #"queue_name": "default",
    # "withmpi": False,
#}


# metadata = {
#
# }
#raspa_options = ParameterData(dict=raspa_options_dict)

# Running

# submit(SeparationWorkChain,
#     structure=structure,
#     raspa_code=raspa_code,
#     raspa_parameters=raspa_parameters,
#     raspa_options=raspa_options
#     )


submit(SeparationWorkChain,
    structure_zeopp=structure_zeopp,
    structure=structure,
    zeopp_code=zeopp_code,
    raspa_code=raspa_code,
    raspa_parameters=raspa_parameters,
    #raspa_options=raspa_options,
    #zeopp_options=zeopp_options,
    zeopp_probe_radius=zeopp_probe_radius_co2_trappe,
    zeopp_atomic_radii=zeopp_atomic_radii_file,
    #label='SeparationWorkChain-test-zeo1',
    )
