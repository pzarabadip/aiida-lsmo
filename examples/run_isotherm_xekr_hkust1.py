import os
from aiida.plugins import DataFactory
from aiida.orm import Code, Dict, Float, Int, Str
from aiida.engine import run, submit
from aiida_lsmo_workflows.isotherm_multi_comp import MultiCompIsothermWorkChain

ParameterData = DataFactory('dict')
SinglefileData = DataFactory('singlefile')
CifData = DataFactory('cif')

# Import the structure
structure = CifData(file=os.path.abspath("./HKUST1.cif"))
structure.label = structure.filename.lower()[:-4]

# Zeopp settings
zeopp_code = Code.get_from_string('zeopp_teslin@teslin')
zeopp_probe_radius_co2_trappe = Float(2.0) #(Angs) It will create 8 pore blocks for test purpose
zeopp_atomic_radii_file = SinglefileData(file=os.path.abspath("./UFF.rad")) # Radius file for the framework

# RASPA settings
raspa_code = Code.get_from_string('raspa_teslin@teslin')

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

raspa_comp = {

    'comp1':{'name':'xenon','mol_fraction':0.2,'radius':1.985, 'mol_def':'TraPPE'},
    'comp2':{'name':'krypton','mol_fraction':0.8,'radius':1.986, 'mol_def':'TraPPE'}

    }


submit(MultiCompIsothermWorkChain,
    structure=structure,
    zeopp_code=zeopp_code,
    raspa_code=raspa_code,
    raspa_parameters=raspa_parameters,
    raspa_isotherm_dynamic=False,
    raspa_isotherm_full=False,
    raspa_comp = raspa_comp,
    selected_pressures=[0.15e5,0.55e5,1.05e5],
    # zeopp_probe_radius=zeopp_probe_radius_co2_trappe,
    zeopp_atomic_radii=zeopp_atomic_radii_file,
    zeopp_accuracy=Str('DEF'),
    # label='SeparationWorkChain',
    )
