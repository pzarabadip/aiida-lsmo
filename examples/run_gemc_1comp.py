from __future__ import print_function

import os
import sys

from aiida.common import NotExistent
from aiida.plugins import DataFactory
from aiida.orm import Code, Dict
from aiida.engine import submit
from aiida_lsmo_workflows.gemc import GEMCWorkChain

ParameterData = DataFactory("dict")

# Reading code information from system argv
if len(sys.argv) != 2:
    print("Usage: test.py <zeopp_code_name> <raspa_code_name>")
    sys.exit(1)

raspa_codename = sys.argv[1]

try:
    raspa_code = Code.get_from_string(raspa_codename)
except NotExistent:
    print("The code '{}' does not exist".format(raspa_codename))
    sys.exit(1)

general_calc_params = Dict(dict={
    "raspa" :{
                "gibbs_volume_prop" : 0.1,
                "verbosity" : 10,
                "cutoff" : 12.0,
                "usecharges": False,
                "simulation_type": "MonteCarlo",
                "system_type" : "Box",
                "additional_cycle" : 500,
                "T_min" : 100,
                "T_max" : 180,
                "dT" : 40,
    },
}
)

raspa_comp = {
    "comp1":{"name":"methane","box_one":150, "box_two":150, "mol_def":"TraPPE", "conv_threshold" : 0.10, "singlebead":True},
    }

raspa_box = {
    "box1" :{"tag":"box_one","box_ax":24,"box_by":24,"box_cz":24,
            "box_alpha":90,"box_beta":90, "box_gamma":90,"temperature":300.0},
    "box2" :{"tag":"box_two","box_ax":24,"box_by":24,"box_cz":24,
            "box_alpha":90,"box_beta":90, "box_gamma":90,"temperature":200.0},
}

raspa_parameters = Dict(
    dict={
        "GeneralSettings": {
            "NumberOfCycles": 100,
            "NumberOfInitializationCycles": 100,
            "PrintEvery": 1000,
            "Forcefield": "GenericMOFs",
        },
        "System": {"box"+str(i+1):{} for i in range(len(list(raspa_box)))},
        "Component": {"comp"+str(i+1):{} for i in range(len(list(raspa_comp)))}
        },
    )

raspa_options = {
    "resources": {
        "num_machines": 1,
        "tot_num_mpiprocs": 1,
    },
    "max_memory_kb": 200000,
    "max_wallclock_seconds": 2 * 60 * 60,
    "withmpi": False,
}

submit(GEMCWorkChain,
    raspa_code=raspa_code,
    raspa_parameters=raspa_parameters,
    raspa_comp = raspa_comp,
    raspa_box = raspa_box,
    general_calc_params=general_calc_params,
    raspa_options=raspa_options,
    metadata={
        "label" : "MultiCompIsothermWorkChain",
        "description" : "Test for GEMC "}
    )
