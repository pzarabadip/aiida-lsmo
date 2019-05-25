# -*- coding: utf-8 -*-
"""Submit a test calculation on localhost.

Usage: verdi run submit.py

Note: This script assumes you have set up computer and code as in README.md.
"""
from __future__ import absolute_import
from __future__ import print_function
import os
import pandas as pd
from aiida.common.example_helpers import test_and_get_code
import aiida_qeq.data.qeq as data
from aiida_qeq.data import DATA_DIR
from aiida.orm import DataFactory

# make sure the "eqeq" binary is in your PATH
code = test_and_get_code('egulp@deneb', expected_code_type='qeq.qeq')

# Prepare input parameters
CifData = DataFactory('cif')
SinglefileData = DataFactory('singlefile')

parameter_file = SinglefileData(
    file=os.path.join(DATA_DIR, data.DEFAULT_PARAM_FILE_NAME))

df=pd.read_csv("../cof_test2/pk_final.csv")

for i in df.index:
    cif_label = df.at[i,'structure']
    if cif_label in ['12061N2', '12062N2', '17060N2', '17061N2', '18070N2', '18121N3', '18122N3']: #failedbecause of bad parsing
        cif = CifData(
            file='/home/daniele/Documents/CoRE-COFs/cifs/{}.cif'.format(cif_label),
            parse_policy='lazy')
        cif.label= cif_label

        # set up calculation
        calc = code.new_calc()
        calc.label = "aiida_qeq Qeq with GMP parameters"
        calc.set_max_wallclock_seconds(3 * 60 * 60)
        calc.set_withmpi(False)
        calc.set_resources({"num_machines": 1, "num_mpiprocs_per_machine": 1})
        calc.use_parameters(parameter_file)
        calc.use_structure(cif)
        calc.store_all()

        # submit calculation and print pk
        calc.submit()
        print("CifData.label= {} (pk= {}), submitted calc.pk= {}".format(
              cif.label, cif.pk, calc.dbnode.pk))
