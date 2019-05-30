from __future__ import print_function
from aiida.orm.calculation.work import WorkCalculation
import sys
import os
import numpy as np
import pandas as pd
from calc_pe.utils import printPE

ParameterData = DataFactory('parameter')
StructureData = DataFactory('structure')
SinglefileData = DataFactory('singlefile')
CifData = DataFactory('cif')

df = pd.read_csv("../cof_test2/pk_final.csv")

debug = False # print pk and if non-porous

for i in df.index:
    structure_label = df.at[i,'structure']
    pe_pk = df.at[i,'pex-co2.pk']
    PE_out_pk = np.nan
    q = QueryBuilder()
    q.append(WorkCalculation, filters={'id':pe_pk}, tag='pe_workcalc')
    q.append(ParameterData, output_of='pe_workcalc', tag='pe_parameterdata')
    q.append(WorkCalculation, output_of='pe_parameterdata', tag='calcpe_workcalc')
    q.append(ParameterData, output_of='calcpe_workcalc', tag='calcpe_parameterdata')
    if len(q.all())>0:
        node = q.all()[0][0]
        PE_out_pk = node.pk
        PE_out_dict = node.get_dict()
        if debug:
            printPE("{} (pk={})".format(structure_label, PE_out_pk),PE_out_dict)
        else:
            printPE(structure_label, PE_out_dict)
    else:
        if debug:
            print("{} (pk={})".format(structure_label, PE_out_pk))
