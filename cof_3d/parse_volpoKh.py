from __future__ import print_function
from aiida.orm.data.structure import StructureData
from aiida.orm.calculation.work import WorkCalculation
import sys
import os
import numpy as np
ParameterData = DataFactory('parameter')
StructureData = DataFactory('structure')
SinglefileData = DataFactory('singlefile')
CifData = DataFactory('cif')

# User settings
last = 0 #select the Nth last result in time
workflow1_label = '3DCOFs-600K-OptAngles'
workflow21_label = 'volpo-Kh-CO2-test1'
workflow22_label = 'volpo-Kh-N2-test1'
with open('3dN.list') as f:
    structure_labels=f.read().splitlines()

ofile = open("parse_volpoKh.out","w+")
for structure_label in structure_labels:
    print("{} | ".format(structure_label),end="",file=ofile)
    for workflow2_label in [workflow21_label,workflow22_label]:
        print("{} ".format(workflow2_label),end="",file=ofile)
        q = QueryBuilder()
        q.append(StructureData, filters={'label': structure_label}, tag='inp_struct')
        q.append(WorkCalculation, filters={'label': workflow1_label}, output_of='inp_struct', tag='wf')
        q.append(CifData, edge_filters={'label': 'output_structure'}, output_of='wf', tag='cif')
        q.append(WorkCalculation, filters={'label':workflow2_label}, output_of='cif', tag='kh')
        q.append(ParameterData, output_of='kh')
        res = q.all()[0][0].get_dict()
        print("POA-vf = {:.2f}, ".format(res['POAV_Volume_fraction']),
              end="",file=ofile)
        try:
          print("Blocking_spheres = {}, Kh = {:.2e} +/- {:.2e} {} | ".format(
                res['number_blocking_spheres'],
                res['henry_coefficient_average'],
                res['henry_coefficient_dev'],
                res['henry_coefficient_units'],
                ),
                end="",file=ofile)
        except:
            print("NOT POROUS {:47s} | ".format(""),end="",file=ofile)
    print(file=ofile)
