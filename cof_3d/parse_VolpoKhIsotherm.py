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
workflow21_label = 'isot-0_co2'
workflow22_label = 'isot-0_n2'
#with open('3dN.list') as f:
#    structure_labels=f.read().splitlines()
structure_labels = ['13180N']
############################################################ Print vol & Kh info
ofile = open("parse_VolpoKhIsotherm.out","w+")
for structure_label in structure_labels:
    print("{} | ".format(structure_label),end="",file=ofile)
    for workflow2_label in [workflow21_label,workflow22_label]:
        print("{} ".format(workflow2_label),end="",file=ofile)
        q = QueryBuilder()
        q.append(StructureData, filters={'label': structure_label}, tag='inp_struct')
        q.append(WorkCalculation, filters={'label': workflow1_label}, output_of='inp_struct', tag='wf1')
        q.append(CifData, edge_filters={'label': 'output_structure'}, output_of='wf1', tag='cif')
        q.append(WorkCalculation, filters={'label': workflow2_label}, output_of='cif', tag='wf2')
        q.append(ParameterData, output_of='wf2')
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
ofile.close()
########################################## Print isotherms for parassitic energy
dir_out="./parse_VolpoKhIsotherm/"
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
for structure_label in structure_labels:
    structure_dir=dir_out+structure_label+"/"
    if not os.path.exists(structure_dir):
        os.makedirs(structure_dir)
    for i,workflow2_label in enumerate([workflow21_label,workflow22_label]):
        gas = ['CO_2','N_2'][i]
        structure_dir_gas=structure_dir+gas+"/"
        if not os.path.exists(structure_dir_gas):
            os.makedirs(structure_dir_gas)
        q = QueryBuilder()
        q.append(StructureData, filters={'label': structure_label}, tag='inp_struct')
        q.append(WorkCalculation, filters={'label': workflow1_label}, output_of='inp_struct', tag='wf1')
        q.append(CifData, edge_filters={'label': 'output_structure'}, output_of='wf1', tag='cif')
        q.append(WorkCalculation, filters={'label':workflow2_label}, output_of='cif', tag='wf2')
        q.append(ParameterData, output_of='wf2')
        res = q.all()[0][0].get_dict()
        ofile = open(structure_dir_gas+"300.csv","w+")
        print("pressure(Pa) loading(mol/kg) HoA(kJ/mol)",file=ofile)
        for i in range(len(res['isotherm_loading'])):
            try:
                p = res['isotherm_loading'][i][0]*1e5 #(Pa)
                q = res['isotherm_loading'][i][1] #(mol/kg)
                h = res['isotherm_enthalpy'][i][1] #(kJ/mol)
                print("{} {} {}".format(p,q,h),file=ofile)
            except:
                pass #NONPOROUS
        ofile.close()
