from __future__ import print_function
from aiida.orm.data.structure import StructureData
from aiida.orm.calculation.work import WorkCalculation
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg') #avoids ssh problems, but can not use .show()
import matplotlib.pyplot as plt
ParameterData = DataFactory('parameter')
StructureData = DataFactory('structure')
SinglefileData = DataFactory('singlefile')
CifData = DataFactory('cif')

# User settings
last = 0 #select the Nth last result in time
workflow1_label = 'test2-smearing'
workflow21_label = 'pe2-co2'
workflow22_label = 'pe2-n2'
with open('../cof_test2/list-smearing.list') as f:
    structure_labels=f.read().splitlines()
structure_labels=structure_labels[:50]
############################################################ Print vol & Kh info
ofile = open("parse_VolpoKhIsotherm.out","w+")

for structure_label in structure_labels:
    print("{} ".format(structure_label),end="",file=ofile)
    for workflow2_label in [workflow21_label,workflow22_label]:
        print("| ".format(structure_label),end="",file=ofile)
        q = QueryBuilder()
        q.append(StructureData, filters={'label': structure_label}, tag='inp_struct')
        q.append(WorkCalculation, filters={'label': workflow1_label}, output_of='inp_struct', tag='wf1')
        q.append(CifData, edge_filters={'label': 'output_structure'}, output_of='wf1', tag='cif')
        q.append(WorkCalculation, filters={'label': workflow2_label}, output_of='cif', tag='wf2')
        if len(q.all())>0:
            pk_work = q.all()[last][0].pk
            print("{} (PK: {} ) ".format(workflow2_label,pk_work),end="",file=ofile)
            q.append(ParameterData, output_of='wf2')
            # Parse results
            try:
                res = q.all()[last][0].get_dict()
            except: # the calculation is absent or not finished
                print("{:74s}".format(""),end="",file=ofile)
                res=None
            if res!=None:
                try:
                    print("POA-vf = {:.2f}, Blocking_spheres = {}, Kh = {:.2e} +/- {:.2e} {} ".format(
                        res['POAV_Volume_fraction'],
                        res['number_blocking_spheres'],
                        res['henry_coefficient_average'],
                        res['henry_coefficient_dev'],
                        res['henry_coefficient_units'],
                        ),
                        end="",file=ofile)
                except:
                    print("NOT POROUS {:62s} ".format(""),end="",file=ofile)
    print(end="\n",file=ofile)
ofile.close()
########################################## Print isotherms for parassitic energy

dir_out="./parse_VolpoKhIsotherm/"
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
for structure_label in structure_labels:
    fig, ax = plt.subplots(figsize=[8, 4.5],nrows=1, ncols=2, sharey=True)
    fig.suptitle('Isotherms for: '+ structure_label)
    ax[0].set(xlabel='Pressure (bar)',
              ylabel='Uptake (mol/kg)')
    ax[1].set(xlabel='Heat of adsorption (kJ/mol)',
              ylabel='Uptake (mol/kg)',)
    porous = False
    for i,workflow2_label in enumerate([workflow21_label,workflow22_label]):
        gas = ['CO_2','N_2'][i]

        q = QueryBuilder()
        q.append(StructureData, filters={'label': structure_label}, tag='inp_struct')
        q.append(WorkCalculation, filters={'label': workflow1_label}, output_of='inp_struct', tag='wf1')
        q.append(CifData, edge_filters={'label': 'output_structure'}, output_of='wf1', tag='cif')
        q.append(WorkCalculation, filters={'label':workflow2_label}, output_of='cif', tag='wf2')
        q.append(ParameterData, output_of='wf2')
        try:
            res = q.all()[last][0].get_dict()
        except: # the calculation is absent or not finished
            break
        if 'isotherm_loading' in res:
            porous = True
            # Parse data into arrays
            p = [a[0] for a in res['isotherm_loading']] #(bar)
            q_avg = [a[1] for a in res['isotherm_loading']] #(mol/kg)
            q_dev = [a[2] for a in res['isotherm_loading']] #(mol/kg)
            h_avg = [a[1] for a in res['isotherm_enthalpy']] #(kJ/mol)
            h_dev = [a[2] for a in res['isotherm_enthalpy']] #(kJ/mol)
            # TRICK: use the enthalpy from widom (energy-RT) which is more accurate that the one at 0.001 bar (and which also is NaN for weakly interacting systems)
            h_avg[0] = res['adsorption_energy_average']-res['temperature']/120.027
            h_dev[0] = res['adsorption_energy_dev']
            # Create directories
            structure_dir=dir_out+structure_label+"/"
            if not os.path.exists(structure_dir):
                os.makedirs(structure_dir)
            structure_dir_gas=structure_dir+gas+"/"
            if not os.path.exists(structure_dir_gas):
                os.makedirs(structure_dir_gas)
            # Print header and info
            ofile = open(structure_dir_gas+"300K.csv","w+")
            print("pressure(Pa) loading(mol/kg) HoA(kJ/mol)",file=ofile)

            for i in range(len(p)):
                print("{} {} {}".format(p[i]*1e5,q_avg[i],h_avg[i]),file=ofile)
            ofile.close()
            # Plot isotherm
            ax[0].errorbar(p,     q_avg, yerr=q_dev, marker ="o", label=gas)
            ax[1].errorbar(h_avg, q_avg, xerr=h_dev, marker ="o")
            # Print Density (converted from g/cm3 to kg/m3)
            ofile = open(structure_dir+"rho.csv","w+")
            print(res['Density']*1000,file=ofile)
            ofile.close()
        else: # the material is non-porous
            break


    if porous:
        ax[0].legend()
        fig.savefig(dir_out+"/"+structure_label+".png",dpi=300)

'''
To compute the parassitic energy:

cd parse_VolpoKhIsotherm
for f in */
  do calPE.py $f coal -datapath ./ > calPE.out
done
'''
