from __future__ import print_function
import pandas as pd
import numpy as np
from aiida_qeq.calculations.qeq import QeqCalculation
from aiida.orm.calculation.work import WorkCalculation

CifData = DataFactory('cif')
ParameterData = DataFactory('parameter')

df = pd.read_csv("../cof_test2/pk_final.csv")
qeqcalc_label = "aiida_qeq Qeq with GMP parameters"
khwork_label = 'qeq-kh-co2'

# Initialize file and print header
khfile = open("./parse_co2kh-qeq.csv", 'w+')
print('cof_label,kH-qeq (mol/kg/Pa),kH-qeq dev (mol/kg/Pa)', file=khfile)

for i in df.index:
    cif_label = df.at[i,'structure']
    dict_khworkres = {'henry_coefficient_average': np.nan,
                      'henry_coefficient_dev': np.nan}
    #Parse and print the qeq-CifData.pk
    q = QueryBuilder()
    q.append(CifData, filters={'label': cif_label}, tag='cif-inp')
    q.append(QeqCalculation, filters={'label': qeqcalc_label},
                             output_of='cif-inp', tag='qeq')
    if len(q.all())==0:
        print("CifData.label= {} failed QeqCalculation missing.".format(
              cif_label,pk_qeqcalc))
    else:
        pk_qeqcalc=q.all()[0][0].pk
        q.append(CifData, output_of='qeq', tag='cif-qeq')
        if len(q.all())==0:
            print("CifData.label= {} failed QeqCalculation.pk= {}".format(
                  cif_label,pk_qeqcalc))
        else:
            # Parse and print Kh values
            pk_cifqeq=q.all()[0][0].pk
            q.append(WorkCalculation, filters={'label': khwork_label},
                                      output_of='cif-qeq', tag='khwork')
            try:
                pk_workcalc =  q.all()[0][0].pk
            except:
                pk_workcalc = 'none'
            q.append(ParameterData, edge_filters={'label': 'results'},
                                    output_of='khwork', tag='khwork_results')
            if len(q.all())==0:
                print("CifData.label= {} with Qeq >>> CifData.pk= {} >>> Failed Kh-calc".format(
                      cif_label,pk_cifqeq))
            else:
                khworkres = q.all()[0][0]
                pk_khworkres = khworkres.pk
                dict_khworkres = khworkres.get_dict()
                if 'henry_coefficient_average' not in dict_khworkres.keys(): #not porous
                     dict_khworkres['henry_coefficient_average']=0
                     dict_khworkres['henry_coefficient_dev']=0
                print("CifData.label= {} with Qeq >>> CifData.pk= {} >>> Kh= {} +/- {} mol/kg/Pa (WorkCalculation.pk = {} ParameterData.pk= {} )".format(
                      cif_label,
                      pk_cifqeq,
                      dict_khworkres['henry_coefficient_average'],
                      dict_khworkres['henry_coefficient_dev'],
                      pk_workcalc,
                      pk_khworkres,
                      ))
    print('{},{},{}'.format(
          cif_label,
          dict_khworkres['henry_coefficient_average'],
          dict_khworkres['henry_coefficient_dev'],
          ), file=khfile)
khfile.close()
