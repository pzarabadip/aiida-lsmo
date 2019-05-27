from __future__ import print_function
import pandas as pd
import numpy as np
from aiida_qeq.calculations.qeq import QeqCalculation
from aiida.orm.calculation.work import WorkCalculation

CifData = DataFactory('cif')
ParameterData = DataFactory('parameter')

df = pd.read_csv("../cof_test2/pk_final.csv")

# Initialize file and print header
khfile = open("./parse_compare_khco2.csv", 'w+')
print("cof_label", end=",", file=khfile)
print("kH-ddec (mol/kg/Pa),kH-ddec_dev,adsE-ddec (kJ/mol),adsE-ddec_dev", end=",", file=khfile)
print("kH-qeq (mol/kg/Pa),kH-qeq_dev,adsE-qeq (kJ/mol),adsE-qeq_dev", end="\n", file=khfile)

for i in df.index:
    cif_label = df.at[i,'structure']
    print("Querying {}".format(cif_label))
    print(cif_label,end=",", file=khfile)
    for j in range(2): # select first ddec and then qeq
        volpokhiso_pk = df.at[i,['pex-co2.pk','qeq-kh-co2.pk'][j]]
        if volpokhiso_pk == 404:
            widom_res_dict = {'henry_coefficient_average': np.nan,
                              'henry_coefficient_dev': np.nan,
                              'adsorption_energy_widom_average': np.nan,
                              'adsorption_energy_widom_dev': np.nan}
        else:
            q = QueryBuilder()
            q.append(WorkCalculation, filters={'id': volpokhiso_pk}, tag='volpokhiso')
            q.append(WorkCalculation, filters={'label': 'RaspaWidom'},
                                      output_of='volpokhiso', tag='widomcalc')
            q.append(ParameterData, edge_filters={'label': 'component_0'},
                                    output_of='widomcalc', tag='widom_res')
            if len(q.all())==0: # not present because not porous
                widom_res_dict = {'henry_coefficient_average': 0,
                                  'henry_coefficient_dev': 0,
                                  'adsorption_energy_widom_average': 0,
                                  'adsorption_energy_widom_dev': 0}
            else:
                widom_res = q.all()[0][0]
                widom_res_pk = widom_res.pk
                widom_res_dict = widom_res.get_dict()

        print('{},{},{},{}'.format(
                  widom_res_dict['henry_coefficient_average'],
                  widom_res_dict['henry_coefficient_dev'],
                  widom_res_dict['adsorption_energy_widom_average'],
                  widom_res_dict['adsorption_energy_widom_dev']
                  ), end="", file=khfile)
        if j==0:
            print(end=",", file=khfile)
        elif j==1:
            print(end="\n", file=khfile)
khfile.close()
