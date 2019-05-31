from __future__ import print_function
from aiida.orm.data.structure import StructureData
from aiida.orm.calculation.work import WorkCalculation
import sys
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg') #avoids ssh problems, but can not use .show()
import matplotlib.pyplot as plt
ParameterData = DataFactory('parameter')
StructureData = DataFactory('structure')
SinglefileData = DataFactory('singlefile')
CifData = DataFactory('cif')

df = pd.read_csv("../cof_test2/pk_final.csv")
############################################################ Print vol & Kh info, and plot isotherms
# Print file
ofile = open("parse_VolpoKhIsotherm.csv","w+")
#plot figure
pmax=30
plotcolor={'co2':'red','n2':'blue'}
plotlabel={'co2':'CO$_2$','n2':'N$_2$'}
dir_out="./parse_VolpoKhIsotherm/"
if not os.path.exists(dir_out):
    os.makedirs(dir_out)

#header
print("cof_label,", end="", file=ofile)
print("co2_pe_label,co2_pe_pk,co2_poavf,co2_blocksph,co2_kh (mol/kg/Pa),co2_dev,", end="", file=ofile)
print("n2_pe_label,n2_pe_pk,n2_poavf,n2_blocksph,n2_kh (mol/kg/Pa),n2_dev", end="\n", file=ofile)

for i in df.index:
    res_found={'co2':False,'n2':False}
    structure_label = df.at[i,'structure']
    #if structure_label in ['05001N2','13180N3','18041N3']:
    if True:
        # print stuff
        print("Querying: {}".format(structure_label), end=" ") # on screen
        print(structure_label,end=",",file=ofile) # on file
        # initialize figure
        fig, ax = plt.subplots(figsize=[7, 4.5],
                               nrows=1,
                               ncols=2,
                               sharey=True,
                               gridspec_kw={"width_ratios":[1.2,0.8],
                                    "wspace":0.0})
        fig.suptitle('Isotherms for: '+ structure_label)
        ax[0].set(xlabel='Pressure (bar)',
                  ylabel='Uptake (mol/kg)')
        ax[1].set(xlabel='Heat of adsorption (kJ/mol)')
        porous = False

        optcif_pk = df.at[i,'CifData']
        #if optcif_pk == 'none': #not to make the QB crash
        #    optcif_pk = 404
        for gas in ['co2','n2']:
            # initialize variables
            pe_pk = 404
            pe_res = None
            poavf=np.nan
            bs=np.nan
            kh=np.nan
            khdev=np.nan
            # query
            pe_label = df.at[i,'pe.label'].replace('xx',gas)
            q = QueryBuilder()
            q.append(CifData, filters={'id': optcif_pk}, tag='optcif')
            q.append(WorkCalculation, filters={'label': pe_label}, output_of='optcif', tag='pewc')
            if len(q.all())>0:  #existing pe workchain
                pe_pk = q.all()[0][0].pk
                q.append(ParameterData, edge_filters={'label': 'results'},
                                        output_of='pewc', tag='pe_res')
                # Parse results
                try:
                    pe_res = q.all()[0][0].get_dict() # this will fail if (1) len(q.all)=0 or (2) last pe WC is not finished yet
                except:
                    pass
                if pe_res!=None: #existing results
                    res_found[gas] = True
                    poavf = pe_res['POAV_Volume_fraction']
                    if 'henry_coefficient_average' in pe_res.keys(): #porous
                        # parse info for printing
                        bs = pe_res['number_blocking_spheres']
                        kh = pe_res['henry_coefficient_average']
                        khdev = pe_res['henry_coefficient_dev']
                        # parse isotherm for plotting
                        p = [a[0] for a in pe_res['isotherm_loading']] #(bar)
                        q_avg = [a[1] for a in pe_res['isotherm_loading']] #(mol/kg)
                        q_dev = [a[2] for a in pe_res['isotherm_loading']] #(mol/kg)
                        h_avg = [a[1] for a in pe_res['isotherm_enthalpy']] #(kJ/mol)
                        h_dev = [a[2] for a in pe_res['isotherm_enthalpy']] #(kJ/mol)
                        # TRICK: use the enthalpy from widom (energy-RT) which is more accurate that the one at 0.001 bar (and which also is NaN for weakly interacting systems)
                        h_avg[0] = pe_res['adsorption_energy_average']-pe_res['temperature']/120.027
                        h_dev[0] = pe_res['adsorption_energy_dev']
                    else: #non porous
                        kh = 0.0
                        khdev = 0.0
                        p = [0, pmax]
                        q_avg = [0, 0]
                        q_dev = [0, 0]
                        h_avg = [0, 0]
                        h_dev = [0, 0]
                    # Plot isotherm
                    ax[0].errorbar(p, q_avg, yerr=q_dev,
                                   marker ="o",
                                   color= plotcolor[gas],
                                   label=plotlabel[gas])
                    ax[1].errorbar(h_avg, q_avg, xerr=h_dev,
                                   marker ="o",
                                   color= plotcolor[gas],)
            print(*[pe_label,pe_pk,poavf,bs,kh,khdev], sep=",", end="", file=ofile)
            if gas=='co2':
                print(end=",",file=ofile)
        # Both gas parsed: (1) print \n to .csv, (2) plot isotherms, (3) print status
        print(end="\n", file=ofile)
        if res_found['co2'] or res_found['n2']:
            ax[0].legend(loc='upper left')
            ax[0].set_xlim([-1,pmax+1])
            ax[0].set_ylim([-1,41])
            ax[1].set_xlim([-51,0+1])
            fig.savefig(dir_out+"/"+structure_label+".png",dpi=300,bbox_inches='tight')
        if res_found['co2'] and res_found['n2']:
            print("... CO2 and N2 parsed correctly!")
        else:
            print("... someting FAILED or STILL RUNNING.")
ofile.close()
