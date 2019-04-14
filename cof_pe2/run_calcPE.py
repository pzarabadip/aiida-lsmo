from __future__ import print_function
from aiida.orm.calculation.work import WorkCalculation
from aiida.orm.data.base import Float, Str
from aiida.work import workfunction as wf
import sys
import os
import numpy as np
import pandas as pd
from calc_pe.utils import printPE

ParameterData = DataFactory('parameter')
StructureData = DataFactory('structure')
SinglefileData = DataFactory('singlefile')
CifData = DataFactory('cif')


@wf
def CalcPE(ResGasCO2, ResGasN2, GasIn, VF, Process, Cp, Yd, ElEff, Opt):
    """ Submit calc_pe calculation using AiiDA, for the CO2 parasitic energy.
    :ParameterData ResGasCO2:
    :ParameterData ResGasN2:
    :Str GasIn:
    :Float VF:
    :Str Process:
    :Float Cp:
    :Float Yd:
    :Str ElEff:
    :Str Opt:
    """
    from calc_pe import mainPE
    T_iso = {}
    iso_df = {}
    for i, ResGas in enumerate([ResGasCO2, ResGasN2]):
        gas = ['CO_2','N_2'][i]
        res = ResGas.get_dict()
        T_iso[gas] = res["temperature"]
        iso_df[gas] = pd.DataFrame(columns = ['pressure(Pa)','loading(mol/kg)','HoA(kJ/mol)'])
        iso_df[gas]['pressure(Pa)'] = [a[0]*1e5 for a in res['isotherm_loading']] # converted from bar to Pa
        iso_df[gas]['loading(mol/kg)'] = [a[1] for a in res['isotherm_loading']]
        iso_df[gas]['HoA(kJ/mol)'] = [a[1] for a in res['isotherm_enthalpy']]
        # TRICK: use the enthalpy from widom (energy-RT) which is more accurate
        #        that the one at 0.001 bar (and which also is NaN for weakly interacting systems)
        iso_df[gas]['HoA(kJ/mol)'].loc[0] = res['adsorption_energy_average']-res['temperature']/120.027

    pe_dict = mainPE(gasin = GasIn.value,        # e.g., "coal"
                      rho = res['Density']*1000, # converted from g/cm3 to kg/m3
                      vf = VF.value,             # e.g., 0.35
                      process = Process.value,   # e.g., "TPSA"
                      cp = Cp.value,             # e.g., 985.0, the average for MOFs in J/kg/K
                      yd = Yd.value,             # e.g., 0.99
                      eleff = ElEff.value,       # e.g., "carnot"
                      opt = Opt.value,           # e.g., "PE"
                      T_iso = T_iso,             # [T_iso_CO2, T_iso_N2]
                      iso_df = iso_df            # [iso_df_CO2, iso_df_N2] df containing pressure (Pa), loading (mol/kg), HeatOfAdsorption (kJ/mol)
                     )
    return ParameterData(dict=pe_dict)

# User settings
last = 0 #select the Nth last result in time
workflow1_label = 'test2-smearing'
workflow21_label = 'pe2-co2'
workflow22_label = 'pe2-n2'
with open('../cof_test2/list-smearing.list') as f:
    structure_labels=f.read().splitlines()

for structure_label in structure_labels[-1]:
    porous = False
    resgas = {}
    for i,workflow2_label in enumerate([workflow21_label,workflow22_label]):
        gas = ['CO_2','N_2'][i]
        q = QueryBuilder()
        q.append(StructureData, filters={'label': structure_label}, tag='inp_struct')
        q.append(WorkCalculation, filters={'label': workflow1_label}, output_of='inp_struct', tag='wf1')
        q.append(CifData, edge_filters={'label': 'output_structure'}, output_of='wf1', tag='cif')
        q.append(WorkCalculation, filters={'label':workflow2_label}, output_of='cif', tag='wf2')
        q.append(ParameterData, output_of='wf2')
        if len(q.all())>0:
            resgas[gas] = q.all()[last][0]
            # print(resgas[gas].pk) #debug
    # Check if the structure is present, porous (has 'isotherm_loading') and compute PE
    if ('CO_2' in resgas.keys() and
        'N_2' in resgas.keys() and
        'isotherm_loading' in resgas['CO_2'].get_dict() and
        'isotherm_loading' in resgas['N_2'].get_dict()):
        porous = True
        ResultPE = CalcPE(ResGasCO2=resgas['CO_2'],
                          ResGasN2=resgas['N_2'],
                          GasIn=Str('coal'),
                          VF=Float(0.35),
                          Process=Str('TPSA'),
                          Cp=Float(985.0),
                          Yd=Float(0.99),
                          ElEff=Str('carnot'),
                          Opt=Str('PE'),
                         )
        printPE(structure_label,ResultPE.get_dict())
    if not porous:
        print("{}: not found or not porous".format(structure_label))
