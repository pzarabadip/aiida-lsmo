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
        iso_df[gas]['pressure(Pa)'] = [a[0]*1e5 for a in res['isotherm_loading']] #(bar>>Pa)
        iso_df[gas]['loading(mol/kg)'] = [a[1] for a in res['isotherm_loading']] #(mol/kg)
        iso_df[gas]['HoA(kJ/mol)'] = [a[1] for a in res['isotherm_enthalpy']] #(kJ/mol)
        # TRICK: use the enthalpy from widom (energy-RT) which is more accurate
        #        that the one at 0.001 bar (and which also is NaN for weakly interacting systems)
        iso_df[gas]['HoA(kJ/mol)'].loc[0] = res['adsorption_energy_average']-res['temperature']/120.027
    pe_dict = mainPE(gasin = GasIn.value,
                      rho = res['Density'],
                      vf = VF.value,
                      process = Process.value,
                      cp = Cp.value,
                      yd = Yd.value,
                      eleff = ElEff.value,
                      opt = Opt.value,
                      T_iso = T_iso,  # [T_iso_CO2, T_iso_N2]
                      iso_df = iso_df # [iso_df_CO2, iso_df_N2]
                     )

    return ParameterData(dict=pe_dict)

# User settings
last = 0 #select the Nth last result in time
workflow1_label = 'test2-0'
workflow21_label = 'pe2-co2'
workflow22_label = 'pe2-n2'
with open('../cof_test2/list-OT.list') as f:
    structure_labels=f.read().splitlines()
# test
workflow1_label = '3DCOFs-600K-OptAngles'
workflow21_label = 'isot-2_co2'
workflow22_label = 'isot-2_n2'
structure_labels=['18122N']

for structure_label in structure_labels:
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
        try:
            resgas[gas] = q.all()[last][0]
        except: # the calculation is absent, not finished or the structure is non-porous
            break
    # Check if the structure is porous (has 'isotherm_loading') and compute PE
    if ('isotherm_loading' in resgas['CO_2'].get_dict()) and ('isotherm_loading' in resgas['N_2'].get_dict()):
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
