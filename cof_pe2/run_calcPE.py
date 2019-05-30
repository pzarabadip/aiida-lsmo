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

df = pd.read_csv("../cof_test2/pk_final.csv")

for i in df.index:
    if (df.at[i,'extension_ok']==0) and (df.at[i,'note']=='still_running'):
        structure_label = df.at[i,'structure']
        porous = False
        resgas = {}
        for gas in ['co2','n2']:
            pe_pk = df.at[i,'pex-xx.pk'.replace('xx',gas)]
            q = QueryBuilder()
            q.append(WorkCalculation, filters={'id':pe_pk}, tag='pe_workcalc')
            q.append(ParameterData, output_of='pe_workcalc')
            if len(q.all())>0:
                resgas[gas] = q.all()[0][0]
            # Check if the structure is present, porous (has 'isotherm_loading') and compute PE
        if ('co2' in resgas.keys() and
            'n2' in resgas.keys() and
            'isotherm_loading' in resgas['co2'].get_dict() and
            'isotherm_loading' in resgas['n2'].get_dict()):
            porous = True
            ResultPE = CalcPE(ResGasCO2=resgas['co2'],
                              ResGasN2=resgas['n2'],
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
