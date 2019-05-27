from __future__ import absolute_import
from __future__ import print_function
import os
import pandas as pd
from aiida.common.example_helpers import test_and_get_code
from aiida.orm import DataFactory, CalculationFactory
from aiida.work import workfunction as wf
ParameterData = DataFactory('parameter')
StructureData = DataFactory('structure')
CifData = DataFactory('cif')
ZeoppCalculation = CalculationFactory('zeopp.network')
NetworkParameters = DataFactory('zeopp.parameters')

global zeopp_prop

def zeopp_list(i,notes,qall):
    """ Make a string out of a dict, in the following style:
    cofid,pk,note,properties
    """
    global zeopp_prop
    if len(qall)==0:
        list = [df.at[i,"structure"], "nopk", "ERROR_missing_output"]
    else:
        d = qall[0][0].get_dict()
        list = [df.at[i,"structure"],qall[0][0].pk,notes]
        for entry in zeopp_prop:
            try:
                list.append(d[entry])
            except:
                list.append("parse_fail")
    return list

# Declare the interested zeopp's properties, then build and print the header
zeopp_prop= [
            "Density",
            "Unitcell_volume",
            "ASA_m^2/cm^3",
            "NASA_m^2/cm^3",
            "AV_cm^3/g",
            "AV_Volume_fraction",
            "POAV_cm^3/g",
            "PONAV_cm^3/g",
            "Number_of_pockets",
            "Number_of_channels",
            "Largest_free_sphere",
            "Largest_included_sphere",
            ]
header_list=["structure","pk_ParDat","notes"]+zeopp_prop
print(*header_list,sep=",")

# Parse the restulst starting from a previous csv
df=pd.read_csv("./pk_final.csv")
for i in range(df.shape[0]):

        # parse Initial
        q = QueryBuilder()
        q.append(StructureData,filters={'id': df.at[i,"StructureData"]})
        q.append(CifData,descendant_of=StructureData)
        q.append(JobCalculation,output_of=CifData)
        q.append(ParameterData,output_of=JobCalculation)
        print(*zeopp_list(i,"unrelaxed",q.all()),sep=",")

        # Parse Optimized
        if df.at[i,"CifData"] != "none":
            q = QueryBuilder()
            q.append(CifData,filters={'id': df.at[i,"CifData"]})
            q.append(JobCalculation,output_of=CifData)
            q.append(ParameterData,output_of=JobCalculation)
            print(*zeopp_list(i,"cp2k_opt",q.all()),sep=",")
        else:
            print(*[df.at[i,"structure"], "nopk", "ERROR_missing_opt_structure"],sep=",")
