from __future__ import print_function
from aiida.orm.data.structure import StructureData
from aiida.orm.calculation.work import WorkCalculation
import sys
import os
import matplotlib
matplotlib.use('Agg') #avoids ssh problems, but can not use .show()
import matplotlib.pyplot as plt
import numpy as np
sys.path.append('/home/daniele/Programs/cp2k_utils/parse_steps')
from parser_utils import print_header, print_steps

# User settings
last = 0 #select the Nth last result in time
workflow1_label = '3DCOFs-600K-OptAngles'
workflow2_label = 'volpo-Kh-CO2-test1'
with open('3dN.list') as f:
    structure_labels=f.read().splitlines()

for structure in structure_labels:
    ofile = open("parse_volpoKh.out","w+")
        #structure = "09000N"
        qb = QueryBuilder()
        qb.append(StructureData, filters={'label': structure}}, tag='structure')
        qb.append(WorkCalculation, filters={'label': workflow1_label}, output_of='structure', tag='wf1')
        qb.append(WorkCalculation,  edge_filters={'label': workflow2_label}, output_of='wf1', tag='wf2')
        if istage == 0:
            qb.append(WorkCalculation, filters={'label':{'==':stage_name[0]}}, output_of='robustgeoopt', tag='dftbase')
        else:
            qb.append(WorkCalculation, filters={'label':{'==':stage_name[istage]}}, output_of='robustgeoopt', tag='stage')
            qb.append(WorkCalculation, filters={'label':{'==':'Cp2kDftBaseWorkChain'}}, output_of='stage', tag='dftbase')
        qb.append(JobCalculation,output_of='dftbase',tag='calc')
        qb.order_by({WorkCalculation:{'ctime':'desc'}})
        try:
            stage_localpath[istage] =  qb.all()[last][0].out.retrieved.get_abs_path()+'/path/aiida.out'
            stage_pk[istage] = str(qb.all()[last][0].pk)
        except:
            pass
        # Print file with pk and local directory
        print('%s\t%s\t%s' %(stage_name[istage],stage_pk[istage],stage_localpath[istage]), file=pkfile)

    pkfile.close()
    # Print on screen the state of the calculation
    completed_stage = 4 - stage_localpath.count("INCOMPLETE")
    if completed_stage == -1:
        mex = "None"
    else:
        mex = stage_name[completed_stage]
    print('%s\tCompleted: %s' %(structure,mex))

    file_out = dir_out + structure + "_" + workflow_label + "_steps.out"
    with open(file_out, 'w+') as fout:
        with redirect_stdout(fout):
            print_header()
            for i in range(completed_stage+1):
                print_steps(stage_localpath[i])
    if completed_stage>0: # at least Stage1 completed
        plot_steps(file_out, structure)
