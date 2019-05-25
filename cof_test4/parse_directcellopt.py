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

from contextlib import contextmanager
@contextmanager
def redirect_stdout(target):
    original = sys.stdout
    sys.stdout = target
    yield
    sys.stdout = original

# User settings
last = 0 #select the Nth last result in time
#workflow_label = 'directcellopt'
workflow_label = 'test4fidis-OT'
with open('list-directcellopt.list') as f:
    structure_labels=f.read().splitlines()
# General settings
dir_out="./parse_steps/"
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
dir_out="./parse_steps/"+workflow_label+"/"
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
# Search for the path for the 4 steps: energy, md, geo_opt, cell_opt. Store them.
for structure in structure_labels:
    # get first the pk of the upper workflow
    localpath=None
    qb = QueryBuilder()
    qb.append(StructureData, filters={'label': {'==':structure}}, tag='structure')
    qb.append(WorkCalculation, filters={'label':{'==':workflow_label}}, output_of='structure', tag='workflow')
    pk_work=str(qb.all()[last][0].pk) #pk of Cp2kGeoOptDdecWorkChain
    qb.append(WorkCalculation, filters={'label':{'==':'Cp2kCellOptWorkChain'}}, output_of='workflow', tag='cellopt')
    qb.append(WorkCalculation, filters={'label':{'==':'Cp2kDftBaseWorkChain'}}, output_of='cellopt', tag='base')
    qb.order_by({WorkCalculation:{'ctime':'desc'}})
    try:
        localpath =  qb.all()[last][0].out.retrieved.get_abs_path()+'/path/aiida.out'
        print('%s\tpk: %s\tloacalpath: %s' %(structure,pk_work,localpath))
        file_out = dir_out + structure + "_" + workflow_label + "_steps.out"
        with open(file_out, 'w+') as fout:
            with redirect_stdout(fout):
                print_header()
                print_steps(localpath)
    except:
        print('%s\tpk: %s\tSTILL RUNNING'%(structure, pk_work))
