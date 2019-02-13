from __future__ import print_function
from aiida.orm.data.structure import StructureData
from aiida.orm.calculation.work import WorkCalculation
import sys
import os
import re
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
workflow_label = 'test1-0'
with open('list-321.list') as f:
    structure_labels=f.read().splitlines()
# Search for the path for the 4 steps: energy, md, geo_opt, cell_opt. Store them.
for structure in structure_labels:
        qb = QueryBuilder()
        qb.append(StructureData, filters={'label': {'==':structure}}, tag='structure')
        qb.append(WorkCalculation, filters={'label':{'==':workflow_label}}, output_of='structure', tag='workflow')
        qb.append(WorkCalculation, filters={'label':{'==':'DdecCp2kChargesWorkChain'}}, output_of='workflow', tag='ddecWork')
        qb.append(WorkCalculation, filters={'label':{'==':'Cp2kDftBaseWorkChain'}}, output_of='ddecWork', tag='dftbase')
        qb.append(JobCalculation,output_of='dftbase',tag='calc')
        qb.order_by({WorkCalculation:{'ctime':'desc'}})
        try:
            path =  qb.all()[last][0].out.retrieved.get_abs_path()+'/path/aiida.out'
        except:
            path='none'

        homo_alpha='x.xxxxxx'
        homo_beta='x.xxxxxx'
        # Print file with pk and local directory
        if not path=='none':
            file = open(path, "r")
            first=True
            for line in file:
                if re.search('HOMO', line):
                    if first:
                        homo_alpha=line.split()[6]
                        first=False
                    else:
                        homo_beta=line.split()[6]
            file.close()
        print('%s %9s %9s %s'%(structure,homo_alpha,homo_beta,path))
