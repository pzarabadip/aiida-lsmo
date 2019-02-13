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

def plot_steps(file_out, structure):
    """ Function to plot the graph of the energy """
    steps = np.genfromtxt(file_out, delimiter="", comments="#",usecols = (0))
    energy = np.genfromtxt(file_out, delimiter="", comments="#",usecols = (1))

    # Take min and max, neglecting spikes
    min_energy = +9999.
    max_energy = -9999.
    spike_thr = 2. #Ha
    for i in range(len(energy)):
        if i==0 or i==len(energy)-1 or \
           abs(energy[i]-energy[i-1])+abs(energy[i]-energy[i+1]) < spike_thr:
           if energy[i]<min_energy:
               min_energy = energy[i]
           if energy[i]>max_energy:
               max_energy = energy[i]

    energy_shifted= [ (x-min_energy) for x in energy ]
    max_energy_shifted = max_energy-min_energy

    fig, ax = plt.subplots(figsize=[8, 4.5])
    ax.set(title='Robust cell optimization of: '+ structure,
           xlabel='Steps',
           ylabel='Energy (Hartree)',
           ylim=[-0.01*max_energy_shifted, +1.01*max_energy_shifted],
           )
    ax.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:,.3f}'))
    #make coloured background
    startindex=[i for i, x in enumerate(steps) if x == 0]
    startindex.append(len(steps))
    if len(startindex) > 2:
        ax.axvspan(startindex[1], startindex[2]-1, ymin=0, ymax=1, color='red', alpha=0.2)
    if len(startindex) > 3:
        ax.axvspan(startindex[2], startindex[3]-1, ymin=0, ymax=1, color='orange', alpha=0.2)
    if len(startindex) > 4:
        ax.axvspan(startindex[3], startindex[4]-1, ymin=0, ymax=1, color='yellow', alpha=0.2)
    if len(startindex) == 6:
        ax.axvspan(startindex[4], startindex[5]-1,   ymin=0, ymax=1, color='green', alpha=0.2)
    #print energy profile
    ax.plot(energy_shifted,color='blue',marker='o',markersize=3,linewidth=1)
    ax.grid()

    fig.savefig(file_out[:-4]+".png",dpi=300)
    plt.close(fig)
    return

# User settings
last = 0 #select the Nth last result in time
workflow_label = 'test1-0'
with open('list-321.list') as f:
    structure_labels=f.read().splitlines()
# structure_labels=['13182N3','14000N2'] #test
# General settings
stage_name = ['Stage0_Energy','Stage1_CellOpt','Stage2_MD','Stage3_GeoOpt','Stage4_CellOpt']
dir_out="./parse_Cp2kGeoOpt/"
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
dir_out="./parse_Cp2kGeoOpt/"+workflow_label+"/"
if not os.path.exists(dir_out):
    os.makedirs(dir_out)

for structure in structure_labels:
    stage_localpath = ["INCOMPLETE"] * 5
    stage_pk = ["INCOMPLETE"] * 5

    # get first the pk of the main workflow Cp2kGeoOptDdecWorkChain
    qb = QueryBuilder()
    qb.append(StructureData, filters={'label': {'==':structure}}, tag='structure')
    qb.append(WorkCalculation, filters={'label':{'==':workflow_label}}, output_of='structure', tag='workflow')
    if len(qb.all())==0:
        pk_work=0
        mex = 'NEVER_STARTED'
    else:
        pk_work=str(qb.all()[last][0].pk)
        pkfile = open(dir_out + structure + "_" + workflow_label + "_pk.out","w+")
        # Search for the path for the 5 steps: energy, cell_opt, md, geo_opt, cell_opt. Store them.
        for istage in range(5):
            if istage == 0 or stage_localpath[istage-1] != "INCOMPLETE":
                qb = QueryBuilder()
                qb.append(StructureData, filters={'label': {'==':structure}}, tag='structure')
                qb.append(WorkCalculation, filters={'label':{'==':workflow_label}}, output_of='structure', tag='workflow')
                qb.append(WorkCalculation, filters={'label':{'==':'Cp2kRobustGeoOptWorkChain'}}, output_of='workflow', tag='robustgeoopt')
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
            file_out = dir_out + structure + "_" + workflow_label + "_steps.out"
            with open(file_out, 'w+') as fout:
                with redirect_stdout(fout):
                    print_header()
                    for i in range(completed_stage+1):
                        try:
                            print_steps(stage_localpath[i])
                        except:
                            mex = "WARNING: Something weird happened reading the steps in cp2k.out"
            if completed_stage>0: #plot_steps crashing with only energy value
                plot_steps(file_out, structure)

    print('%s\tpk: %s\tCompleted: %s' %(structure,pk_work,mex))
