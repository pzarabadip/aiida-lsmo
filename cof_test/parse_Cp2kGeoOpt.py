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

    min_energy = min(energy)
    energy_shifted= [ (x-min_energy) for x in energy ]

    fig, ax = plt.subplots(figsize=[8, 4.5])
    ax.set(ylabel='Energy (Hartree)',
           xlabel='Steps',
           title='Robust cell optimization of: '+ structure)
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

#seairch for the path for the 4 steps: energy, md, geo_opt, cell_opt. Store them.
#workflow_name = 'Test18CofsNoFixAngles'
workflow_name = 'Test18CofsFixAngles50MD'
stage_name = ['Stage0_Energy','Stage1_CellOpt','Stage2_MD','Stage3_GeoOpt','Stage4_CellOpt']
last = 0 #select the Nth last result in time
#with open('cif_list.txt') as finp:
#    labels=finp.read().splitlines()
structure_labels = [
'05000N',
'05001N',
'07000N',
'07001N',
'07002N',
'07010N',
'07011N',
'07012N',
'07013N',
'08000N',
'08010N',
'08011N',
'08012N',
'08013N',
'08020N',
'08030N',
'09000N',
'09010N',
'13000N',
'13010N',
'13011N',
'13012N',
'13020N',
'13030N',
'13040N',
'13050N',
'13051N',
'13060N',
'13070N',
'13071N',
'13072N',
'13073N',
'13074N',
'13075N',
'13076N',
'13077N',
'13078N',
'13080N',
'13090N',
'13100N',
'13101N',
'13110N',
'13120N',
'13121N',
'13122N',
'13123N',
'13130N',
'13140N',
'13141N',
'13142N',
'13150N',
'13160N',
'13161N',
'13170N',
'13180N',
'13181N',
'13182N',
]

dir_out="./parse_Cp2kGeoOpt/"
if not os.path.exists(dir_out):
    os.makedirs(dir_out)

for structure in structure_labels:
    stage_localpath = ["INCOMPLETE"] * 5
    stage_pk = ["INCOMPLETE"] * 5
    pkfile = open(dir_out + structure + "_" + workflow_name + "_pk.out","w+")
    for istage in range(5):
        if istage == 0 or stage_localpath[istage-1] != "INCOMPLETE":
            qb = QueryBuilder()
            qb.append(StructureData, filters={'label': {'==':structure}}, tag='structure')
            qb.append(WorkCalculation, filters={'label':{'==':workflow_name}}, output_of='structure', tag='workflow')
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
    print('%s\tCompleted: %s' %(structure,mex))

    file_out = dir_out + structure + "_" + workflow_name + "_steps.out"
    with open(file_out, 'w+') as fout:
        with redirect_stdout(fout):
            print_header()
            for i in range(completed_stage+1):
                print_steps(stage_localpath[i])
    if completed_stage>0: # at least Stage1 completed
        plot_steps(file_out, structure)
