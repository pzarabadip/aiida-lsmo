from __future__ import print_function
from aiida.orm.calculation.work import WorkCalculation
import sys
import os
import re
import matplotlib
import pandas as pd
matplotlib.use('Agg') #avoids ssh problems, but can not use .show()
import matplotlib.pyplot as plt
import numpy as np
sys.path.append('/home/daniele/Programs/cp2k_utils/parse_steps')
from parser_utils import print_header, print_steps
ParameterData = DataFactory('parameter')
StructureData = DataFactory('structure')
CifData = DataFactory('cif')

from contextlib import contextmanager
@contextmanager
def redirect_stdout(target):
    original = sys.stdout
    sys.stdout = target
    yield
    sys.stdout = original

def get_bandgap(stepsfile):
    '''Given a cp2k.out file where the band gap is compute, returns alpha and beta bandgap.'''
    bg = [np.nan,np.nan]
    file = open(stepsfile, "r")
    first=True
    for line in file:
        if re.search('HOMO', line):
            if first:
                bg[0]=float(line.split()[6])
                first=False
            else:
                bg[1]=float(line.split()[6])
    return bg

def get_startindex(l):
    '''Take a list (of steps), and decide starting indexes and final number of steps.'''
    startindex=[]
    for i in range(len(l)):
        if i==0:
            startindex.append(0)
        else:
            if l[i]<=l[i-1]:
                startindex.append(i)
    startindex.append(len(l))
    return startindex

def cellopt_check(stepsfile):
    '''Given the _steps.out file, evaluates the geometry convergence'''
    data = pd.read_csv(stepsfile,sep=' ')
    if data['#step'].iloc[-1]==1000:
        mex = 'cellopt_reached_1000'
    else:
        #Look for the min energy value, getting rid of weird low energy data
        for i in range(10):
            min0=data['energy(Ha)'].sort_values().iloc[i]
            min1=data['energy(Ha)'].sort_values().iloc[i+1]
            if min1-min0<0.1: #Tollerance (Ha) to decide if a min is weird
                break
        ediff=data['energy(Ha)'].iloc[-1]-min0
        ethr = 0.005 #Tollerance (Ha) to decide if the difference is relevant
        if ediff>ethr:
            mex = 'energy[-1]_>_min(energy)+{0}Ha_ediff={1:.3f}Ha'.format(ethr,ediff)
        else:
            mex = 'EVERYTHING_FINE'
    return mex

def plot_steps(stepsfile, structure):
    """ Function to plot the graph of the energy """

    steps = np.genfromtxt(stepsfile, delimiter="", comments="#",usecols = (0))
    energy = np.genfromtxt(stepsfile, delimiter="", comments="#",usecols = (1))

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
    startindex=get_startindex(steps)
    if len(startindex) > 2: #Print only Stage1_CellOpt
        ax.axvspan(startindex[1], startindex[2]-1, ymin=0, ymax=1, color='red', alpha=0.2)
    if len(startindex) > 3: #Print also Stage2_MD
        ax.axvspan(startindex[2], startindex[3]-1, ymin=0, ymax=1, color='orange', alpha=0.2)
    if len(startindex) > 4: #Print also Stage3_CellOpt
        ax.axvspan(startindex[3], startindex[4]-1, ymin=0, ymax=1, color='green', alpha=0.2)
    #print energy profile
    ax.plot(energy_shifted,color='blue',marker='o',markersize=3,linewidth=1)
    ax.grid()

    fig.savefig(stepsfile[:-4]+".png",dpi=300)
    plt.close(fig)
    return

# User settings ************************************************************************************
last = 1 #select the Nth last result in time
workflow_label = 'test2-smearing3'
with open('list-smearing.list') as f:
    structure_labels=f.read().splitlines()
structure_labels=['18081N2'] #

# General settings
stage_name = ['Stage0_Energy','Stage1_CellOpt','Stage2_MD','Stage3_CellOpt']
dir_out="./parse_Cp2kCellOpt/"
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
dir_out="./parse_Cp2kCellOpt/"+workflow_label+"/"
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
# Print header on screen
print('Structure  pk_Cp2kCellOptDdecWorkChain  Completed       min(BandGap)  pk_CifDDEC  CellOpt_check ')
for structure in structure_labels:
    stage_localpath = ["INCOMPLETE"] * 4
    stage_pk = ["INCOMPLETE"] * 4
    pk_CifDDEC = 'none'
    # get first the pk of the main workflow Cp2kCellOptDdecWorkChain and of the final cif
    qb = QueryBuilder()
    qb.append(StructureData, filters={'label': {'==':structure}}, tag='structure')
    qb.append(WorkCalculation, filters={'label':{'==':workflow_label}}, output_of='structure', tag='workflow')
    if len(qb.all())==0:
        pk_work=0
        mex = 'NEVER_STARTED'
        mexcheck='never_started'
        minbg = np.inf
    else:
        pk_work=str(qb.all()[last][0].pk)
        pkfile = open(dir_out + structure + "_" + workflow_label + "_pk.out","w+")
        try: #get the final cif
            qb.append(CifData, edge_filters={'label': 'output_structure'},output_of='workflow')
            qb.order_by({WorkCalculation:{'ctime':'desc'}})
            pk_CifDDEC = str(qb.all()[last][0].pk)
        except:
            pass

        # Search for the path for the 5 steps: energy, cell_opt, md, geo_opt, cell_opt. Store them.
        for istage in range(4):
            if istage == 0 or stage_localpath[istage-1] != "INCOMPLETE": #stop if the previous stage was incomplete
                qb = QueryBuilder()
                qb.append(StructureData, filters={'label': {'==':structure}}, tag='structure')
                qb.append(WorkCalculation, filters={'label':{'==':workflow_label}}, output_of='structure', tag='workflow')
                qb.append(WorkCalculation, filters={'label':{'==':'Cp2kRobustCellOptWorkChain'}}, output_of='workflow', tag='robustcellopt')
                if istage == 0:
                    qb.append(WorkCalculation, filters={'label':{'==':stage_name[0]}}, output_of='robustcellopt', tag='dftbase')
                else:
                    qb.append(WorkCalculation, filters={'label':{'==':stage_name[istage]}}, output_of='robustcellopt', tag='stage')
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
        completed_stage = 3 - stage_localpath.count("INCOMPLETE")
        if completed_stage == -1:
            mex = "None"
        else:
            mex = stage_name[completed_stage]
        bgfile = open(dir_out + structure + "_" + workflow_label + "_bg.out","w+")
        minbg = np.inf
        stepsfile = dir_out + structure + "_" + workflow_label + "_steps.out"
        with open(stepsfile, 'w+') as fout:
            with redirect_stdout(fout):
                print_header()
                for i in range(completed_stage+1):
                    if os.path.exists(stage_localpath[i]):
                        # Print steps for cp2k.out to _steps.out w/redirect
                        print_steps(stage_localpath[i])
                        # Print bandgap and store the min bg excluding the Stage0 where the geometry can be weird
                        bandgap = get_bandgap(stage_localpath[i])
                        if i>0 and min(bandgap)< minbg:
                            minbg = min(bandgap)
                        print('%s\t%f\t%f'%(stage_name[i],bandgap[0],bandgap[1]), file=bgfile)
                    else:
                        mex = "WARNING_cp2k.out_missing"
                bgfile.close()
        # Plot graph only if at least Stage1 is finished
        if completed_stage>=1:
            plot_steps(stepsfile, structure)
        # Check the energy convergence of the whole wf
        if completed_stage==3:
           mexcheck=cellopt_check(stepsfile)
        else:
           mexcheck='still-running_or_crashed'
    #Print info on screem
    print('%-10s %-28s %-15s %-13.3f %-11s %-s' %(structure,pk_work,mex,minbg,pk_CifDDEC,mexcheck))
