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

def directcellopt_check(stepsfile):
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
            diff=data['energy(Ha)'].iloc[-1]-data['energy(Ha)'].iloc[0]
            steps=data['#step'].iloc[-1]
            mex = 'Ediff={:.3f}Ha_Steps={:d}'.format(diff,steps)
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
    ax.set(title='Final cell optimization of: '+ structure,
           xlabel='Steps',
           ylabel='Energy (Hartree)',
           ylim=[-0.01*max_energy_shifted, +1.01*max_energy_shifted],
           )
    ax.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:,.3f}'))
    #print energy profile
    ax.plot(energy_shifted,color='blue',marker='o',markersize=3,linewidth=1)
    ax.grid()

    fig.savefig(stepsfile[:-4]+".png",dpi=300)
    plt.close(fig)
    return

# User settings
last = 0 #select the Nth last result in time
for calc in ['ot','smearing']:
    if calc=='ot':
        workflow_label = 'test3-OT'
        prevWorkflow = 'test2-0'
        with open('../cof_test2/list-OT.list') as f:
            structure_labels=f.read().splitlines()
        structure_labels=structure_labels[140:]
    elif calc=='smearing':
        workflow_label = 'test3-smearing'
        prevWorkflow = 'test2-smearing'
        with open('../cof_test2/list-smearing.list') as f:
            structure_labels=f.read().splitlines()

    # General settings
    dir_out="./parse_DirectCp2kCellOpt/"
    if not os.path.exists(dir_out):
        os.makedirs(dir_out)
    dir_out="./parse_DirectCp2kCellOpt/"+workflow_label+"/"
    if not os.path.exists(dir_out):
        os.makedirs(dir_out)
    # Print header on screen
    print('{:<9s} {:<12s} {:<10s} {:<8s} {:<7s} {:<30s} {:<10s} {:s}'.format('Structure','pk_WorkChain','pk_Cellopt','BG-alpha','BG-beta','CellOpt_check','pk_CifDDEC','dir_cellopt'))
    for structure in structure_labels:
        # Clear variables
        localpath_CellOpt = 'none'
        pk_DirectCellOptDdec='none'
        pk_CellOpt ='none'
        bandgap = [np.nan,np.nan]
        cellopt_check='none'
        pk_CifData = 'none'
        # Start QB: first if WC started, then if ended
        qb = QueryBuilder()
        qb.append(StructureData, filters={'label': structure}, tag='inp_struct')
        qb.append(WorkCalculation, filters={'label':prevWorkflow},output_of='inp_struct', tag='wf')
        qb.append(CifData, edge_filters={'label': 'output_structure'},output_of='wf',tag='cifdata')
        qb.append(StructureData, descendant_of='cifdata', tag='opt_struct')
        qb.append(WorkCalculation, filters={'label':{'==':workflow_label}}, output_of='opt_struct', tag='wf2')
        qb.order_by({WorkCalculation:{'ctime':'desc'}})
        if len(qb.all())>0:
            pk_DirectCellOptDdec=str(qb.all()[last][0].pk)
            qb.append(WorkCalculation, filters={'label':{'==':'Cp2kCellOptWorkChain'}}, output_of='wf2', tag='cellopt')
            qb.append(WorkCalculation, filters={'label':{'==':'Cp2kDftBaseWorkChain'}}, output_of='cellopt', tag='dftbase')
            qb.append(JobCalculation,output_of='dftbase',tag='calc')
            qb.order_by({WorkCalculation:{'ctime':'desc'}})
            try:
                localpath_CellOpt =  qb.all()[last][0].out.retrieved.get_abs_path()+'/path/aiida.out'
                pk_CellOpt = str(qb.all()[last][0].pk)
            except:
                localpath_CellOpt=None
                pass
            # If WC ended, parse stuff
            if localpath_CellOpt!='none':
                stepsfile = dir_out + structure + "_" + workflow_label + "_steps.out"
                with open(stepsfile, 'w+') as fout:
                    with redirect_stdout(fout):
                        print_header()
                        if os.path.exists(localpath_CellOpt):
                            # Print steps for cp2k.out to _steps.out w/redirect
                            print_steps(localpath_CellOpt)
                            # Print bandgap and store the min bg excluding the Stage0 where the geometry can be weird
                            bandgap = get_bandgap(localpath_CellOpt)
                        else:
                            bandgap = [404,404] #Missing file
                try:
                    # Plot graph only if at least Stage1 is finished
                    plot_steps(stepsfile, structure)
                    # Check the energy convergence of the whole wf
                    cellopt_check=directcellopt_check(stepsfile)
                except:
                    cellopt_check='WARNING'
                # See if the final DDEC cif exist
                qb.append(CifData, edge_filters={'label': 'output_structure'},output_of='wf2',tag='cifdata2')
                try:
                    pk_CifData = str(qb.all()[last][0].pk)
                except:
                    pass

        #Print info on screen (for localpath print the directory, not the file aiida.out)
        if localpath_CellOpt[-9:]=='aiida.out':
            localpath_CellOpt=localpath_CellOpt[:-9]
        print('{:<9s} {:<12s} {:<10s} {:<8.3f} {:<7.3f} {:<30s} {:<10s} {:s}'.format(structure,pk_DirectCellOptDdec,pk_CellOpt,bandgap[0],bandgap[1],cellopt_check,pk_CifData,localpath_CellOpt))
