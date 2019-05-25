from __future__ import print_function
from aiida.orm.calculation.work import WorkCalculation
import sys
import os
import re
import glob
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
            min0=data['energy(eV/atom)'].sort_values().iloc[i]
            min1=data['energy(eV/atom)'].sort_values().iloc[i+1]
            if min1-min0<0.1: #Tollerance (Ha) to decide if a min is weird
                break
        ediff=data['energy(eV/atom)'].iloc[-1]-min0
        ethr = 0.005 #Tollerance (Ha) to decide if the difference is relevant
        if ediff>ethr:
            mex = 'energy[-1]_>_min(energy)+{0}eV/atom_ediff={1:.3f}eV/atom'.format(ethr,ediff)
        else:
            diff=data['energy(eV/atom)'].iloc[-1]-data['energy(eV/atom)'].iloc[0]
            steps=data['#step'].iloc[-1]
            mex = 'Ediff={:.3f}eV/atom_Steps={:d}'.format(diff,steps)
    return mex

def plot_steps(stepsfile,structure_label,workflow_label):
    """ Function to plot the graph of the energy """
    df = pd.read_csv(stepsfile,sep=' ')
    steps = df['#step'].tolist()
    energy = df['energy(eV/atom)'].tolist()

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
    ax.set(title='{} {}'.format(structure_label,workflow_label),
           xlabel='Steps',
           ylabel='Relative Energy (eV/atom)',
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
last=0
structure_labels = [s.split('/')[-1].split('.')[0] for s in  glob.glob("./cifs/*.cif")]
workflow_label = 'mof_directcellopt_ddec_batch1'
dir_out="./out_{}/".format(workflow_label)
if not os.path.exists(dir_out):
    os.makedirs(dir_out)

# Print header on screen
print('{:<15s} {:<12s} {:<10s} {:<8s} {:<7s} {:<30s} {:<10s} {:s}'.format('Structure','pk_WorkChain','pk_Cellopt','BG-alpha','BG-beta','CellOpt_check','pk_CifDDEC','dir_cellopt'))
for structure_label in structure_labels:
    # Clear variables
    localpath_CellOpt = 'none'
    pk_DirectCellOptDdec='none'
    pk_CellOpt ='none'
    bandgap = [np.nan,np.nan]
    cellopt_check='none'
    pk_cifdata = 'none'
    # Start QB: first if WC started, then if ended
    qb = QueryBuilder() # WARNING: This QB appears to be VERY slow, taking many seconds every time qb.all() is called!
    qb.append(StructureData, filters={'label': structure_label}, tag='inp_struct')
    qb.append(WorkCalculation, filters={'label':workflow_label},output_of='inp_struct', tag='wf')
    qb.order_by({WorkCalculation:{'ctime':'desc'}})
    if len(qb.all())>0:
        pk_DirectCellOptDdec=str(qb.all()[last][0].pk)
        qb.append(WorkCalculation, filters={'label':{'==':'Cp2kCellOptWorkChain'}}, output_of='wf', tag='cellopt')
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
            stepsfile = "{}/{}_steps.out".format(dir_out,structure_label)
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
                plot_steps(stepsfile, structure_label, workflow_label)
                # Check the energy convergence of the whole wf
                cellopt_check=directcellopt_check(stepsfile)
            except:
                cellopt_check='WARNING'
            # See if the final DDEC cif exist and print it
            qb.append(CifData, edge_filters={'label': 'output_structure'},output_of='wf',tag='cifdata')
            try:
                cifdata = qb.all()[last][0]
                pk_cifdata = str(cifdata.pk)
                cifile = open('./{}/{}_cellopt.cif'.format(dir_out,structure_label),'w+')
                print(cifdata.values,file=cifile)
                cifile.close()
            except:
                pass

    #Print info on screen (for localpath print the directory, not the file aiida.out)
    if localpath_CellOpt[-9:]=='aiida.out':
        localpath_CellOpt=localpath_CellOpt[:-9]
    print('{:<15s} {:<12s} {:<10s} {:<8.3f} {:<7.3f} {:<30s} {:<10s} {:s}'.format(structure_label,pk_DirectCellOptDdec,pk_CellOpt,bandgap[0],bandgap[1],cellopt_check,pk_cifdata,localpath_CellOpt))
