from __future__ import print_function
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

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

def plot_steps(stepsfile1,stepsfile2, structure, pngdir):
    """ Function to plot the graph of the energy """

    steps = np.genfromtxt(stepsfile1, delimiter="", comments="#",usecols = (0))
    energy = np.genfromtxt(stepsfile1, delimiter="", comments="#",usecols = (1))
    energy2 = np.genfromtxt(stepsfile2, delimiter="", comments="#",usecols = (1))

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
    energy2_shifted= [ (x-min_energy) for x in energy2 ]
    max_energy_shifted = max_energy-min_energy
    # ylim[ymin,ymax] are a little complicate but effective to leave some margin!
    ymax=+1.01*max_energy_shifted
    ymin=energy2_shifted[-1]-0.01*(max_energy_shifted-energy2_shifted[-1])

    fig, ax = plt.subplots(figsize=[8, 4.5])
    ax.set(title='Robust cell optimization of: '+ structure,
           xlabel='Steps',
           ylabel='Energy (Hartree)',
           ylim=[ymin,ymax],
           )
    ax.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:,.3f}'))
    #make coloured background
    startindex=get_startindex(steps)
    ax.axvspan(startindex[1], startindex[2]-1, ymin=0, ymax=1, color='red', alpha=0.2)
    ax.axvspan(startindex[2], startindex[3]-1, ymin=0, ymax=1, color='orange', alpha=0.2)
    ax.axvspan(startindex[3], startindex[4]-1, ymin=0, ymax=1, color='green', alpha=0.2)
    ax.axvspan(startindex[4], startindex[4]+len(energy2), ymin=0, ymax=1, color='cyan', alpha=0.2)
    #print energy profile
    ax.plot(energy_shifted+energy2_shifted,
            color='blue',marker='o',markersize=3,linewidth=1)
    ax.grid()

    fig.savefig("{}/{}.png".format(pngdir,structure),dpi=300)
    plt.close(fig)
    return

# Main
if True: #choose OT or choose smearing***********************************************************************************
    workflow_label = 'test3-OT'
    prevWorkflow = 'test2-0'
    with open('../cof_test2/list-OT.list') as f:
        structure_labels=f.read().splitlines()
else:
    workflow_label = 'test3-smearing'
    prevWorkflow = 'test2-smearing'
    with open('../cof_test2/list-smearing.list') as f:
        structure_labels=f.read().splitlines()

pngdir="./parse_plot_Robust+Direct/"
if not os.path.exists(pngdir):
    os.makedirs(pngdir)

for cof in structure_labels:
    stepsfile1="../cof_test2/parse_Cp2kCellOpt/{}/{}_{}_steps.out".format(prevWorkflow,cof,prevWorkflow)
    stepsfile2="../cof_test3/parse_DirectCp2kCellOpt/{}/{}_{}_steps.out".format(workflow_label,cof,workflow_label)
    steps1_exists=os.path.exists(stepsfile1)
    steps2_exists=os.path.exists(stepsfile2)
    print('Plotting {}: file1 {}, file2 {}'.format(cof,steps1_exists,steps2_exists),end='')
    if steps1_exists and steps2_exists:
        try:
            plot_steps(stepsfile1, stepsfile2, cof, pngdir)
        except:
            print(' *WARNING*',end='')
    print('')
