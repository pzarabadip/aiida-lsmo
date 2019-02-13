import pandas as pd
import matplotlib.pyplot as plt
import os.path

with open('list-321.list') as f:
    ids=f.read().splitlines()

#Look what I'm checking from the final mex
for id in ids:
    steps_file='parse_Cp2kGeoOpt/test1-0/{}_test1-0_steps.out'.format(id)
    if not os.path.exists(steps_file):
        mex = 'not_started'
    else:
        data = pd.read_csv(steps_file,sep=' ')
        if data['#step'].value_counts()[0]<5: #count the number of "0" in #steps
            mex = 'not_finished'
        else:
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
    print('%s %s'%(id,mex))
