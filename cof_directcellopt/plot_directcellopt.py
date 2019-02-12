import pandas as pd
import matplotlib.pyplot as plt

with open('list-directcellopt.list') as f:
    ids=f.read().splitlines()

labels=['RobustCellOpt','DirectCellOpt','Rand+DirectCellOpt']
for id in ids:
    file0 = '../cof_test1/parse_Cp2kGeoOpt/test1-0/{}_test1-0_steps.out'.format(id)
    data0 = pd.read_csv(file0,sep=' ')
    file1 = './parse_steps/directcellopt/{}_directcellopt_steps.out'.format(id)
    data1 = pd.read_csv(file1,sep=' ')
    file2 = './parse_steps/directcellopt_randomized/{}_directcellopt_randomized_steps.out'.format(id)
    data2 = pd.read_csv(file2,sep=' ')
    plt.plot(data0['energy'], label=labels[0], linewidth = 4)
    plt.plot(data1['energy(Ha)'], label=labels[1], linewidth = 3)
    plt.plot(data2['energy(Ha)'], label=labels[2], linewidth = 2)
    plt.title(id)
    plt.grid()
    plt.legend()
    #plt.show()
    plt.savefig("./parse_steps/{}compare.png".format(id))
    plt.clf()
