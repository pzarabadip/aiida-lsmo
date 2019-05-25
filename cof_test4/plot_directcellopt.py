import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

ids=[
    ['05001N2','OT'],
    ['07000N2','OT'],
    ['07001N2','smearing'],
    ['07002N2','OT'],
    ['07010N3','OT'],
    ['08000N3','OT'],
    ['08011N2','OT'],
    ['08012N2','OT'],
    ['08013N2','OT'],
    ['08020N2','OT'],
    ['08030N2','OT'],
    ['09010N2','OT'],
]


labels=['CellOpt+AIMD+CellOpt','Direct CellOpt']
for id in ids:
    file0 = '../cof_test2/parse_Cp2kCellOpt/test2-0_list-{}.list/{}_test2-0_steps.out'.format(id[1],id[0])
    data0 = pd.read_csv(file0,sep=' ')
    file1 = './parse_steps/test4fidis-OT/{}_test4fidis-OT_steps.out'.format(id[0])
    data1 = pd.read_csv(file1,sep=' ')
    #file2 = './parse_steps/directcellopt_randomized/{}_directcellopt_randomized_steps.out'.format(id)
    #data2 = pd.read_csv(file2,sep=' ')
    fig, ax = plt.subplots()
    ax.plot(data0['energy(eV/atom)'], label=labels[0], linewidth = 4)
    ax.plot(data1['energy(eV/atom)'], label=labels[1], linewidth = 3)
    #plt.plot(data2['energy(Ha)'], label=labels[2], linewidth = 2)
    ax.set_title(id[0])
    ax.grid()
    ax.set_xlabel("Optimization step")
    ax.set_ylabel("Absolute energy (eV/atom)")
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    #plt.legend()
    #plt.show()
    plt.savefig("./parse_steps/{}_compare.png".format(id[0]),bbox_inches='tight',dpi=500)
    plt.clf()
