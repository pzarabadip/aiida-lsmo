import pandas as pd
from aiida_qeq.calculations.qeq import QeqCalculation

CifData = DataFactory('cif')

df = pd.read_csv("../cof_test2/pk_final.csv")
wc_label = "aiida_qeq Qeq with GMP parameters"

for i in df.index:
    cif_label = df.at[i,'structure']
    q = QueryBuilder()
    q.append(CifData, filters={'label': cif_label}, tag='cif-inp')
    q.append(QeqCalculation, filters={'label': wc_label}, output_of='cif-inp', tag='qeq')
    if len(q.all())>0:
        pk_qeqcalc=q.all()[0][0].pk
    else:
        pk_cifqeq='missing'
    q.append(CifData, output_of='qeq', tag='cif-qeq')
    if len(q.all())>0:
        pk_cifqeq=q.all()[0][0].pk
        print("CifData.label= {} with Qeq >>> CifData.pk= {} ".format(
              cif_label,pk_cifqeq))
    else:
        print("CifData.label= {} failed QeqCalculation.pk= {}".format(
              cif_label,pk_qeqcalc))
