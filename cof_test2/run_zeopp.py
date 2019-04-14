from __future__ import absolute_import
from __future__ import print_function
import os
import pandas as pd
from aiida.common.example_helpers import test_and_get_code
from aiida.orm import DataFactory, CalculationFactory
from aiida.work import workfunction as wf
ParameterData = DataFactory('parameter')
StructureData = DataFactory('structure')
CifData = DataFactory('cif')
ZeoppCalculation = CalculationFactory('zeopp.network')
NetworkParameters = DataFactory('zeopp.parameters')

@wf
def StructureData2CifData(structuredata):
    """Converts StructureData to ASE, enforcing a filename with .cif extension."""
    from aiida.orm.data.cif import cif_from_ase, pycifrw_from_cif, ase_loops
    from aiida.common.utils import HiddenPrints
    import tempfile
    cif = cif_from_ase(structuredata.get_ase())
    cifdata = CifData()
    with tempfile.NamedTemporaryFile(suffix='.cif') as f:
        with HiddenPrints():
            f.write(pycifrw_from_cif(cif, loops=ase_loops).WriteOut())
        f.flush()
        cifdata.set_file(f.name)

    cifdata.description = "Converted from StructureData"
    return cifdata

# Get and Test the code
zeopp_code = test_and_get_code('zeopp@deneb', expected_code_type='zeopp.network')

# Prepare input parameters
d = {
    'ha': True,                    # Using high accuracy (mandatory!)
    'res': True,                   # Max included, free and incl in free sphere
    'chan': 1.2,                   # Small probe to identify channels
    'sa': [1.86, 1.86, 100000],    # Nitrogen probe to compute surface
    'vol': [0.0, 0.0, 1000000],    # Geometric pore volume
    'volpo': [1.86, 1.86, 100000], # Nitrogen probe to compute PO pore volume
    'psd': [1.2, 1.2, 10000]       # Small probe to compute the pore size distr
}
parameters = NetworkParameters(dict=d)

# Read the list of PKs
df=pd.read_csv("./pk_final.csv")

# Import the structures, set up the calculation and run it
for pk in df.StructureData: # Starting geometries
    cifdata = StructureData2CifData(load_node(int(pk)))
#for pk in df.CifData: # Optimized geometries
#  if pk != 'none':
#    cifdata = load_node(int(pk))
    calc = zeopp_code.new_calc()
    calc.label = "aiida_zeopp std calculation"
    calc.description = "Standard geometry characterization for porous materials"
    calc.set_max_wallclock_seconds(24 * 60 * 60)
    calc.set_withmpi(False)
    calc.set_resources({"num_machines": 1})
    calc.use_parameters(parameters)
    calc.use_structure(cifdata)

    # Optional: use radii file
    use_radii=False
    if use_radii:
        SinglefileData = DataFactory('singlefile')
        atomic_radii = SinglefileData(file=os.path.join(tests.TEST_DIR, 'UFFrad'))
        calc.use_atomic_radii(atomic_radii)

    calc.store_all()
    calc.submit()
    print("submitted calculation: calc_pk= {} cif_pk= {}".format(calc.dbnode.pk,cifdata.pk))
