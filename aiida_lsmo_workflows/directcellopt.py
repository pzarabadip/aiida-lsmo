from aiida.orm import CalculationFactory, DataFactory
from aiida.orm.code import Code
from aiida.orm.data.base import Float
from aiida.work import workfunction as wf
from aiida.work.run import submit
from aiida.work.workchain import WorkChain, ToContext, Outputs
from copy import deepcopy

# subworkflows
from aiida_cp2k.workflows import Cp2kCellOptWorkChain
from aiida_ddec.workflows import DdecCp2kChargesWorkChain

# data objects
ArrayData = DataFactory('array')
CifData = DataFactory('cif')
ParameterData = DataFactory('parameter')
RemoteData = DataFactory('remote')
StructureData = DataFactory('structure')

def dict_merge(dct, merge_dct):
    """ Taken from https://gist.github.com/angstwad/bf22d1822c38a92ec0a9
    Recursive dict merge. Inspired by :meth:``dict.update()``, instead of
    updating only top-level keys, dict_merge recurses down into dicts nested
    to an arbitrary depth, updating keys. The ``merge_dct`` is merged into
    ``dct``.
    :param dct: dict onto which the merge is executed
    :param merge_dct: dct merged into dct
    :return: None
    """
    import collections
    for k, v in merge_dct.iteritems():
        if (k in dct and isinstance(dct[k], dict)
                and isinstance(merge_dct[k], collections.Mapping)):
            dict_merge(dct[k], merge_dct[k])
        else:
            dct[k] = merge_dct[k]

spin = {
        "H"  : 0.0 / 2.0,
        "Li" : 0.0 / 2.0,  # oxidation state +1
        "Be" : 0.0 / 2.0,  # oxidation state +2
        "B"  : 0.0 / 2.0,
        "C"  : 0.0 / 2.0,
        "N"  : 0.0 / 2.0,
        "O"  : 0.0 / 2.0,
        "F"  : 0.0 / 2.0,
        "Na" : 0.0 / 2.0,  # oxidation state +1
        "Mg" : 0.0 / 2.0,  # oxidation state +2
        "Al" : 0.0 / 2.0,  # oxidation state +3
        "Si" : 0.0 / 2.0,
        "P"  : 0.0 / 2.0,
        "S"  : 0.0 / 2.0,
        "Cl" : 0.0 / 2.0,
        "K"  : 0.0 / 2.0,  # oxidation state +1
        "Ca" : 0.0 / 2.0,  # oxidation state +2
        "Sc" : 0.0 / 2.0,  # oxidation state +1
        "Ti" : 0.0 / 2.0,  # oxidation state +4
        "V"  : 3.0 / 2.0,  # oxidation state +2, high spin
        "Cr" : 4.0 / 2.0,  # oxidation state +2, high spin
        "Mn" : 5.0 / 2.0,  # oxidation state +2, high spin
        "Fe" : 4.0 / 2.0,  # oxidation state +2, high spin
        "Co" : 3.0 / 2.0,  # oxidation state +2, high spin
        "Ni" : 2.0 / 2.0,  # oxidation state +2, high spin
        "Cu" : 1.0 / 2.0,  # oxidation state +2, high spin
        "Zn" : 0.0 / 2.0,  # oxidation state +2
        "Zr" : 0.0 / 2.0,  # oxidation state +4
        }

@wf
def guess_multiplicity(structure):
    multiplicity = 1
    all_atoms = structure.get_ase().get_chemical_symbols()
    for key, value in spin.iteritems():
        multiplicity += all_atoms.count(key) * value * 2.0
    multiplicity = int(round(multiplicity))
    multiplicity_dict = {'FORCE_EVAL': {'DFT': {'MULTIPLICITY' :multiplicity}}}
    if multiplicity != 1:
        multiplicity_dict['FORCE_EVAL']['DFT']['UKS'] = True
    return ParameterData(dict=multiplicity_dict)

@wf
def multiply_unit_cell (struct, threshold):
    """Resurns the multiplication factors (tuple of 3 int) for the cell vectors
    that are needed to respect: min(perpendicular_width) > threshold."""
    from math import cos, sin, sqrt, pi
    import numpy as np
    # angle between vectors
    def angle(v1,v2):
        return np.arccos(np.dot(v1,v2) / (np.linalg.norm(v1)*np.linalg.norm(v2)))

    threshold = threshold.value / 2.0

    a = np.linalg.norm(struct.cell[0])
    b = np.linalg.norm(struct.cell[1])
    c = np.linalg.norm(struct.cell[2])

    alpha = angle(struct.cell[1], struct.cell[2])
    beta = angle(struct.cell[0], struct.cell[2])
    gamma = angle(struct.cell[0], struct.cell[1])

    # first step is computing cell parameters according to  https://en.wikipedia.org/wiki/Fractional_coordinates
    # Note: this is the algoriGthm implemented in Raspa (framework.c/UnitCellBox). There also is a simpler one but it is less robust.
    v = sqrt(1-cos(alpha)**2-cos(beta)**2-cos(gamma)**2+2*cos(alpha)*cos(beta)*cos(gamma))
    cell=np.zeros((3,3))
    cell[0,:] = [a, 0, 0]
    cell[1,:] = [b*cos(gamma), b*sin(gamma),0]
    cell[2,:] = [c*cos(beta), c*(cos(alpha)-cos(beta)*cos(gamma))/(sin(gamma)),c*v/sin(gamma)]
    cell=np.array(cell)

    # diagonalizing the cell matrix: note that the diagonal elements are the perpendicolar widths because ay=az=bz=0
    diag = np.diag(cell)
    repeat = tuple(int(i) for i in np.ceil(threshold/diag*2.))
    return StructureData(ase=struct.get_ase().repeat(repeat)).store()

class Cp2kDirectCellOptWorkChain(WorkChain):
    """Workchain that for a given matherial will optimize the structure
    and compute ddec charges."""
    @classmethod
    def define(cls, spec):
        super(Cp2kDirectCellOptWorkChain, cls).define(spec)

        # structure, adsorbant, pressures
        spec.input('structure', valid_type=StructureData)
        spec.input("min_cell_size", valid_type=Float, default=Float(10))

        # cp2k
        spec.input('cp2k_code', valid_type=Code)
        spec.input('cp2k_parameters', valid_type=ParameterData, default=ParameterData(dict={}))
        spec.input("_cp2k_options", valid_type=dict, default=None, required=False)
        spec.input('cp2k_parent_folder', valid_type=RemoteData, default=None, required=False)

        # settings
        spec.input('_guess_multiplicity', valid_type=bool, default=False)

        # workflow
        spec.outline(
            cls.run_cellopt,
            cls.parse_cellopt,
            cls.return_results,
        )

        # specify the outputs of the workchain
        spec.output('output_structure', valid_type=StructureData, required=False)

    def run_cellopt(self):
        """Run CELL_OPT calculation."""
        self.ctx.structure = multiply_unit_cell(self.inputs.structure, self.inputs.min_cell_size)
        cp2k_parameters = deepcopy(self.inputs.cp2k_parameters.get_dict())

       # Expand the unit cell so that: min(perpendicular_width) > threshold
        inputs = {
            'code'          : self.inputs.cp2k_code,
            'structure'     : self.ctx.structure,
            #'min_cell_size' : self.inputs.min_cell_size, not necessary!
            '_options'      : self.inputs._cp2k_options,
            '_label'        : "Cp2kCellOptWorkChain",
        }

        # Trying to guess the multiplicity of the system
        if self.inputs._guess_multiplicity:
            self.report("Guessing multiplicity")
            dict_merge(cp2k_parameters, guess_multiplicity(self.ctx.structure).get_dict())

        # Take parameters
        self.ctx.cp2k_parameters = ParameterData(dict=cp2k_parameters)
        inputs['parameters'] = self.ctx.cp2k_parameters

        # Create the calculation process and launch it
        running = submit(Cp2kCellOptWorkChain, **inputs)
        self.report("pk: {} | Running Cp2kCellOptWorkChainn to optimize cell and geometry".format(running.pid))
        return ToContext(cellopt_calc=Outputs(running))

    def parse_cellopt(self):
        """Extract optimized structure and put it into self.ctx.structure"""
        self.ctx.structure = self.ctx.cellopt_calc['output_structure']

    def return_results(self):
        self.out('output_structure', self.ctx.structure)
        self.report("Cp2kGeoOptDdecWorkChain is completed.")
