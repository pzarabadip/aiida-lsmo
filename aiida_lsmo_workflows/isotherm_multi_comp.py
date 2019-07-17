"""
Multi Component Isotherm WorkChain
"""
from __future__ import absolute_import
import os
import six
from copy import deepcopy

from aiida.common import AttributeDict
from aiida.plugins import CalculationFactory, DataFactory
from aiida.orm import Code, Dict, Float, Int, List, Str
from aiida.engine import submit
from aiida.engine import ToContext, WorkChain, workfunction, if_, while_, append_


from aiida_lsmo_workflows.utils.multiply_unitcell import multiply_unit_cell
from aiida_lsmo_workflows.utils.pressure_points import choose_pressure_points

CifData = DataFactory('cif')
ParameterData = DataFactory('dict')
SinglefileData = DataFactory('singlefile')
RaspaCalculation = CalculationFactory('raspa')
ZeoppCalculation = CalculationFactory('zeopp.network')
NetworkParameters = DataFactory('zeopp.parameters')


class MultiCompIsothermWorkChain(WorkChain):
    """
    The MultiCompIsothermWorkChain is designed to perform zeo++ and
    RASPA calculations for multi-component mixtures.
    """

    @classmethod
    def define(cls, spec):
        super(MultiCompIsothermWorkChain, cls).define(spec)
        """
        Define workflow specification.
        This is the most important method of a Workchain, which defines the
        inputs it takes, the logic of the execution and the outputs
        that are generated in the process.
        """
        # General inputs
        spec.input('structure', valid_type=CifData, required=True, help='Input structure in cif format')
        # Zeopp inputs
        spec.input('zeopp_code', valid_type=Code, required=False)
        # It is being provided for each component so we can remove it from here.
        # spec.input("zeopp_probe_radius", valid_type=Float, required=False)
        spec.input("zeopp_atomic_radii", valid_type=SinglefileData, required=False)

        # Raspa inputs
        spec.input("raspa_code", valid_type=Code, required=False)
        spec.input("raspa_parameters", valid_type=ParameterData, required=False)
        spec.input("raspa_isotherm_dynamic", valid_type=bool, default=False, required=False, non_db=True)
        spec.input("raspa_isotherm_full", valid_type=bool, default=False, required=False, non_db=True)
        spec.input("raspa_usecharges", valid_type=bool, default=False, required=False, non_db=True)
        spec.input("raspa_charge_cif", valid_type=bool, default=False, required=False, non_db=True)

        spec.input("raspa_cutoff", valid_type=Float, default=Float(12.0), required=False)
        # spec.input("raspa_minKh", valid_type=Float, default=Float(1e-10), required=False)
        # spec.input("raspa_minKh_sel", valid_type=Float, default=Float(5.0), required=False)
        spec.input("raspa_verbosity", valid_type=Int, default=Int(10), required=False)
        spec.input("raspa_widom_cycle_mult", valid_type=Int, default=Int(10), required=False)
        spec.input("raspa_num_of_components", valid_type=Int, default=Int(2), required=False)
        spec.input_namespace("raspa_comp", valid_type=dict,required=False, dynamic=True)

        spec.input("raspa_pressure_min", valid_type=Float, default=Float(0.1), required=False)
        spec.input("raspa_pressure_max", valid_type=Float, default=Float(1.0), required=False)
        # TODO: Here we need to decide for different compoentns
        spec.input("raspa_molsatdens", valid_type=Float, required=False)
        spec.input("raspa_gcmc_press_precision", valid_type=Float, required=False)
        spec.input("raspa_gcmc_press_maxstep", valid_type=Float, required=False)
        spec.input("selected_pressures", valid_type=list, required=False, non_db=True)

        # Zeopp Extra Inputs
        spec.input("zeopp_accuracy", valid_type=Str, default=Str('DEF'), required=False)
        spec.input("zeopp_block_samples_A3", valid_type=Int, default=Int(100), required=False) #100 samples / Ang^3: accurate for all the structures
        spec.input("zeopp_volpo_samples_UC", valid_type=Int, default=Int(100), required=False) #100k samples, may need more for structures bigger
        spec.input("zeopp_lcd_max", valid_type=Float, default=Float(15.0), required=False)
        spec.input("zeopp_pld_min", valid_type=Float, default=Float(3.9), required=False)

        # Workflow
        spec.outline(
            cls.setup,
            cls.run_pore_dia_zeopp,
            if_(cls.should_run_zeopp_full)(
                cls.run_zeopp_full,
                cls.inspect_zeopp_calc,
                cls.init_raspa_widom,
                cls.run_raspa_widom,
                if_(cls.should_run_gcmc)(
                    cls.init_raspa_gcmc,
                    while_(cls.should_run_another_gcmc)(
                        cls.run_raspa_gcmc,
                        cls.parse_raspa_gcmc,
                    ),
                ),
            ),
            cls.return_results,
            )


        # to be returned
        spec.outputs.dynamic = True

    def setup(self):
        """
        Initialize variables and setup screening protocol!
        """

        self.ctx.zeopp_options = {
            "resources": {
                "num_machines": 1,
                "tot_num_mpiprocs": 1,
            },
            "max_memory_kb" : 2000000,
            "max_wallclock_seconds": 1 * 30 * 60,
            "withmpi": False,
        }

        self.ctx.raspa_options = {
            "resources": {
                "num_machines": 1,
                "tot_num_mpiprocs": 1,
            },
            "max_memory_kb": 200000,
            "max_wallclock_seconds": 2 * 60 * 60,
            "withmpi": False,
        }


        self.ctx.raspa_comp = AttributeDict(self.inputs.raspa_comp)

        #TODO: Adding zeopp parameters too.
        self.ctx.raspa_pressure_min = self.inputs.raspa_pressure_min
        self.ctx.raspa_pressure_max = self.inputs.raspa_pressure_max


    def run_pore_dia_zeopp(self):
        """
        It performs the zeopp pore diameter calculation at high accuracy
        """
        # Required inputs
        params = {
            'res' : [self.inputs.structure.label + '.res'],
            'ha'  : self.inputs.zeopp_accuracy.value
        }

        parameters = NetworkParameters(dict=params)

        inputs = {
            'code'      : self.inputs.zeopp_code,
            'structure' : self.inputs.structure,
            'parameters': parameters,
            'metadata'  :{
                'options':  self.ctx.zeopp_options,
                'label'  : 'Pore Diameter Calculation',
            }
        }

        # Use default zeopp atomic radii only if a .rad file is not specified
        try:
            inputs['atomic_radii'] = self.inputs.zeopp_atomic_radii
            self.report("Zeopp will use atomic radii from the .rad file")
        except:
            self.report("Zeopp will use default atomic radii")

        # Creating the calculation process and submit it
        res = self.submit(ZeoppCalculation, **inputs)
        self.report("pk: {} | Running zeo++ pore diameter calculation".format(res.pk))
        return ToContext(zeopp_res=res)

    def should_run_zeopp_full(self):
        """
        It uses largest included sphere (Di or LCD) and largest free sphere
        (Df or PLD) as pre-screenig descriptors to pass or reject the
        structure.
        Note: The default setting is LCD < 7.0 A and PLD > 3.9 A.
        """
        lcd_lim = self.inputs.zeopp_lcd_max.value
        pld_lim = self.inputs.zeopp_pld_min.value
        lcd_current = self.ctx.zeopp_res.outputs.output_parameters.get_dict()['Largest_included_sphere']
        pld_current = self.ctx.zeopp_res.outputs.output_parameters.get_dict()['Largest_free_sphere']

        if (lcd_current < lcd_lim) and (pld_current > pld_lim):
            self.report("{} is a suitable structure for further investigation".format(self.inputs.structure.label))
            return True
        else:
            self.report("{} does not look like promising: stop".format(self.inputs.structure.label))
            return False

    def run_zeopp_full(self):
        """
        It calculated the surface area, pore volume, and
        possible block pockets for pre-selected materials.
        """

        for key, value in self.ctx.raspa_comp.items():
            if key in list(self.inputs.raspa_comp):
                comp_name = value.name
                probe_radius = value.radius

                params = {
                        'sa'   : [probe_radius,
                                  probe_radius,
                                  self.inputs.zeopp_block_samples_A3.value,
                                  self.inputs.structure.label + '_' + comp_name + '.sa'],
                        'volpo': [probe_radius,
                                  probe_radius,
                                  self.inputs.zeopp_volpo_samples_UC.value,
                                  self.inputs.structure.label + '_' + comp_name + '.volpo'],
                        'block': [probe_radius,
                                  self.inputs.zeopp_block_samples_A3.value,
                                  self.inputs.structure.label + '_' + comp_name + '.block'],
                        'ha'   :  self.inputs.zeopp_accuracy.value
                }

                # Required inputs
                inputs = {
                    'code'      : self.inputs.zeopp_code,
                    'structure' : self.inputs.structure,
                    'parameters': NetworkParameters(dict=params).store(),
                    'metadata'  :{
                        'options':  self.ctx.zeopp_options,
                        'label'  : 'Zeo++ Calculation',
                    }
                }

                try:
                    inputs['atomic_radii'] = self.inputs.zeopp_atomic_radii
                    self.report("Zeopp will use atomic radii from the .rad file")
                except:
                    self.report("Zeopp will use default atomic radii")

                # Creating the calculation process and submit it
                zeopp_full = self.submit(ZeoppCalculation, **inputs)
                zeopp_label = 'zeopp_{}'.format(comp_name)
                self.report("pk: {} | Running Zeo++ calculation using {} probe".format(zeopp_full.pk,comp_name))
                self.to_context(**{zeopp_label: zeopp_full})


    def inspect_zeopp_calc(self):
        for key, value in self.ctx.raspa_comp.items():
            if key in list(self.inputs.raspa_comp):
                comp_name = value.name
                zeopp_label = 'zeopp_{}'.format(comp_name)
                self.report("Checking if {} job finished OK".format(zeopp_label))
                assert self.ctx[zeopp_label].is_finished_ok

    # def should_run_zeopp_block(self):
    #     """
    #     Base on the outcomes of SaVolpo calculation here we decide if we need to
    #     calcualte block pockets or no.
    #     """
    #     for key, value in self.ctx.raspa_comp.items():
    #         if key in list(self.inputs.raspa_comp):
    #             comp_name = value.name
    #
    #             sv_label = 'sv_{}'.format(comp_name)
    #
    #             output_sa = self.ctx[sv_label].outputs.output_parameters.get_dict()
    #             ASA = output_sa['ASA_A^2']
    #             NASA = output_sa['NASA_A^2']
    #
    #             if (NASA > 0.0) and (ASA > 0.0):
    #                 self.report("We need to calcualte the block pocket file for {}!".format(comp_name))
    #                 self.to_context(**{sv_label: True})
    #                 # return True
    #             elif (NASA == 0.0) and (ASA > 0.0):
    #                 self.report("All pores are accessible to {}!".format(comp_name))
    #                 self.to_context(**{sv_label: False})
    #                 # return False
    #             else:
    #                 self.report("This structure does not present any accessibility to! : stop")
    #                 pass

    # def run_block_zeopp(self):
    #     """
    #     It calculates the block pockets only on structures with non-accessible void fraction.
    #     """
    #     for key, value in self.ctx.raspa_comp.items():
    #         if key in list(self.inputs.raspa_comp):
    #             comp_name = value.name
    #             probe_radius = value.radius
    #
    #
    #
    #             params = {
    #                     'block': [probe_radius,
    #                               self.inputs.zeopp_block_samples_A3.value,
    #                               self.inputs.structure.label + '_' + comp_name + '.block'],
    #                     'ha'   : self.inputs.zeopp_accuracy.value
    #             }
    #
    #             # Required inputs
    #             inputs = {
    #                 'code'      : self.inputs.zeopp_code,
    #                 'structure' : self.inputs.structure,
    #                 'parameters': NetworkParameters(dict=params).store(),
    #                 'metadata'  :{
    #                     'options':  self.ctx.zeopp_options,
    #                     'label'  : 'Pore Block Calculation',
    #                 }
    #             }
    #
    #             try:
    #                 inputs['atomic_radii'] = self.inputs.zeopp_atomic_radii
    #                 self.report("Zeopp will use atomic radii from the .rad file")
    #             except:
    #                 self.report("Zeopp will use default atomic radii")
    #
    #
    #             # Creating the calculation process and submit it
    #             block = self.submit(ZeoppCalculation, **inputs)
    #             block_label = 'block_{}'.format(comp_name)
    #             self.report("pk: {} | Running Pore Block calculation using {} probe".format(block.pk,comp_name))
    #             #return ToContext(zeopp_sv=append_(sv))
    #             self.to_context(**{block_label: block})

    # def inspect_block_calc(self):
    #     for key, value in self.ctx.raspa_comp.items():
    #         if key in list(self.inputs.raspa_comp):
    #             comp_name = value.name
    #             block_label = 'block_{}'.format(comp_name)
    #             # self.report("Checking if {} job finished OK".format(sv_label))
    #             assert self.ctx[block_label].is_finished_ok

    def init_raspa_widom(self):
        """
        It generates the ParameterData for a Raspa and use the blocking spheres if there is any.
        With current plugin we can have more than one component to be used for Widom insertion.
        """

        # Create a deepcopy of the user parameters, to modify before submission
        self.ctx.raspa_parameters = deepcopy(self.inputs.raspa_parameters.get_dict())

        for key, value in self.ctx.raspa_comp.items():
            if key in list(self.inputs.raspa_comp):
                comp_name = value.name
                mol_def = value.mol_def
                bp_label = '_'.join((self.inputs.structure.label,comp_name))
                self.ctx.raspa_parameters['Component'][comp_name] = self.ctx.raspa_parameters['Component'].pop(key)
                self.ctx.raspa_parameters["Component"][comp_name]["MoleculeDefinition"] = mol_def
                self.ctx.raspa_parameters["Component"][comp_name]["WidomProbability"] = 1.0
                self.ctx.raspa_parameters["Component"][comp_name]["BlockPocketsFileName"] = bp_label

        ucs = multiply_unit_cell(self.inputs.structure, self.inputs.raspa_cutoff.value)
        self.ctx.raspa_parameters['GeneralSettings']['UnitCells'] = "{} {} {}".format(ucs[0], ucs[1], ucs[2])

        # Turn on charges if requested
        if self.inputs.raspa_usecharges:
            self.ctx.raspa_parameters['GeneralSettings']['ChargeMethod'] = "Ewald"
            self.ctx.raspa_parameters['GeneralSettings']['EwaldPrecision'] = 1e-6
            if self.inputs.raspa_charge_cif:
                self.ctx.raspa_parameters['GeneralSettings']['UseChargesFromCIFFile'] = "yes"
            else:
                self.ctx.raspa_parameters['GeneralSettings']['ChargeFromChargeEquilibration'] = "yes"
        else:
            self.ctx.raspa_parameters['GeneralSettings']['ChargeMethod'] = "None"

        # CORRECT the settings to have only Widom insertion
        self.ctx.raspa_parameters["GeneralSettings"]["SimulationType"] = "MonteCarlo"
        self.ctx.raspa_parameters["GeneralSettings"]["NumberOfInitializationCycles"] = 0
        self.ctx.raspa_parameters["GeneralSettings"]["NumberOfCycles"] = self.inputs.raspa_widom_cycle_mult.value * self.inputs.raspa_parameters.get_dict()["GeneralSettings"]["NumberOfCycles"]
        self.ctx.raspa_parameters["GeneralSettings"]["PrintPropertiesEvery"] = int(self.ctx.raspa_parameters["GeneralSettings"]["NumberOfCycles"] / self.inputs.raspa_verbosity.value)
        self.ctx.raspa_parameters["GeneralSettings"]["PrintEvery"] = int(1e6) #never
        return

    def run_raspa_widom(self):
        """
        It runs a Widom calculation in Raspa.
        """
        # Create the inputs dictionary
        inputs = {
            'framework'  : {self.inputs.structure.label : self.inputs.structure},
            'code'       : self.inputs.raspa_code,
            'parameters' : ParameterData(dict=self.ctx.raspa_parameters).store(),
            'metadata'  :{
                'options':  self.ctx.raspa_options,
                'label'  : 'RASPA Widom Calculation',
            }
        }

        inputs['block_pocket'] = {}
        self.ctx.number_blocking_spheres = {}

        for key, value in self.ctx.raspa_comp.items():
            if key in list(self.inputs.raspa_comp):
                comp_name = value.name
                zeopp_label = 'zeopp_{}'.format(comp_name)
                bp_label = '_'.join((self.inputs.structure.label,comp_name))
                bp_path = os.path.join(self.ctx[zeopp_label].outputs.retrieved._repository._get_base_folder().abspath, bp_label + '.block')

                with open(bp_path, 'r') as block_file:
                    self.ctx.number_blocking_spheres[comp_name] = int(block_file.readline().strip())
                    if self.ctx.number_blocking_spheres[comp_name] > 0:
                        # bp_comp = SinglefileData(file=bp_path)
                        inputs['block_pocket'][bp_label] = self.ctx[zeopp_label].outputs.block

                        self.report("({}) Blocking spheres are present for ({}) and used for Raspa".format(self.ctx.number_blocking_spheres[comp_name],comp_name))
                    else:
                        self.report("No blocking spheres found for ({})".format(comp_name))

        # Create the calculation process and launch it
        widom = self.submit(RaspaCalculation, **inputs)
        self.report("pk: {} | Running Raspa Widom for the Henry coefficient".format(widom.pk))

        return ToContext(raspa_widom=widom)

    def should_run_gcmc(self):
        """
        It decides if we should calculate loadings. The decision should be taken
        on two criteria:
        1 - Selectivity inversion: The selectivity toward desired component should
        be above unity. However, we can have a tolerance (user-defined) for not missing
        structures around unity.

        2 - Selectivity above a user-defined threshold.
        """
        output = self.ctx.raspa_widom.outputs.output_parameters.get_dict()
        #TODO: make it more comprehensive to flag possible selectivity inversion.
        # Desired compound should be set to comp1 in the input.
        Kh_desired_comp = output[self.inputs.structure.label]['components'][self.ctx.raspa_comp.comp1.name]['henry_coefficient_average']
        for key, value in self.ctx.raspa_comp.items():
            if key in list(self.inputs.raspa_comp):
                comp_name = value.name
                if comp_name != self.ctx.raspa_comp.comp1.name:
                    Kh_comp = output[self.inputs.structure.label]['components'][comp_name]['henry_coefficient_average']
                    if Kh_comp < Kh_desired_comp:
                        self.report("Ideal selectivity larger than the threshold: compute binary loading")
                        return True
                    else:
                        self.report("Ideal selectivty either below threshold or inverted!: don't compute isotherm")
                        return False

        # TODO: using count and len() for more than two components.


    def init_raspa_gcmc(self):
        """
        Initialization of input for raspa gcmc.
        It generates a pressure list based on our decision.
        1 - Full isotherm simulation based on Daniele work.
        2 - Partial isotherm simulation.
        """
        # Initializate counter and set restart to None
        self.ctx.current_p_index = 0
        # self.ctx.restart_raspa_calc = None

        if (self.inputs.raspa_isotherm_dynamic) and (self.inputs.raspa_isotherm_full):
            # Estimate the total loading qsat and choose the pressure points
            satDens = self.inputs.raspa_molsatdens.value #(mol/l)
            poreVol = self.ctx.zeopp_sv.outputs.output_parameters.get_dict()['POAV_cm^3/g'] #(cm3/g = l/kg)
            self.ctx.estimated_qsat = satDens * poreVol
            self.ctx.pressures = choose_pressure_points(
                # TODO: It currently takes the greatest henry coefficient.
                Kh = self.ctx.kh_for_iso, #(mol/kg/Pa)
                qst = self.ctx.estimated_qsat, #(mol/kg_frame)
                dpa = self.inputs.raspa_gcmc_press_precision.value, #(kg*Pa/mol)
                dpmax = self.inputs.raspa_gcmc_press_maxstep.value, #(Pa)
                pmax = self.inputs.raspa_pressure_max.value, #(Pa)
                pmin = self.inputs.raspa_pressure_min.value, #(Pa)
                dynamic = self.inputs.raspa_isotherm_dynamic,
                full = self.inputs.raspa_isotherm_full,
                )
        elif (not self.inputs.raspa_isotherm_dynamic) and (self.inputs.raspa_isotherm_full):
            self.ctx.pressures = choose_pressure_points(
                pmin = self.self.inputs.raspa_pressure_min.value,
                dpa = self.inputs.raspa_gcmc_press_precision.value,
                pmax = self.inputs.raspa_pressure_max.value,
                dynamic = self.inputs.raspa_isotherm_dynamic,
                full = self.inputs.raspa_isotherm_full,
            )
        else:
            self.ctx.pressures = choose_pressure_points(
                selected_pressures = self.inputs.selected_pressures,
                dynamic = self.inputs.raspa_isotherm_dynamic,
                full = self.inputs.raspa_isotherm_full,
            )

        self.report("{} number of points are chosen for isotherm construction".format(len(self.ctx.pressures)))

        # CORRECT the parameters to perform GCMC
        self.ctx.raspa_parameters["GeneralSettings"]["NumberOfInitializationCycles"] = self.inputs.raspa_parameters.get_dict()["GeneralSettings"]["NumberOfInitializationCycles"]
        self.ctx.raspa_parameters["GeneralSettings"]["NumberOfCycles"] = self.inputs.raspa_parameters.get_dict()["GeneralSettings"]["NumberOfCycles"]
        self.ctx.raspa_parameters["GeneralSettings"]["PrintPropertiesEvery"] = int(1e6) #never
        self.ctx.raspa_parameters["GeneralSettings"]["PrintEvery"] = int(self.ctx.raspa_parameters["GeneralSettings"]["NumberOfCycles"]/self.inputs.raspa_verbosity.value)
        #TODO: Using for loop here to prevent repeating everthing.

        for key, value in self.ctx.raspa_comp.items():
            if key in list(self.inputs.raspa_comp):
                comp_name = value.name
                mol_frac = value.mol_fraction
                self.ctx.raspa_parameters["Component"][comp_name]["MolFraction"] = float(mol_frac)
                del self.ctx.raspa_parameters["Component"][comp_name]["WidomProbability"]
                self.ctx.raspa_parameters["Component"][comp_name]["TranslationProbability"] = 0.5
                self.ctx.raspa_parameters["Component"][comp_name]["RotationProbability"] = 0.5
                self.ctx.raspa_parameters["Component"][comp_name]["ReinsertionProbability"] = 0.5
                self.ctx.raspa_parameters["Component"][comp_name]["SwapProbability"] = 1.0
        return

    def should_run_another_gcmc(self):
        """
        We run another raspa calculation only if the current iteration is
        smaller than the total number of pressures we want to compute.
        """
        return self.ctx.current_p_index < len(self.ctx.pressures)


    def run_raspa_gcmc(self):
        """
        It runs a binary GCMC calculation at T and lower pressure.
        """

        pressure = self.ctx.pressures[self.ctx.current_p_index]
        self.ctx.raspa_parameters['System']['hkust1']['ExternalPressure'] = pressure
        # Create the input dictionary
        inputs = {
            'framework'  : {self.inputs.structure.label:self.inputs.structure},
            'code'       : self.inputs.raspa_code,
            'parameters' : ParameterData(dict=self.ctx.raspa_parameters).store(),
            'metadata'  :{
                'options':  self.ctx.raspa_options,
                'label'  : 'RASPA GCMC Calculation',
            }
        }

        inputs['block_pocket'] ={}

        for key, value in self.ctx.raspa_comp.items():
            if key in list(self.inputs.raspa_comp):
                comp_name = value.name
                zeopp_label = 'zeopp_{}'.format(comp_name)
                bp_label = '_'.join((self.inputs.structure.label,comp_name))
                if self.ctx.number_blocking_spheres[comp_name] > 0:
                    # bp_comp = self.ctx[zeopp_label].outputs.block
                    inputs['block_pocket'][bp_label] = self.ctx[zeopp_label].outputs.block

        # Create the calculation process and launch it
        gcmc = self.submit(RaspaCalculation, **inputs)
        self.report("pk: {} | Running Raspa GCMC at p(bar)={:.3f} ({} of {})".format(gcmc.pk, pressure/1e5, self.ctx.current_p_index+1, len(self.ctx.pressures)))
        return ToContext(raspa_gcmc=gcmc)

    def parse_raspa_gcmc(self):
        """
        Extract the pressure and loading average of the last completed raspa
        calculation.
        """
        pressure = self.ctx.raspa_parameters['System'][self.inputs.structure.label]['ExternalPressure']/1e5
        output_gcmc = self.ctx.raspa_gcmc.outputs.output_parameters.get_dict()

        if self.ctx.current_p_index == 0:
            self.ctx.loading = {}
            for key, value in self.ctx.raspa_comp.items():
                if key in list(self.inputs.raspa_comp):
                    comp_name = value.name
                    self.ctx.loading[comp_name] = []


        for key, value in self.ctx.raspa_comp.items():
            if key in list(self.inputs.raspa_comp):
                comp_name = value.name
                conv1 = output_gcmc[self.inputs.structure.label]['components'][comp_name]['conversion_factor_molec_uc_to_mol_kg']
                loading_average_comp = conv1 * output_gcmc[self.inputs.structure.label]['components'][comp_name]['loading_absolute_average']
                loading_dev_comp = conv1 * output_gcmc[self.inputs.structure.label]['components'][comp_name]['loading_absolute_dev']
                self.ctx.loading[comp_name].append([pressure, loading_average_comp, loading_dev_comp])

        # Update counter and parent folder for restart
        self.ctx.current_p_index += 1
        # self.ctx.restart_raspa_calc = self.ctx.raspa_gcmc['retrieved_parent_folder']
        return

    # This part also needs a revision! I should double check the parsers.
    def return_results(self):
        """
        Attach the results to the output.
        """
        result_dict = {}
        result_dict['loading'] = {}
        result_dict['isotherm_loading_header'] = ['Pressure(bar)', 'Loading_average(mol/kg)', 'Loading_deviation(mol/kg)']
        result_dict['henry_coefficient_average'] = {}
        result_dict['henry_coefficient_dev'] = {}
        result_dict['adsorption_energy_average'] = {}
        result_dict['adsorption_energy_dev'] = {}
        result_dict['isotherm_loading'] = {}
        result_dict['isotherm_enthalpy'] = {}
        result_dict['number_blocking_spheres'] = {}
        result_dict['POAV_Volume_fraction'] = {}
        result_dict['PONAV_Volume_fraction'] = {}
        result_dict['POAV'] = {}
        result_dict['GASA'] = {}
        result_dict['VASA'] = {}
        result_dict['NGASA'] = {}
        result_dict['NVASA'] = {}
        result_dict['Channel_surface_area'] = {}
        result_dict['Pocket_surface_area'] = {}
        result_dict['number_blocking_spheres'] = {}
        result_dict['Density_unit'] = "g/cm^3"
        result_dict['POAV_unit'] = "cm^3/g"
        result_dict['GASA_unit'] = "m^2/g"
        result_dict['VASA_unit'] = "m^2/cm^3"
        result_dict['NGASA_unit'] = "m^2/g"
        result_dict['NVASA_unit'] = "ASA_m^2/cm^3"
        result_dict['Channel_surface_area_unit'] = "A^2"
        result_dict['Pocket_surface_area_unit'] = "A^2"


        # Zeopp section
        # TODO: Making it to return results for different components.

        for key, value in self.ctx.raspa_comp.items():
            if key in list(self.inputs.raspa_comp):
                comp_name = value.name
                zeopp_label = 'zeopp_{}'.format(comp_name)
                output_zeo = self.ctx[zeopp_label].outputs.output_parameters.get_dict()
                result_dict['POAV_Volume_fraction'][comp_name] = output_zeo['POAV_Volume_fraction']
                result_dict['PONAV_Volume_fraction'][comp_name] = output_zeo['PONAV_Volume_fraction']
                result_dict['POAV'][comp_name] = output_zeo['POAV_cm^3/g']
                result_dict['GASA'][comp_name] = output_zeo['ASA_m^2/g']
                result_dict['VASA'][comp_name] = output_zeo['ASA_m^2/cm^3']
                result_dict['NGASA'][comp_name] = output_zeo['NASA_m^2/g']
                result_dict['NVASA'][comp_name] = output_zeo['NASA_m^2/cm^3']
                result_dict['Channel_surface_area'][comp_name] = output_zeo['Channel_surface_area_A^2']
                result_dict['Pocket_surface_area'][comp_name] = output_zeo['Pocket_surface_area_A^2']
                result_dict['number_blocking_spheres'][comp_name] = self.ctx.number_blocking_spheres[comp_name]
                # TODO: Take it out of the loop!
                result_dict['Density'] = output_zeo['Density']

        # RASPA Section
        output_widom = self.ctx.raspa_widom.outputs.output_parameters.get_dict()
        result_dict['temperature'] = self.ctx.raspa_parameters["System"][self.inputs.structure.label]["ExternalTemperature"]
        result_dict['temperature_unit'] = "K"
        result_dict['henry_coefficient_units'] = output_widom[self.inputs.structure.label]['components'][self.ctx.raspa_comp.comp1.name]['henry_coefficient_units']
        result_dict['adsorption_energy_units'] = output_widom[self.inputs.structure.label]['components'][self.ctx.raspa_comp.comp1.name]['adsorption_energy_widom_units']

        # result_dict['conversion_factor_molec_uc_to_cm3stp_cm3'] = output_gcmc['hkust1']['components'][self.inputs.raspa_comp1.value]['conversion_factor_molec_uc_to_cm3stp_cm3']
        # result_dict['conversion_factor_molec_uc_to_gr_gr'] = output_gcmc['hkust1']['components'][self.inputs.raspa_comp1.value]['conversion_factor_molec_uc_to_gr_gr']
        # result_dict['conversion_factor_molec_uc_to_mol_kg'] = output_gcmc['hkust1']['components'][self.inputs.raspa_comp1.value]['conversion_factor_molec_uc_to_mol_kg']


        for key, value in self.ctx.raspa_comp.items():
            if key in list(self.inputs.raspa_comp):
                comp_name = value.name
                comp = output_widom[self.inputs.structure.label]['components'][comp_name]
                result_dict['henry_coefficient_average'][comp_name] = comp['henry_coefficient_average']
                result_dict['henry_coefficient_dev'][comp_name] = comp['henry_coefficient_dev']
                result_dict['adsorption_energy_average'][comp_name] = comp['adsorption_energy_widom_average']
                result_dict['adsorption_energy_dev'][comp_name] = comp['adsorption_energy_widom_dev']
                result_dict['isotherm_loading'][comp_name] = self.ctx.loading[comp_name]
                result_dict['number_blocking_spheres'][comp_name] = self.ctx.number_blocking_spheres[comp_name]

        self.out("results", ParameterData(dict=result_dict).store())
        # self.out('blocking_spheres', self.ctx.zeopp_block['block'])
        self.report("Workchain completed successfully")
        return
# EOF
