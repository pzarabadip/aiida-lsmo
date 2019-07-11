"""
Selectivity and Working Capacity WorkChain
"""
from __future__ import absolute_import
import os
import six
from copy import deepcopy

from aiida.common import AttributeDict
from aiida.plugins import CalculationFactory, DataFactory
from aiida.orm import Code, Dict, Float, Int, List, Str
from aiida.engine import submit
from aiida.engine import ToContext, WorkChain, workfunction, if_, while_

from aiida_lsmo_workflows.utils.multiply_unitcell import multiply_unit_cell
from aiida_lsmo_workflows.utils.pressure_points import choose_pressure_points

CifData = DataFactory('cif')
ParameterData = DataFactory('dict')
SinglefileData = DataFactory('singlefile')
RaspaCalculation = CalculationFactory('raspa')
ZeoppCalculation = CalculationFactory('zeopp.network')
NetworkParameters = DataFactory('zeopp.parameters')

class SeparationWorkChain(WorkChain):
    """
    The SeparationWorkChain is designed to calculate loadings
    of a multi-component gas mixture at two pressure points
    and return the selectivity and working capacity/deliverable_capacity.
    """

    @classmethod
    def define(cls, spec):
        super(SeparationWorkChain, cls).define(spec)
        """Define workflow specification.

        This is the most important method of a Workchain, which defines the
        inputs it takes, the logic of the execution and the outputs
        that are generated in the process.
        """
        # General inputs
        spec.input('structure', valid_type=CifData, required=True, help='Input structure in cif format')
        # Zeopp inputs
        spec.input('zeopp_code', valid_type=Code, required=False)
        spec.input("zeopp_probe_radius", valid_type=Float, required=False)
        spec.input("zeopp_atomic_radii", valid_type=SinglefileData, required=False)


        # Raspa inputs
        spec.input("raspa_code", valid_type=Code, required=False)
        spec.input("raspa_parameters", valid_type=ParameterData, required=False)
        spec.input("raspa_isotherm_dynamic", valid_type=bool, default=False, required=False, non_db=True)
        spec.input("raspa_isotherm_full", valid_type=bool, default=False, required=False, non_db=True)
        #spec.input("raspa_usecharges", valid_type=bool, default=False, required=False)
        spec.input("raspa_cutoff", valid_type=Float, default=Float(12.0), required=False)
        spec.input("raspa_minKh", valid_type=Float, default=Float(1e-10), required=False)
        spec.input("raspa_minKh_sel", valid_type=Float, default=Float(5.0), required=False)
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

        #Workflow
        # spec.outline(
        #     cls.setup,
        # #     cls.run_pore_dia_zeopp,
        #     cls.init_raspa_widom,
        #     cls.run_raspa_widom,
        #     if_(cls.should_run_gcmc),(
        #         cls.init_raspa_gcmc,
        #         while_(cls.should_run_another_gcmc)(
        #            cls.run_raspa_gcmc),
        #            ),
        #     # cls.run_raspa_gcmc,
        # #     cls.run_raspa_gcmc_low,
        # # )
        spec.outline(
            cls.setup,
            cls.run_pore_dia_zeopp,
            #cls.init_raspa_widom,
            #cls.run_raspa_widom,
            if_(cls.should_run_zeopp_volpo)(
                cls.run_savolpo_zeopp,
                 # if_(cls.should_run_zeopp_block)(
                    # cls.run_block_zeopp,
                 # ),
                cls.init_raspa_widom,
                cls.run_raspa_widom),
                if_(cls.should_run_gcmc)(
                    cls.init_raspa_gcmc,
                    while_(cls.should_run_another_gcmc)(
                            cls.run_raspa_gcmc,
                            cls.parse_raspa_gcmc,
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
            "max_wallclock_seconds": 1 * 30 * 60,
            "withmpi": False,
        }

        self.ctx.raspa_options = {
            "resources": {
                "num_machines": 1,
                "tot_num_mpiprocs": 1,
            },
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

    def should_run_zeopp_volpo(self):
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
            self.report("This is a suitable structure for further investigation")
            return True
        else:
            self.report("It does not look like promising: stop")
            return False

    def run_savolpo_zeopp(self):
        """
        It calculated the surface area and pore volume for pre-selected materials.
        """
        params = {
                'sa'   : [self.inputs.zeopp_probe_radius.value,
                          self.inputs.zeopp_probe_radius.value,
                          self.inputs.zeopp_block_samples_A3.value,
                          self.inputs.structure.label + '.sa'],
                'volpo': [self.inputs.zeopp_probe_radius.value,
                          self.inputs.zeopp_probe_radius.value,
                          self.inputs.zeopp_volpo_samples_UC.value,
                          self.inputs.structure.label + '.volpo'],
                'ha'   :  self.inputs.zeopp_accuracy.value
        }

        # Required inputs
        inputs = {
            'code'      : self.inputs.zeopp_code,
            'structure' : self.inputs.structure,
            'parameters': NetworkParameters(dict=params).store(),
            'metadata'  :{
                'options':  self.ctx.zeopp_options,
                'label'  : 'SaVolpo Calculation',
            }
        }

        try:
            inputs['atomic_radii'] = self.inputs.zeopp_atomic_radii
            self.report("Zeopp will use atomic radii from the .rad file")
        except:
            self.report("Zeopp will use default atomic radii")

        # Creating the calculation process and submit it
        sv = self.submit(ZeoppCalculation, **inputs)
        self.report("pk: {} | Running zeo++ block pocket calculation".format(sv.pk))
        return ToContext(zeopp_sv=sv)

    def should_run_zeopp_block(self):
        """
        Base on the outcomes of SaVolpo calculation here we decide if we need to
        calcualte block pockets or no.
        """

        ASA = self.ctx.zeopp_sv.get_outgoing().get_node_by_label('output_parameters').get_dict()['ASA_A^2']
        NASA = self.ctx.zeopp_sv.get_outgoing().get_node_by_label('output_parameters').get_dict()['NASA_A^2']

        if (NASA > 0.0) and (ASA > 0.0):
            self.report("We need to calcualte the block pocket file!")
            return True
        elif (NASA == 0.0) and (ASA > 0.0):
            self.report("All pores is accessible!")
            return False
        else:
            self.report("This structure is full of inaccessibility! : stop")
            pass

    def run_block_zeopp(self):
        """
        It calculates the block pockets only on structures with non-accessible void fraction.
        """

        self.ctx.block_filename = '_'.join((self.inputs.structure.label, self.inputs.raspa_comp1.value))

        params = {
                'block': [self.inputs.zeopp_probe_radius.value,
                          self.inputs.zeopp_block_samples_A3.value,
                          self.ctx.block_filename + '.block'],
                'ha'   : self.inputs.zeopp_accuracy.value
        }

        # Required inputs
        inputs = {
            'code'      : self.inputs.zeopp_code,
            'structure' : self.inputs.structure,
            'parameters': NetworkParameters(dict=params).store(),
            'metadata'  :{
                'options':  self.ctx.zeopp_options,
                'label'  : 'Pore Block Calculation',
            }
        }

        try:
            inputs['atomic_radii'] = self.inputs.zeopp_atomic_radii
            self.report("Zeopp will use atomic radii from the .rad file")
        except:
            self.report("Zeopp will use default atomic radii")

        # Creating the calculation process and submit it
        block = self.submit(ZeoppCalculation, **inputs)
        self.report("pk: {} | Running zeo++ block pocket calculation".format(block.pk))
        return ToContext(zeopp_block=block)

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
                self.ctx.raspa_parameters['Component'][comp_name] = self.ctx.raspa_parameters['Component'].pop(key)
                self.ctx.raspa_parameters["Component"][comp_name]["WidomProbability"] = 1.0

        ucs = multiply_unit_cell(self.inputs.structure, self.inputs.raspa_cutoff.value)
        self.ctx.raspa_parameters['GeneralSettings']['UnitCells'] = "{} {} {}".format(ucs[0], ucs[1], ucs[2])

        # Adding block files:
        # We check if

        # Turn on charges if requested
        # if self.inputs.raspa_usecharges:
        #     self.ctx.raspa_parameters['GeneralSettings']['ChargeMethod'] = "Ewald"
        #     self.ctx.raspa_parameters['GeneralSettings']['EwaldPrecision'] = 1e-6
        #     self.ctx.raspa_parameters['GeneralSettings']['UseChargesFromCIFFile'] = "yes"
        # else:
        #     self.ctx.raspa_parameters['GeneralSettings']['ChargeMethod'] = "None"

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
        # TODO: make the label of structure variable.
        inputs = {
            'framework'  : {self.inputs.structure.label : self.inputs.structure},
            'code'       : self.inputs.raspa_code,
            'parameters' : ParameterData(dict=self.ctx.raspa_parameters).store(),
            'metadata'  :{
                'options':  self.ctx.raspa_options,
                'label'  : 'RASPA Widom Calculation',
            }
        }

        # bp1_path = os.path.join(self.ctx.zeopp_block.outputs.retrieved._repository._get_base_folder().abspath, self.ctx.block_filename + '.block')
        # block_pocket_comp1 = SinglefileData(file=bp1_path)
        #
        # inputs['block_pocket'] ={}
        # inputs['block_pocket'][self.ctx.block_filename] = {}
        # inputs['block_pocket'][self.ctx.block_filename] = block_pocket_comp1

        # Check if there are blocking spheres (reading the header of the file) and use them for Raspa
        # Here, we iterate over compnents and find suitable one if it exists.
        # This section needs to be modified nicely.

        #with open(self.ctx.zeopp_block['block'].get_abs_path() + '/path/' + \
        #          self.ctx.zeopp_block['block'].get_folder_list()[0]) as f:
        #    self.ctx.number_blocking_spheres = int(f.readline().strip())
        #if self.ctx.number_blocking_spheres > 0:
        #    inputs['block_component_0'] = self.ctx.zeopp['block']
        #    self.report("Blocking spheres ({}) are present and used for Raspa".format(self.ctx.number_blocking_spheres))
        #else:
        #    self.report("No blocking spheres found")

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
        Kh_desired_comp = output['hkust1']['components'][self.ctx.raspa_comp.comp1.name]['henry_coefficient_average']
        for key, value in self.ctx.raspa_comp.items():
            if key in list(self.inputs.raspa_comp):
                comp_name = value.name
                if comp_name != self.ctx.raspa_comp.comp1.name:
                    Kh_comp = output['hkust1']['components'][comp_name]['henry_coefficient_average']
                    if Kh_comp < Kh_desired_comp:
                        self.report("Ideal selectivity larger than the threshold: compute binary loading")
                        return True
                    else:
                        self.report("Ideal selectivty either below threshold or inverted!: don't compute isotherm")
                        return False

        # TODO: using count and len() for more than two components.

        # kh = {}
        # kh.update({comp_name:(output['hkust1']['components'][comp_name]['henry_coefficient_average'])/})

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
                # It currently takes the greatest henry coefficient.
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
            'framework'  : {'hkust1':self.inputs.structure},
            #'hkust1'     : self.inputs.structure_raspa,
            'code'       : self.inputs.raspa_code,
            'parameters' : ParameterData(dict=self.ctx.raspa_parameters).store(),
            'metadata'  :{
                'options':  self.ctx.raspa_options,
                'label'  : 'RASPA GCMC Low Pressure',
            }
        }

        # Check if there are poket blocks to be loaded
        # Needs to be revised.
        #if self.ctx.number_blocking_spheres > 0:
        #    inputs['block_component_0'] = self.ctx.zeopp['block']


        # Create the calculation process and launch it
        gcmc = self.submit(RaspaCalculation, **inputs)
        self.report("pk: {} | Running Raspa GCMC at p(bar)={:.3f} ({} of {})".format(gcmc.pk, pressure/1e5, self.ctx.current_p_index+1, len(self.ctx.pressures)))
        return ToContext(raspa_gcmc=gcmc)

    def parse_raspa_gcmc(self):
        """
        Extract the pressure and loading average of the last completed raspa
        calculation.
        """
        pressure = self.ctx.raspa_parameters['System']['hkust1']['ExternalPressure']/1e5
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
                conv1 = output_gcmc['hkust1']['components'][comp_name]['conversion_factor_molec_uc_to_mol_kg']
                loading_average_comp = conv1 * output_gcmc['hkust1']['components'][comp_name]['loading_absolute_average']
                loading_dev_comp = conv1 * output_gcmc['hkust1']['components'][comp_name]['loading_absolute_dev']
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


        # Zeopp section
        output_zeo = self.ctx.zeopp_sv.outputs.output_parameters.get_dict()
        result_dict['Density'] = output_zeo['Density']
        result_dict['Density_unit'] = "g/cm^3"
        result_dict['POAV_Volume_fraction'] = output_zeo['POAV_Volume_fraction']
        result_dict['PONAV_Volume_fraction'] = output_zeo['PONAV_Volume_fraction']
        result_dict['POAV'] = output_zeo['POAV_cm^3/g']
        result_dict['POAV_unit'] = "cm^3/g"
        result_dict['GASA'] = output_zeo['ASA_m^2/g']
        result_dict['GASA_unit'] = "m^2/g"
        result_dict['VASA'] = output_zeo['ASA_m^2/cm^3']
        result_dict['VASA_unit'] = "m^2/cm^3"
        result_dict['NGASA'] = output_zeo['NASA_m^2/g']
        result_dict['NGASA_unit'] = "m^2/g"
        result_dict['NVASA'] = output_zeo['NASA_m^2/cm^3']
        result_dict['NVASA_unit'] = "ASA_m^2/cm^3"
        result_dict['Channel_surface_area'] = output_zeo['Channel_surface_area_A^2']
        result_dict['Channel_surface_area_unit'] = "A^2"
        result_dict['Pocket_surface_area'] = output_zeo['Pocket_surface_area_A^2']
        result_dict['Pocket_surface_area_unit'] = "A^2"


        #TODO needs to be fixed!
        # try:
        #     result_dict['number_blocking_spheres'] = self.ctx.number_blocking_spheres
        # except AttributeError:
        #     pass

        

        output_widom = self.ctx.raspa_widom.outputs.output_parameters.get_dict()
        result_dict['temperature'] = self.ctx.raspa_parameters["System"]['hkust1']["ExternalTemperature"]
        result_dict['temperature_unit'] = "K"
        result_dict['henry_coefficient_units'] = output_widom['hkust1']['components'][self.ctx.raspa_comp.comp1.name]['henry_coefficient_units']
        result_dict['adsorption_energy_units'] = output_widom['hkust1']['components'][self.ctx.raspa_comp.comp1.name]['adsorption_energy_widom_units']

        # result_dict['conversion_factor_molec_uc_to_cm3stp_cm3'] = output_gcmc['hkust1']['components'][self.inputs.raspa_comp1.value]['conversion_factor_molec_uc_to_cm3stp_cm3']
        # result_dict['conversion_factor_molec_uc_to_gr_gr'] = output_gcmc['hkust1']['components'][self.inputs.raspa_comp1.value]['conversion_factor_molec_uc_to_gr_gr']
        # result_dict['conversion_factor_molec_uc_to_mol_kg'] = output_gcmc['hkust1']['components'][self.inputs.raspa_comp1.value]['conversion_factor_molec_uc_to_mol_kg']


        for key, value in self.ctx.raspa_comp.items():
            if key in list(self.inputs.raspa_comp):
                comp_name = value.name
                # Widom results
                comp = output_widom['hkust1']['components'][comp_name]
                result_dict['henry_coefficient_average'][comp_name] = comp['henry_coefficient_average']
                result_dict['henry_coefficient_dev'][comp_name] = comp['henry_coefficient_dev']
                result_dict['adsorption_energy_average'][comp_name] = comp['adsorption_energy_widom_average']
                result_dict['adsorption_energy_dev'][comp_name] = comp['adsorption_energy_widom_dev']
                result_dict['isotherm_loading'][comp_name] = self.ctx.loading[comp_name]

        self.out("results", ParameterData(dict=result_dict).store())
        # self.out('blocking_spheres', self.ctx.zeopp_block['block'])
        self.report("Workchain completed successfully")
        return
# EOF
