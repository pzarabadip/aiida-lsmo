"""
Multi Component Isotherm WorkChain
"""
from __future__ import absolute_import
import os
import six
from copy import deepcopy

from aiida.common import AttributeDict
from aiida.plugins import CalculationFactory, DataFactory
from aiida.orm import Code, Dict, Float, Int, List, Str, load_node
from aiida.engine import submit
from aiida.engine import ToContext, WorkChain, workfunction, if_, while_, append_


from aiida_lsmo_workflows.utils.multiply_unitcell import multiply_unit_cell
from aiida_lsmo_workflows.utils.pressure_points import choose_pressure_points

CifData = DataFactory("cif")
ParameterData = DataFactory("dict")
SinglefileData = DataFactory("singlefile")
FolderData = DataFactory('folder')
RaspaCalculation = CalculationFactory("raspa")
ZeoppCalculation = CalculationFactory("zeopp.network")
NetworkParameters = DataFactory("zeopp.parameters")

# Lambda function taken from (https://stackoverflow.com/a/36977549)
# to make report nicer by using ordinary numbers.
ordinal = lambda n: "%d%s"%(n,{1:"st",2:"nd",3:"rd"}.get(n if n<20 else n%10,"th"))

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
        spec.input("structure", valid_type=CifData, required=True, help="Input structure in cif format")
        spec.input("general_calc_params", valid_type=Dict, required=False)

        # Zeopp inputs
        spec.input("zeopp_code", valid_type=Code, required=False)
        spec.input("zeopp_atomic_radii", valid_type=SinglefileData, required=False)

        # Raspa inputs
        spec.input("raspa_code", valid_type=Code, required=False)
        spec.input("raspa_parameters", valid_type=ParameterData, required=False)
        spec.input_namespace("raspa_comp", valid_type=dict,required=False, dynamic=True)
        spec.input('retrieved_parent_folder', valid_type=FolderData, required=False)

        # Scheduler options.
        spec.input_namespace("zeopp_options", required=False, dynamic=True, non_db=True)
        spec.input_namespace("raspa_options", required=False, dynamic=True, non_db=True)

        # Workflow
        spec.outline(
            cls.setup,
            cls.run_pore_dia_zeopp,
            if_(cls.should_run_zeopp_full)(
                cls.run_zeopp_full,
                cls.inspect_zeopp_calc,
                if_(cls.should_run_widom)(
                    cls.init_raspa_widom,
                    cls.run_raspa_widom),
                    if_(cls.should_run_gcmc)(
                        cls.init_raspa_gcmc,
                        while_(cls.should_run_another_gcmc)(
                            while_(cls.gcmc_conv_check)(
                                cls.run_raspa_gcmc,
                            ),
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
        # Set the number of runs counter to zero.
        self.ctx.nruns = 0

        # Create AttributeDict for RASPA components dictionary.
        self.ctx.raspa_comp = AttributeDict(self.inputs.raspa_comp)

        # Create context of general calculations parameters.
        self.ctx.general_calc_params = self.inputs.general_calc_params

        # Scheduler options
        self.ctx.zeopp_options = self.inputs.zeopp_options
        self.ctx.raspa_options = self.inputs.raspa_options

    def run_pore_dia_zeopp(self):
        """
        It performs the zeopp pore diameter calculation.
        """
        # Required inputs
        ha_flag = self.ctx.general_calc_params["zeopp"]["accuracy"]
        params = {
            "res" : [self.inputs.structure.label + ".res"],
            "ha"  : ha_flag
        }

        parameters = NetworkParameters(dict=params)

        inputs = {
            "code"      : self.inputs.zeopp_code,
            "structure" : self.inputs.structure,
            "parameters": parameters,
            "metadata"  :{
                "options":  self.ctx.zeopp_options,
                "label"  : "Pore Diameter Calculation",
                "description" : "Pore Diameter Calculation for <{}> at <{}> accuracy level".format(self.inputs.structure.label,ha_flag)
            }
        }

        # Use default zeopp atomic radii only if a .rad file is not specified
        try:
            inputs["atomic_radii"] = self.inputs.zeopp_atomic_radii
            self.report("Zeopp will use atomic radii from the .rad file")
        except:
            self.report("Zeopp will use default atomic radii")

        # Creating the calculation process and submit it
        res = self.submit(ZeoppCalculation, **inputs)
        self.report("pk: <{}> | Running zeo++ pore diameter calculation".format(res.pk))
        return ToContext(zeopp_res=res)

    def should_run_zeopp_full(self):
        """
        It uses largest included sphere (Di or LCD) and largest free sphere
        (Df or PLD) as pre-screenig descriptors to pass or reject the
        structure.
        """

        lcd_lim = self.ctx.general_calc_params["zeopp"]["lcd_max"]
        pld_lim = self.ctx.general_calc_params["zeopp"]["pld_min"]
        lcd_current = self.ctx.zeopp_res.outputs.output_parameters.get_dict()["Largest_included_sphere"]
        pld_current = self.ctx.zeopp_res.outputs.output_parameters.get_dict()["Largest_free_sphere"]

        if (lcd_current < lcd_lim) and (pld_current > pld_lim):
            self.report("<{}> is a suitable structure for further investigation".format(self.inputs.structure.label))
            return True
        else:
            self.report("<{}> does not look like promising: stop".format(self.inputs.structure.label))
            return False

    def run_zeopp_full(self):
        """
        It calculated the surface area, pore volume, and
        possible block pockets for pre-selected materials.
        """
        ha_flag = self.ctx.general_calc_params["zeopp"]["accuracy"]
        # Iterating over provided components.
        for key, value in self.ctx.raspa_comp.items():
            if key in list(self.inputs.raspa_comp):
                comp_name = value.name
                probe_radius = value.radius

                params = {
                        "sa"   : [probe_radius,
                                  probe_radius,
                                  self.ctx.general_calc_params["zeopp"]["sa_samples"],
                                  self.inputs.structure.label + "_" + comp_name + ".sa"],
                        "volpo": [probe_radius,
                                  probe_radius,
                                  self.ctx.general_calc_params["zeopp"]["volpo_samples"],
                                  self.inputs.structure.label + "_" + comp_name + ".volpo"],
                        "block": [probe_radius,
                                  self.ctx.general_calc_params["zeopp"]["block_samples"],
                                  self.inputs.structure.label + "_" + comp_name + ".block"],
                        "ha"   :  ha_flag
                }

                # Required inputs
                inputs = {
                    "code"      : self.inputs.zeopp_code,
                    "structure" : self.inputs.structure,
                    "parameters": NetworkParameters(dict=params).store(),
                    "metadata"  :{
                        "options":  self.ctx.zeopp_options,
                        "label"  : "Zeo++ SaVolpoBlock Calculation",
                        "description" : "Zeo++ SaVolpoBlock Calculation for <{}> using <{}_{}> probe at <{}> accuracy level".format(self.inputs.structure.label,comp_name,probe_radius,ha_flag)
                    }
                }

                try:
                    inputs["atomic_radii"] = self.inputs.zeopp_atomic_radii
                    self.report("Zeopp will use atomic radii from the .rad file")
                except:
                    self.report("Zeopp will use default atomic radii")

                # Creating the calculation process and submit it
                zeopp_full = self.submit(ZeoppCalculation, **inputs)
                zeopp_label = "zeopp_{}".format(comp_name)
                self.report("pk: <{}> | Running Zeo++ SaVolpoBlock calculation using <{}> probe".format(zeopp_full.pk,comp_name))
                self.to_context(**{zeopp_label: zeopp_full})

    def inspect_zeopp_calc(self):
        for key, value in self.ctx.raspa_comp.items():
            if key in list(self.inputs.raspa_comp):
                comp_name = value.name
                zeopp_label = "zeopp_{}".format(comp_name)
                self.report("Checking if <{}> job finished successfully".format(zeopp_label))
                assert self.ctx[zeopp_label].is_finished_ok

    def should_run_widom(self):
        """
        Submit widom calculation only if there is some accessible volume for
        desired component.
        """
        zeopp_label = "zeopp_{}".format(self.ctx.raspa_comp.comp1.name)
        output = self.ctx[zeopp_label].outputs.output_parameters.get_dict()
        POAV_desired_comp = output["POAV_Volume_fraction"]

        if POAV_desired_comp > 1e-5:
            self.report("Found accessible pore volume for component of interest: continue")
            return True
        else:
            self.report("NOT Found any accessible pore volume: stop")
            return False

    def init_raspa_widom(self):
        """
        It generates the ParameterData for RASPA and use the blocking spheres if there is any.
        With current plugin we can have more than one component to be used for Widom insertion.
        """

        # Create a deepcopy of the user parameters, to modify before submission
        self.ctx.raspa_parameters = deepcopy(self.inputs.raspa_parameters.get_dict())
        # Create key for the framework in the System subdict.
        # TODO: Improve it to handle more than one system?
        self.ctx.raspa_parameters["System"][self.inputs.structure.label] = {}

        # Iterating over components and update the input component section.
        for key, value in self.ctx.raspa_comp.items():
            if key in list(self.inputs.raspa_comp):
                comp_name = value.name
                mol_def = value.mol_def
                bp_label = "_".join((self.inputs.structure.label,comp_name))
                self.ctx.raspa_parameters["Component"][comp_name] = self.ctx.raspa_parameters["Component"].pop(key)
                self.ctx.raspa_parameters["Component"][comp_name]["MoleculeDefinition"] = mol_def
                self.ctx.raspa_parameters["Component"][comp_name]["WidomProbability"] = 1.0

        # Obtain the unit cell multiplications and update the corresponding input section.
        ucs = multiply_unit_cell(self.inputs.structure, self.ctx.general_calc_params["raspa"]["cutoff"])
        self.ctx.raspa_parameters["System"][self.inputs.structure.label]["UnitCells"] = "{} {} {}".format(ucs[0], ucs[1], ucs[2])

        # Turn on charges if requested
        if self.ctx.general_calc_params["raspa"]["usecharges"]:
            self.ctx.raspa_parameters["GeneralSettings"]["ChargeMethod"] = "Ewald"
            self.ctx.raspa_parameters["GeneralSettings"]["EwaldPrecision"] = 1e-6
            # Setting use charge from cif or EQeq
            if self.ctx.general_calc_params["raspa"]["charge_from_cif"]:
                self.ctx.raspa_parameters["GeneralSettings"]["UseChargesFromCIFFile"] = "yes"
            else:
                self.ctx.raspa_parameters["GeneralSettings"]["ChargeFromChargeEquilibration"] = "yes"
        else:
            self.ctx.raspa_parameters["GeneralSettings"]["ChargeMethod"] = "None"

        # Updating input for Widom particle insertion.
        # TODO: Improve number of cycles setting.
        self.ctx.raspa_parameters["System"][self.inputs.structure.label]["type"] = self.ctx.general_calc_params["raspa"]["system_type"]
        self.ctx.raspa_parameters["System"][self.inputs.structure.label]["ExternalTemperature"] = self.ctx.general_calc_params["raspa"]["temperature"]
        self.ctx.raspa_parameters["GeneralSettings"]["SimulationType"] = self.ctx.general_calc_params["raspa"]["simulation_type"]
        self.ctx.raspa_parameters["GeneralSettings"]["CutOff"] = self.ctx.general_calc_params["raspa"]["cutoff"]
        self.ctx.raspa_parameters["GeneralSettings"]["NumberOfInitializationCycles"] = 0
        self.ctx.raspa_parameters["GeneralSettings"]["NumberOfCycles"] = self.ctx.general_calc_params["raspa"]["widom_cycle_mult"] * self.inputs.raspa_parameters.get_dict()["GeneralSettings"]["NumberOfCycles"]
        self.ctx.raspa_parameters["GeneralSettings"]["PrintPropertiesEvery"] = int(self.ctx.raspa_parameters["GeneralSettings"]["NumberOfCycles"] / self.ctx.general_calc_params["raspa"]["verbosity"])
        self.ctx.raspa_parameters["GeneralSettings"]["PrintEvery"] = int(1e6) #never
        return

    def run_raspa_widom(self):
        """
        It runs a Widom calculation in RASPA.
        """
        # Create the inputs dictionary
        inputs = {
            "framework"  : {self.inputs.structure.label : self.inputs.structure},
            "code"       : self.inputs.raspa_code,
            "metadata"  :{
                "options":  self.ctx.raspa_options,
                "label"  : "RASPA Widom Calculation",
                "description" : "RASPA Widom Calculation for <{}>".format(self.inputs.structure.label)
            }
        }

        # Creating dictionary for block_pocket
        inputs["block_pocket"] = {}
        self.ctx.number_blocking_spheres = {}

        # Iterating over components and update input for using block pockets.
        for key, value in self.ctx.raspa_comp.items():
            if key in list(self.inputs.raspa_comp):
                comp_name = value.name
                zeopp_label = "zeopp_{}".format(comp_name)
                bp_label = "_".join((self.inputs.structure.label,comp_name))
                bp_path = os.path.join(self.ctx[zeopp_label].outputs.retrieved._repository._get_base_folder().abspath, bp_label + ".block")

                with open(bp_path, "r") as block_file:
                    self.ctx.number_blocking_spheres[comp_name] = int(block_file.readline().strip())
                    if self.ctx.number_blocking_spheres[comp_name] > 0:
                        self.ctx.raspa_parameters["Component"][comp_name]["BlockPocketsFileName"] = {}
                        self.ctx.raspa_parameters["Component"][comp_name]["BlockPocketsFileName"][self.inputs.structure.label] = bp_label

                        inputs["block_pocket"][bp_label] = self.ctx[zeopp_label].outputs.block

                        self.report("<{}> Blocking spheres are present for <{}> and used for Raspa".format(self.ctx.number_blocking_spheres[comp_name],comp_name))
                    else:
                        self.report("No blocking spheres found for <{}>".format(comp_name))

        # Storing parameters, and submittig the process.
        inputs["parameters"] = ParameterData(dict=self.ctx.raspa_parameters).store()
        widom = self.submit(RaspaCalculation, **inputs)
        self.report("pk: <{}> | Running Raspa Widom insertion".format(widom.pk))
        return ToContext(raspa_widom=widom)

    def should_run_gcmc(self):
        """
        If the henry coefficient of main component would be higher
        than a pre-defined threshold, it will return True.
        """
        # Getting Kh of desired component and its user-defined threshold.
        output = self.ctx.raspa_widom.outputs.output_parameters.get_dict()
        self.ctx.Kh_desired_comp = output[self.inputs.structure.label]["components"][self.ctx.raspa_comp.comp1.name]["henry_coefficient_average"]
        Kh_min = self.ctx.general_calc_params["raspa"]["kh_min"]

        # TODO: Improving for checking against other components too.
        if self.ctx.Kh_desired_comp > Kh_min:
            self.report("Kh_{} <{}> is greater than <{}>: compute the GCMC".format(self.ctx.raspa_comp.comp1.name,self.ctx.Kh_desired_comp,Kh_min))
            return True
        else:
            self.report("Kh_{} <{}> is less than <{}>: DO NOT compute the GCMC".format(self.ctx.raspa_comp.comp1.name,self.ctx.Kh_desired_comp,Kh_min))
            return False

    def init_raspa_gcmc(self):
        """
        Initialization of input for raspa gcmc.
        It generates a pressure list based on our decision.
        1 - Full and Dyanmic isotherm Simulation:
        This part is based on Daniele work and generates the pressure points List
        using the provided input of pore volume and density.
        2 - Full and non-dynamic isotherm simulation:
        It generates the pressure list betweem with equal distance between a min and maximum pressure.
        3 - Selecte: It just takes user defined pressure list and does the rest.
        """
        # Initializate counter and set restart to None
        self.ctx.current_p_index = 0
        self.ctx.restart_raspa_calc = None

        # Case 1
        if (self.ctx.general_calc_params["raspa"]["isotherm_dynamic"]) and (self.ctx.general_calc_params["raspa"]["isotherm_full"]):
            # Estimate the total loading qsat and choose the pressure points
            satDens = self.ctx.general_calc_params["raspa"]["molsatdens"] #(mol/l)
            # TODO: Currently only support single-comp isotherm.
            poreVol = self.ctx["zeopp_" + self.ctx.raspa_comp.comp1.name].outputs.output_parameters.get_dict()["POAV_cm^3/g"] #(cm3/g = l/kg)
            self.ctx.estimated_qsat = satDens * poreVol
            self.ctx.pressures = choose_pressure_points(
                Kh = self.ctx.Kh_desired_comp, #(mol/kg/Pa)
                qsat = self.ctx.estimated_qsat, #(mol/kg_frame)
                dpa = self.ctx.general_calc_params["raspa"]["dpa"], #(kg*Pa/mol)
                dpmax = self.ctx.general_calc_params["raspa"]["dpmax"], #(Pa)
                pmax = self.ctx.general_calc_params["raspa"]["pressure_max"], #(Pa)
                pmin = self.ctx.general_calc_params["raspa"]["pressure_min"], #(Pa)
                dynamic = self.ctx.general_calc_params["raspa"]["isotherm_dynamic"],
                full = self.ctx.general_calc_params["raspa"]["isotherm_full"],
            )
        # Case 2
        elif (not self.ctx.general_calc_params["raspa"]["isotherm_dynamic"]) and (self.ctx.general_calc_params["raspa"]["isotherm_full"]):
            self.ctx.pressures = choose_pressure_points(
                pmin = self.ctx.general_calc_params["raspa"]["pressure_min"],
                dpa = self.ctx.general_calc_params["raspa"]["dpa"],
                pmax = self.ctx.general_calc_params["raspa"]["pressure_max"],
                dynamic = self.ctx.general_calc_params["raspa"]["isotherm_dynamic"],
                full = self.ctx.general_calc_params["raspa"]["isotherm_full"],
            )
        # Case 3
        else:
            self.ctx.pressures = choose_pressure_points(
                selected_pressures = self.ctx.general_calc_params["raspa"]["selected_pressures"],
                dynamic = self.ctx.general_calc_params["raspa"]["isotherm_dynamic"],
                full = self.ctx.general_calc_params["raspa"]["isotherm_full"],
            )

        self.report("<{}> points are chosen for GCMC calculation".format(len(self.ctx.pressures)))

        # CORRECT the parameters to perform GCMC
        self.ctx.raspa_parameters["GeneralSettings"]["NumberOfInitializationCycles"] = self.inputs.raspa_parameters.get_dict()["GeneralSettings"]["NumberOfInitializationCycles"]
        self.ctx.raspa_parameters["GeneralSettings"]["NumberOfCycles"] = self.inputs.raspa_parameters.get_dict()["GeneralSettings"]["NumberOfCycles"]
        self.ctx.raspa_parameters["GeneralSettings"]["PrintPropertiesEvery"] = int(1e6) #never
        self.ctx.raspa_parameters["GeneralSettings"]["PrintEvery"] = int(self.ctx.raspa_parameters["GeneralSettings"]["NumberOfCycles"]/self.ctx.general_calc_params["raspa"]["verbosity"])

        # Iterating over components and update the input.
        # TODO: Make the MC move probabilities user-defined.
        for key, value in self.ctx.raspa_comp.items():
            if key in list(self.inputs.raspa_comp):
                comp_name = value.name
                mol_frac = value.mol_fraction
                singlebead = value.singlebead
                self.ctx.raspa_parameters["Component"][comp_name]["MolFraction"] = float(mol_frac)
                del self.ctx.raspa_parameters["Component"][comp_name]["WidomProbability"]
                self.ctx.raspa_parameters["Component"][comp_name]["TranslationProbability"] = 0.5
                # Only adds RotationProbability move if it is not singlebead model.
                if not singlebead:
                    self.ctx.raspa_parameters["Component"][comp_name]["RotationProbability"] = 0.5
                self.ctx.raspa_parameters["Component"][comp_name]["ReinsertionProbability"] = 0.5
                self.ctx.raspa_parameters["Component"][comp_name]["SwapProbability"] = 1.0
                self.ctx.raspa_parameters["Component"][comp_name]["IdentityChangeProbability"] = 1.0
                self.ctx.raspa_parameters["Component"][comp_name]["NumberOfIdentityChanges"] = len(list(self.inputs.raspa_comp))
                self.ctx.raspa_parameters["Component"][comp_name]["IdentityChangesList"] = [i for i in range(len(list(self.inputs.raspa_comp)))]
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
        # Updateing the ExternalPressure.
        pressure = self.ctx.pressures[self.ctx.current_p_index]
        self.ctx.raspa_parameters["System"][self.inputs.structure.label]["ExternalPressure"] = pressure
        # Create the input dictionary
        inputs = {
            "framework"  : {self.inputs.structure.label:self.inputs.structure},
            "code"       : self.inputs.raspa_code,
            "parameters" : ParameterData(dict=self.ctx.raspa_parameters).store(),
            "metadata"  :{
                "options":  self.ctx.raspa_options,
                "label"  : "RASPA GCMC Calculation",
                "description" : "<{}> RASPA GCMC Calculation at p(bar)={:.3f}".format(ordinal(self.ctx.nruns+1),pressure/1e5)
            }
        }

        # Iterating over components and update the block pocket section.
        inputs["block_pocket"] ={}
        for key, value in self.ctx.raspa_comp.items():
            if key in list(self.inputs.raspa_comp):
                comp_name = value.name
                zeopp_label = "zeopp_{}".format(comp_name)
                bp_label = "_".join((self.inputs.structure.label,comp_name))
                if self.ctx.number_blocking_spheres[comp_name] > 0:
                    inputs["block_pocket"][bp_label] = self.ctx[zeopp_label].outputs.block

        # Handling the retrive and usage of restart feature.
        if self.ctx.restart_raspa_calc is not None:
           inputs['retrieved_parent_folder'] = self.ctx.restart_raspa_calc

        # Create the calculation process and launch it
        gcmc = self.submit(RaspaCalculation, **inputs)
        self.report("pk: <{}> | Running Raspa GCMC at p(bar)={:.3f} <{} of {}>".format(gcmc.pk, pressure/1e5, self.ctx.current_p_index+1, len(self.ctx.pressures)))
        # Increase the number of runs counter.
        self.ctx.nruns += 1
        return ToContext(raspa_gcmc=gcmc)

    def gcmc_conv_check(self):
        """
        It checks if the GCMC calculation is converged to a desired level
        of accuracy. The threshold should be passed as general_calc_params.
        What do we need?
        """

        # Workaround for the first calculation at each pressure.
        if self.ctx.nruns == 0:
            converged = [False]
        # Worksround for the 2nd and more calculation.
        else:
            # Iterating over components and check the convergences.
            converged = []
            output_gcmc = self.ctx.raspa_gcmc.outputs.output_parameters.get_dict()
            self.ctx.restart_raspa_calc = self.ctx.raspa_gcmc.outputs['retrieved']
            for key, value in self.ctx.raspa_comp.items():
                if key in list(self.inputs.raspa_comp):
                    comp_name = value.name
                    conv_threshold = value.conv_threshold
                    conv1 = output_gcmc[self.inputs.structure.label]["components"][comp_name]["conversion_factor_molec_uc_to_mol_kg"]
                    loading_average_comp = conv1 * output_gcmc[self.inputs.structure.label]["components"][comp_name]["loading_absolute_average"]
                    loading_dev_comp = conv1 * output_gcmc[self.inputs.structure.label]["components"][comp_name]["loading_absolute_dev"]

                    # TODO: It can happen at very low pressure.
                    # Solution: We can increase the pressure index if it happens.
                    if loading_average_comp == 0:
                        converged.append(False)
                        self.report("<{}> is not adsorbed at current simulation setup".format(comp_name))
                    else:
                        error = round((loading_dev_comp / loading_average_comp), 2)
                        if error <= conv_threshold:
                            self.report("Error(%) <{:.1f}> less than <{:.1f}> | GCMC simulation is converged for <{}>".format(error*100, conv_threshold*100, comp_name))
                            converged.append(True)
                        else:
                            self.report("Error(%) <{:.1f}> greater than <{:.1f}> | GCMC simulation IS NOT converged for <{}>".format(error*100, conv_threshold*100, comp_name))
                            converged.append(False)


        # Get pressure for having in the self.report
        pressure = self.ctx.pressures[self.ctx.current_p_index]
        additional_cycle = self.ctx.general_calc_params['raspa']['additional_cycle']
        # First run.
        if self.ctx.nruns == 0:
            self.ctx.raspa_parameters["GeneralSettings"]["NumberOfInitializationCycles"] = self.ctx.raspa_parameters["GeneralSettings"]["NumberOfInitializationCycles"]
            self.ctx.raspa_parameters["GeneralSettings"]["NumberOfCycles"] = self.ctx.raspa_parameters["GeneralSettings"]["NumberOfCycles"]
            self.report("<{}> GCMC simulation at p(bar)=<{:.3}>".format(ordinal(self.ctx.nruns+1),pressure/1e5))
        # Second and more runs which are not converged.
        elif (not all(converged)) and self.ctx.nruns != 0:
            self.ctx.raspa_parameters["GeneralSettings"]["NumberOfInitializationCycles"] = self.ctx.raspa_parameters["GeneralSettings"]["NumberOfInitializationCycles"] + additional_cycle
            self.ctx.raspa_parameters["GeneralSettings"]["NumberOfCycles"] = self.ctx.raspa_parameters["GeneralSettings"]["NumberOfCycles"] + additional_cycle
            ParameterData(dict=self.ctx.raspa_parameters).store()
            new_init = int(self.ctx.raspa_parameters["GeneralSettings"]["NumberOfInitializationCycles"])
            new_prod = int(self.ctx.raspa_parameters["GeneralSettings"]["NumberOfCycles"])
            self.report("<{}> GCMC simulation at p(bar)=<{:.3}>".format(ordinal(self.ctx.nruns+1),pressure/1e5))
            self.report("Increasing MC cycles to Initialization<{:d}> -- Production<{:d}>".format(new_init,new_prod))
        # Converged run.
        else:
            self.report("GCMC simulation at p(bar)=<{:.3}> is converged after <{}> run".format(pressure/1e5,ordinal(self.ctx.nruns)))
            self.report("Resetting MC cycles and run counter")
            self.ctx.nruns = 0
            init = int(self.inputs.raspa_parameters.get_dict()["GeneralSettings"]["NumberOfInitializationCycles"])
            prod = int(self.inputs.raspa_parameters.get_dict()["GeneralSettings"]["NumberOfCycles"])
            self.ctx.raspa_parameters["GeneralSettings"]["NumberOfInitializationCycles"] = init
            self.ctx.raspa_parameters["GeneralSettings"]["NumberOfCycles"] = prod
            ParameterData(dict=self.ctx.raspa_parameters).store()

        return not all(converged)

    def parse_raspa_gcmc(self):
        """
        Extract the pressure and loading average of the last completed and converged
        RASPA calculation.
        """
        pressure = self.ctx.raspa_parameters["System"][self.inputs.structure.label]["ExternalPressure"]/1e5
        output_gcmc = self.ctx.raspa_gcmc.outputs.output_parameters.get_dict()

        # Creating the loading empty dictionary only at first run.
        if self.ctx.current_p_index == 0:
            self.ctx.loading = {}
            for key, value in self.ctx.raspa_comp.items():
                if key in list(self.inputs.raspa_comp):
                    comp_name = value.name
                    self.ctx.loading[comp_name] = []

        # Iterating over components and append the loadings and
        # error bars to the dictionary.
        for key, value in self.ctx.raspa_comp.items():
            if key in list(self.inputs.raspa_comp):
                comp_name = value.name
                conv1 = output_gcmc[self.inputs.structure.label]["components"][comp_name]["conversion_factor_molec_uc_to_mol_kg"]
                loading_average_comp = conv1 * output_gcmc[self.inputs.structure.label]["components"][comp_name]["loading_absolute_average"]
                loading_dev_comp = conv1 * output_gcmc[self.inputs.structure.label]["components"][comp_name]["loading_absolute_dev"]
                self.ctx.loading[comp_name].append([pressure, loading_average_comp, loading_dev_comp])

        # Simulation is converged that we are parsing. Increase the pressure index.
        self.ctx.current_p_index += 1
        return

    def return_results(self):
        """
        Attach the results to the output.
        """
        # Create empty results dictionary.
        result_dict = {}

        # Zeopp section
        try:
            # Output of Di and Df
            result_dict["Largest_included_sphere"] = self.ctx.zeopp_res.outputs.output_parameters.get_dict()["Largest_included_sphere"]
            result_dict["Largest_free_sphere"] = self.ctx.zeopp_res.outputs.output_parameters.get_dict()["Largest_free_sphere"]
            # Creating the needed keys within the result dictionary for zeo++ results.
            result_dict["POAV_Volume_fraction"] = {}
            result_dict["PONAV_Volume_fraction"] = {}
            result_dict["POAV"] = {}
            result_dict["GASA"] = {}
            result_dict["VASA"] = {}
            result_dict["NGASA"] = {}
            result_dict["NVASA"] = {}
            result_dict["Channel_surface_area"] = {}
            result_dict["Pocket_surface_area"] = {}
            result_dict["Density_unit"] = "g/cm^3"
            result_dict["POAV_unit"] = "cm^3/g"
            result_dict["GASA_unit"] = "m^2/g"
            result_dict["VASA_unit"] = "m^2/cm^3"
            result_dict["NGASA_unit"] = "m^2/g"
            result_dict["NVASA_unit"] = "ASA_m^2/cm^3"
            result_dict["Channel_surface_area_unit"] = "A^2"
            result_dict["Pocket_surface_area_unit"] = "A^2"
            result_dict["number_blocking_spheres"] = {}
            # Iterating over components and extract results.
            for key, value in self.ctx.raspa_comp.items():
                if key in list(self.inputs.raspa_comp):
                    comp_name = value.name
                    zeopp_label = "zeopp_{}".format(comp_name)
                    output_zeo = self.ctx[zeopp_label].outputs.output_parameters.get_dict()
                    result_dict["POAV_Volume_fraction"][comp_name] = output_zeo["POAV_Volume_fraction"]
                    result_dict["PONAV_Volume_fraction"][comp_name] = output_zeo["PONAV_Volume_fraction"]
                    result_dict["POAV"][comp_name] = output_zeo["POAV_cm^3/g"]
                    result_dict["GASA"][comp_name] = output_zeo["ASA_m^2/g"]
                    result_dict["VASA"][comp_name] = output_zeo["ASA_m^2/cm^3"]
                    result_dict["NGASA"][comp_name] = output_zeo["NASA_m^2/g"]
                    result_dict["NVASA"][comp_name] = output_zeo["NASA_m^2/cm^3"]
                    result_dict["Channel_surface_area"][comp_name] = output_zeo["Channel_surface_area_A^2"]
                    result_dict["Pocket_surface_area"][comp_name] = output_zeo["Pocket_surface_area_A^2"]
                    result_dict["number_blocking_spheres"][comp_name] = self.ctx.number_blocking_spheres[comp_name]
                    result_dict["Density"] = output_zeo["Density"]
        except AttributeError:
            pass

        # RASPA widom Section
        try:
            # Getting the output parameters of widom calculation.
            output_widom = self.ctx.raspa_widom.outputs.output_parameters.get_dict()
            # Creating the needed keys within the result dictionary for RASPA Widom results.
            result_dict["henry_coefficient_average"] = {}
            result_dict["henry_coefficient_dev"] = {}
            result_dict["adsorption_energy_average"] = {}
            result_dict["adsorption_energy_dev"] = {}
            result_dict["temperature"] = self.ctx.raspa_parameters["System"][self.inputs.structure.label]["ExternalTemperature"]
            result_dict["temperature_unit"] = "K"
            result_dict["henry_coefficient_units"] = output_widom[self.inputs.structure.label]["components"][self.ctx.raspa_comp.comp1.name]["henry_coefficient_units"]
            result_dict["adsorption_energy_units"] = output_widom[self.inputs.structure.label]["components"][self.ctx.raspa_comp.comp1.name]["adsorption_energy_widom_units"]
            # Iterating over components and extract results.
            for key, value in self.ctx.raspa_comp.items():
                if key in list(self.inputs.raspa_comp):
                    comp_name = value.name
                    comp = output_widom[self.inputs.structure.label]["components"][comp_name]
                    mol_frac = value.mol_fraction
                    result_dict["henry_coefficient_average"][comp_name] = comp["henry_coefficient_average"]
                    result_dict["henry_coefficient_dev"][comp_name] = comp["henry_coefficient_dev"]
                    result_dict["adsorption_energy_average"][comp_name] = comp["adsorption_energy_widom_average"]
                    result_dict["adsorption_energy_dev"][comp_name] = comp["adsorption_energy_widom_dev"]
        except AttributeDict:
            pass

        # RASPA Section
        try:
            # Getting the output parameters of converged GCMC calculation.
            output_gcmc = self.ctx.raspa_gcmc.outputs.output_parameters.get_dict()
            # Creating the needed keys within the result dictionary for RASPA Widom results.
            result_dict["conversion_factor_molec_uc_to_cm3stp_cm3"] = {}
            result_dict["conversion_factor_molec_uc_to_gr_gr"] = {}
            result_dict["conversion_factor_molec_uc_to_mol_kg"] = {}
            result_dict["loading"] = {}
            result_dict["isotherm_loading_header"] = ["Pressure(bar)", "Loading_average(mol/kg)", "Loading_deviation(mol/kg)"]
            result_dict["isotherm_loading"] = {}
            result_dict["mol_fraction"] = {}
            result_dict["isotherm_enthalpy"] = {}
            # Iterating over components and extract results.
            for key, value in self.ctx.raspa_comp.items():
                if key in list(self.inputs.raspa_comp):
                    comp_name = value.name
                    mol_frac = value.mol_fraction
                    result_dict["conversion_factor_molec_uc_to_cm3stp_cm3"][comp_name] = output_gcmc[self.inputs.structure.label]["components"][comp_name]["conversion_factor_molec_uc_to_cm3stp_cm3"]
                    result_dict["conversion_factor_molec_uc_to_gr_gr"][comp_name] = output_gcmc[self.inputs.structure.label]["components"][comp_name]["conversion_factor_molec_uc_to_gr_gr"]
                    result_dict["conversion_factor_molec_uc_to_mol_kg"][comp_name] = output_gcmc[self.inputs.structure.label]["components"][comp_name]["conversion_factor_molec_uc_to_mol_kg"]
                    result_dict["mol_fraction"][comp_name] = output_gcmc[self.inputs.structure.label]["components"][comp_name]["mol_fraction"]
                    result_dict["isotherm_loading"][comp_name] = self.ctx.loading[comp_name]
        except AttributeError:
            pass

        # Blocking Spheres Section
        for key, value in self.ctx.raspa_comp.items():
            if key in list(self.inputs.raspa_comp):
                comp_name = value.name
                zeopp_label = "zeopp_{}".format(comp_name)
                bp_label = "_".join((self.inputs.structure.label,comp_name))
                # Only keeping non-empty block pocket files.
                if self.ctx.number_blocking_spheres[comp_name] > 0:
                    self.out("blocking_spheres_{}".format(comp_name), self.ctx[zeopp_label].outputs.block)

        # Finalizing the results and report!
        self.out("results", ParameterData(dict=result_dict).store())
        self.report("Workchain completed successfully! | Result Dict is <{}>".format(self.outputs["results"].pk))
        return
# EOF
