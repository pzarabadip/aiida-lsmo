"""
GEMC WorkChain
"""
from __future__ import absolute_import
import os
import six
from copy import deepcopy

from aiida.common import AttributeDict
from aiida.plugins import CalculationFactory, DataFactory
from aiida.orm import Code, Dict
from aiida.engine import submit
from aiida.engine import ToContext, WorkChain, if_, while_, append_


from aiida_lsmo_workflows.utils.multiply_unitcell import multiply_unit_cell
from aiida_lsmo_workflows.utils.temperature_points import choose_temp_points

ParameterData = DataFactory("dict")
FolderData = DataFactory('folder')
RaspaCalculation = CalculationFactory("raspa")

# Lambda function taken from (https://stackoverflow.com/a/36977549)
# to make report nicer by using ordinary numbers.
ordinal = lambda n: "%d%s"%(n,{1:"st",2:"nd",3:"rd"}.get(n if n<20 else n%10,"th"))

class GEMCWorkChain(WorkChain):
    """
    The GEMC Worchain is designed to construct VLCC.
    """

    @classmethod
    def define(cls, spec):
        super(GEMCWorkChain, cls).define(spec)
        """
        Define workflow specification.
        This is the most important method of a Workchain, which defines the
        inputs it takes, the logic of the execution and the outputs
        that are generated in the process.
        """
        # General inputs
        spec.input("general_calc_params", valid_type=Dict, required=False)

        # Raspa inputs
        spec.input("raspa_code", valid_type=Code, required=False)
        spec.input("raspa_parameters", valid_type=ParameterData, required=False)
        spec.input_namespace("raspa_comp", valid_type=dict,required=False, dynamic=True)
        spec.input_namespace("raspa_box", valid_type=dict,required=False, dynamic=True)
        spec.input('retrieved_parent_folder', valid_type=FolderData, required=False)

        # Scheduler options.
        spec.input_namespace("raspa_options", required=False, dynamic=True, non_db=True)

        # Workflow
        spec.outline(
            cls.setup,
            cls.init_raspa_gemc,
            while_(cls.should_run_another_gemc)(
                while_(cls.gemc_conv_check)(
                        while_(cls.is_final_box_ok)(
                                cls.run_raspa_gemc
                        ),
                ),
                cls.parse_raspa_gemc,
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
        self.ctx.done = False
        # Create AttributeDict for RASPA components dictionary.
        self.ctx.raspa_comp = AttributeDict(self.inputs.raspa_comp)
        self.ctx.raspa_box = AttributeDict(self.inputs.raspa_box)

        # Create context of general calculations parameters.
        self.ctx.general_calc_params = self.inputs.general_calc_params

        # Setup initial the temperature list and initial temperature and index.
        self.ctx.temperature = choose_temp_points(
                                self.ctx.general_calc_params['raspa']['T_min'],
                                self.ctx.general_calc_params['raspa']['T_max'],
                                self.ctx.general_calc_params['raspa']['dT'],
                                )
        self.ctx.init_T_index = self.ctx.temperature[0].index(self.ctx.temperature[1])
        self.ctx.max_T_index = self.ctx.temperature[0].index(self.ctx.temperature[0][-1])
        self.ctx.min_T_index = self.ctx.temperature[0].index(self.ctx.temperature[0][0])

        self.ctx.uphill = True
        self.ctx.downhill = False

        self.report("<{}> number of temperature points are chosen for GEMC".format(len(self.ctx.temperature[0])))
        self.report("Starting from <{}>K, then toward maximum <{}>, and finally minimum <{}>".format(self.ctx.temperature[1],self.ctx.temperature[0][-1],self.ctx.temperature[0][0]))
        # Scheduler options
        self.ctx.raspa_options = self.inputs.raspa_options

    def init_raspa_gemc(self):
        """
        Initializing RASPA parameters for GEMC.
        """
        self.ctx.current_T_index = self.ctx.init_T_index
        # Create a deepcopy of the user parameters, to modify before submission
        self.ctx.raspa_parameters = deepcopy(self.inputs.raspa_parameters.get_dict())

        # Dealing with System section:
        for key, value in self.ctx.raspa_box.items():
            box_tag = value.tag
            box_ax = value.box_ax
            box_by = value.box_by
            box_cz = value.box_cz
            box_alpha = value.box_alpha
            box_beta = value.box_beta
            box_gamma = value.box_gamma
            # box_temperature = value.temperature
            self.ctx.raspa_parameters["System"][box_tag] = self.ctx.raspa_parameters["System"].pop(key)
            self.ctx.raspa_parameters["System"][box_tag]["type"] = "Box"
            self.ctx.raspa_parameters["System"][box_tag]["BoxLengths"] = "{} {} {}".format(box_ax, box_by, box_cz)
            self.ctx.raspa_parameters["System"][box_tag]["BoxAngles"] = "{} {} {}".format(box_alpha, box_beta, box_gamma)
            # self.ctx.raspa_parameters["System"][box_tag]["ExternalTemperature"] = self.ctx.temperature[1]

        # Dealing with component sections.
        for key, value in self.ctx.raspa_comp.items():
            comp_name = value.name
            mol_def = value.mol_def
            singlebead = value.singlebead
            number_box_one = value.box_one
            number_box_two = value.box_two
            self.ctx.raspa_parameters["Component"][comp_name] = self.ctx.raspa_parameters["Component"].pop(key)
            self.ctx.raspa_parameters["Component"][comp_name]["MoleculeDefinition"] = mol_def
            self.ctx.raspa_parameters["Component"][comp_name]["GibbsSwapProbability"] = 1.0
            self.ctx.raspa_parameters["Component"][comp_name]["CreateNumberOfMolecules"] = {}
            self.ctx.raspa_parameters["Component"][comp_name]["TranslationProbability"] = 0.5
            if not singlebead:
                self.ctx.raspa_parameters["Component"][comp_name]["RotationProbability"] = 0.5
            # self.ctx.raspa_parameters["Component"][comp_name]["ReinsertionProbability"] = 0.5
            self.ctx.raspa_parameters["Component"][comp_name]["CreateNumberOfMolecules"]["box_one"] = number_box_one
            self.ctx.raspa_parameters["Component"][comp_name]["CreateNumberOfMolecules"]["box_two"] = number_box_two

        # Turn on charges if requested
        if self.ctx.general_calc_params["raspa"]["usecharges"]:
            self.ctx.raspa_parameters["GeneralSettings"]["ChargeMethod"] = "Ewald"
            self.ctx.raspa_parameters["GeneralSettings"]["EwaldPrecision"] = 1e-6
        else:
            self.ctx.raspa_parameters["GeneralSettings"]["ChargeMethod"] = "None"

        # Dealing with GeneralSettings section.
        self.ctx.raspa_parameters["GeneralSettings"]["SimulationType"] = self.ctx.general_calc_params["raspa"]["simulation_type"]
        self.ctx.raspa_parameters["GeneralSettings"]["GibbsVolumeChangeProbability"] = self.ctx.general_calc_params["raspa"]["gibbs_volume_prop"]
        self.ctx.raspa_parameters["GeneralSettings"]["CutOff"] = self.ctx.general_calc_params["raspa"]["cutoff"]
        self.ctx.raspa_parameters["GeneralSettings"]["NumberOfInitializationCycles"] = self.inputs.raspa_parameters.get_dict()["GeneralSettings"]["NumberOfInitializationCycles"]
        self.ctx.raspa_parameters["GeneralSettings"]["NumberOfCycles"] = self.inputs.raspa_parameters.get_dict()["GeneralSettings"]["NumberOfCycles"]
        self.ctx.raspa_parameters["GeneralSettings"]["PrintPropertiesEvery"] = int(self.ctx.raspa_parameters["GeneralSettings"]["NumberOfCycles"] / self.ctx.general_calc_params["raspa"]["verbosity"])
        self.ctx.raspa_parameters["GeneralSettings"]["PrintEvery"] = int(1e6) #never
        return

    def run_raspa_gemc(self):
        """
        It runs a Widom calculation in RASPA.
        """
        temperature = self.ctx.temperature[0][self.ctx.current_T_index]
        for key, value in self.ctx.raspa_box.items():
            box_tag = value.tag
            self.ctx.raspa_parameters["System"][box_tag]["ExternalTemperature"] = temperature
        # Create the inputs dictionary
        inputs = {
            "code"       : self.inputs.raspa_code,
            "metadata"  :{
                "options":  self.ctx.raspa_options,
                "label"  : "RASPA GEMC Calculation",
                "description" : "RASPA GEMC Calculation",
            }
        }
        # Handling the retrive and usage of restart feature.
        # if self.ctx.restart_raspa_calc is not None:
           # inputs['retrieved_parent_folder'] = self.ctx.restart_raspa_calc

        # Storing parameters, and submittig the process.
        inputs["parameters"] = ParameterData(dict=self.ctx.raspa_parameters).store()
        gemc = self.submit(RaspaCalculation, **inputs)
        self.report("pk: <{}> | Running Raspa GEMC calculation".format(gemc.pk))
        self.ctx.nruns += 1
        return ToContext(raspa_gemc=gemc)

    def is_final_box_ok(self):
        """
        Here, we inspect if the box BoxLengths still satisfy the
        minimum image convetion or not!
        If not, increase the BoxLengths and restart.

        """
        # 1st run.
        if self.ctx.nruns == 0:
            min_image_stat = [False]
        # Worksround for the 2nd and more calculation.
        else:
            min_image_stat = []
            output_gemc = self.ctx.raspa_gemc.outputs.output_parameters.get_dict()
            for key, value in self.ctx.raspa_box.items():
                tag = value.tag
                box_ax = output_gemc[tag]['general']['box_ax']
                box_by = output_gemc[tag]['general']['box_by']
                box_cz = output_gemc[tag]['general']['box_cz']
                if box_ax > 2 * self.ctx.general_calc_params["raspa"]["cutoff"]:
                    min_image_stat.append(True)
                else:
                    min_image_stat.append(False)
                if box_by > 2 * self.ctx.general_calc_params["raspa"]["cutoff"]:
                    min_image_stat.append(True)
                else:
                    min_image_stat.append(False)
                if box_cz > 2 * self.ctx.general_calc_params["raspa"]["cutoff"]:
                    min_image_stat.append(True)
                else:
                    min_image_stat.append(False)

            if all(min_image_stat):
                self.report("The final avergae box dimension still satisfies the minimum image convetion")
            else:
                # Update the box
                for key, value in self.ctx.raspa_box.items():
                    box_tag = value.tag
                    box_ax = value.box_ax + 2.0
                    box_by = value.box_by + 2.0
                    box_cz = value.box_cz + 2.0
                    self.ctx.raspa_parameters["System"][box_tag]["BoxLengths"] = "{} {} {}".format(box_ax, box_by, box_cz)
                ParameterData(dict=self.ctx.raspa_parameters).store()
                self.report("Ooops! You need to increase the box lengths.")
                self.report("Increasing the cell dimension")

        return not all(min_image_stat)

    def should_run_another_gemc(self):
        """
        We run another RASPA GEMC calculation only if the current temperature
        is smaller than the maximum temperature.
        """
        return not self.ctx.done

    def gemc_conv_check(self):
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
            output_gemc = self.ctx.raspa_gemc.outputs.output_parameters.get_dict()
            # self.ctx.restart_raspa_calc = self.ctx.raspa_gcmc.outputs['retrieved']
            for key, value in self.ctx.raspa_box.items():
                for k, v in self.ctx.raspa_comp.items():
                    comp_name = v.name
                    conv_threshold = v.conv_threshold
                    tag = value.tag

                    adsorbate_density_average = output_gemc[tag]["components"][comp_name]["adsorbate_density_average"]
                    adsorbate_density_dev = output_gemc[tag]["components"][comp_name]["adsorbate_density_dev"]
                    error = round((adsorbate_density_dev / adsorbate_density_average), 2)
                    if error <= conv_threshold:
                        self.report("Error(%) <{:.1f}> less than <{:.1f}> | GEMC simulation is converged for <{}>".format(error*100, conv_threshold*100, comp_name))
                        converged.append(True)
                    else:
                        self.report("Error(%) <{:.1f}> greater than <{:.1f}> | GEMC simulation IS NOT converged for <{}>".format(error*100, conv_threshold*100, comp_name))
                        converged.append(False)

        # Get pressure for having in the self.report
        temperature = self.ctx.temperature[0][self.ctx.current_T_index]
        additional_cycle = self.ctx.general_calc_params['raspa']['additional_cycle']
        # First run.
        if self.ctx.nruns == 0:
            self.ctx.raspa_parameters["GeneralSettings"]["NumberOfInitializationCycles"] = self.ctx.raspa_parameters["GeneralSettings"]["NumberOfInitializationCycles"]
            self.ctx.raspa_parameters["GeneralSettings"]["NumberOfCycles"] = self.ctx.raspa_parameters["GeneralSettings"]["NumberOfCycles"]
            self.report("<{}> GEMC simulation at T(K)=<{:.1f}>".format(ordinal(self.ctx.nruns+1),temperature))
        # Second and more runs which are not converged.
        elif (not all(converged)) and self.ctx.nruns != 0:
            self.ctx.raspa_parameters["GeneralSettings"]["NumberOfInitializationCycles"] = self.ctx.raspa_parameters["GeneralSettings"]["NumberOfInitializationCycles"] + additional_cycle
            self.ctx.raspa_parameters["GeneralSettings"]["NumberOfCycles"] = self.ctx.raspa_parameters["GeneralSettings"]["NumberOfCycles"] + additional_cycle
            ParameterData(dict=self.ctx.raspa_parameters).store()
            new_init = int(self.ctx.raspa_parameters["GeneralSettings"]["NumberOfInitializationCycles"])
            new_prod = int(self.ctx.raspa_parameters["GeneralSettings"]["NumberOfCycles"])
            self.report("<{}> GEMC simulation at T(K)=<{:.1f}>".format(ordinal(self.ctx.nruns+1),temperature))
            self.report("Increasing MC cycles to Initialization<{:d}> -- Production<{:d}>".format(new_init,new_prod))
        # Converged run.
        else:
            self.report("GEMC simulation at T(K)=<{:.1f}> is converged after <{}> run".format(temperature,ordinal(self.ctx.nruns)))
            self.report("Resetting MC cycles and run counter")
            self.ctx.nruns = 0
            init = int(self.inputs.raspa_parameters.get_dict()["GeneralSettings"]["NumberOfInitializationCycles"])
            prod = int(self.inputs.raspa_parameters.get_dict()["GeneralSettings"]["NumberOfCycles"])
            self.ctx.raspa_parameters["GeneralSettings"]["NumberOfInitializationCycles"] = init
            self.ctx.raspa_parameters["GeneralSettings"]["NumberOfCycles"] = prod
            ParameterData(dict=self.ctx.raspa_parameters).store()

        return not all(converged)

    def parse_raspa_gemc(self):
        """
        Extract the pressure and loading average of the last completed and converged
        RASPA calculation.
        """
        temperature = self.ctx.temperature[0][self.ctx.current_T_index]
        output_gemc = self.ctx.raspa_gemc.outputs.output_parameters.get_dict()
        self.report("{} Current -- {} Initial".format(self.ctx.current_T_index,self.ctx.init_T_index))

        # Creating the loading empty dictionary only at first run.
        if self.ctx.current_T_index == self.ctx.init_T_index:
            self.report("creating the dictionary")
            self.ctx.density = {}
            for key, value in self.ctx.raspa_box.items():
                for k,v in self.ctx.raspa_comp.items():
                    comp_name = v.name
                    tag = value.tag
                    self.ctx.density[comp_name]= {}
                    self.ctx.density[comp_name][tag] = []
                    self.report("creating nested dictionary")

        # Iterating over components and append the loadings and
        # error bars to the dictionary.
        for key, value in self.ctx.raspa_box.items():
            for k,v in self.ctx.raspa_comp.items():
                comp_name = v.name
                tag = value.tag
                adsorbate_density_average = output_gemc[tag]["components"][comp_name]["adsorbate_density_average"]
                adsorbate_density_dev = output_gemc[tag]["components"][comp_name]["adsorbate_density_dev"]
                self.ctx.density[comp_name][tag].append([temperature,adsorbate_density_average,adsorbate_density_dev])

        # Simulation is converged that we are parsing. Increase the pressure index.
        if self.ctx.uphill:
            self.ctx.current_T_index += 1
            if self.ctx.current_T_index > self.ctx.max_T_index:
                self.ctx.uphill = False
                self.ctx.downhill = True
                self.ctx.current_T_index = self.init_T_index - 1
        elif self.ctx.downhill:
            self.ctx.current_T_index -= 1
            if self.ctx.current_T_index < self.ctx.min_T_index:
                self.ctx.done = True
        return

    def return_results(self):
        """
        Attach the results to the output.
        """
        # Create empty results dictionary.
        result_dict = {}

        # RASPA Section
        try:
            # Getting the output parameters of converged GCMC calculation.
            output_gemc = self.ctx.raspa_gemc.outputs.output_parameters.get_dict()
            result_dict["density"] = {}
            result_dict["header"] = ["Temperature(K)", "Density_average(kg/m^3)", "Density_deviation(kg/m^3)"]
            # Creating the needed keys within the result dictionary for RASPA Widom results.
            for key, value in self.ctx.raspa_box.items():
                for k,v in self.ctx.raspa_comp.items():
                    comp_name = v.name
                    tag = value.tag
                    result_dict["density"][comp_name] = {}
                    result_dict["density"][comp_name][tag] = {}

            # Iterating over components and extract results.
            for key, value in self.ctx.raspa_box.items():
                for k,v in self.ctx.raspa_comp.items():
                    comp_name = v.name
                    tag = value.tag
                    result_dict["adsorbate_density"][comp_name][tag] = self.ctx.density[comp_name][tag]

        except AttributeError:
            pass


        # Finalizing the results and report!
        self.out("results", ParameterData(dict=result_dict).store())
        self.report("Workchain completed successfully! | Result Dict is <{}>".format(self.outputs["results"].pk))
        return
# EOF
