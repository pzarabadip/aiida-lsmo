"""
Smart High_Throughput Screening WorkChain
"""
from __future__ import absolute_import
import os
import six

# AiiDA modules
from aiida.plugins import CalculationFactory, DataFactory, WorkflowFactory
from aiida.orm import Dict, Float, Int, List, Str, SinglefileData
from aiida.engine import calcfunction
from aiida.engine import ToContext, WorkChain, if_, while_
from aiida_lsmo.utils import aiida_dict_merge, check_resize_unit_cell

RaspaBaseWorkChain = WorkflowFactory('raspa.base')

# Defining DataFactory and CalculationFactory
CifData = DataFactory("cif")
ZeoppParameters = DataFactory("zeopp.parameters")

ZeoppCalculation = CalculationFactory("zeopp.network")


@calcfunction
def get_components_dict(mixture_name):
    """Get a Dict from the provided pre-defined mixture"""
    import ruamel.yaml as yaml
    thisdir = os.path.dirname(os.path.abspath(__file__))
    yamlfile = os.path.join(thisdir, "isotherm_data", mixture_name.value + ".yaml")
    with open(yamlfile, 'r') as stream:
        mixture_dict = yaml.safe_load(stream)
    return Dict(dict=mixture_dict)

@calcfunction
def get_atomic_radii(htsparam):
    """Get {forcefield}.rad as SinglefileData form workchain/isotherm_data"""
    thisdir = os.path.dirname(os.path.abspath(__file__))
    fullfilename = htsparam['forcefield'] + ".rad"
    return SinglefileData(file=os.path.join(thisdir, "isotherm_data", fullfilename))

@calcfunction
def get_zeopp_parameters(components, htsparams, probrad):
    """Get the ZeoppParameters from the inputs of the workchain"""
    # probe_rad = components[comp.value]['probrad']
    probe_rad = probrad.value
    param_dict = {
        'ha': 'DEF',
        'res': True,
        'sa' : [probe_rad, probe_rad, htsparams['zeopp_sa_samples']],
        'volpo': [probe_rad, probe_rad, htsparams['zeopp_volpo_samples']],
        'block': [probe_rad, htsparams['zeopp_block_samples']],
    }
    return ZeoppParameters(dict=param_dict)

@calcfunction
def get_geometric_output(zeopp_out):
    """Return the geometric_output Dict from Zeopp results, including Qsat and is_porous"""
    geometric_output = zeopp_out.get_dict()
    geometric_output.update({
        'is_porous': geometric_output["POAV_A^3"] > 0.000
    })
    return Dict(dict=geometric_output)

@calcfunction
def choose_pressure_points(inp_param):
    """If 'presure_list' is not provide, model the isotherm as single-site langmuir and return the most important
    pressure points to evaluate for an isotherm, in a List.
    """
    if inp_param["pressure_list"]:
        pressure_points = inp_param["pressure_list"]
    else:
        raise Exception('Pressue list is not provided!')
    return List(list=pressure_points)


# Deafault parameters
HTSPARAMETERS_DEFAULT = Dict(
    dict={  #TODO: create IsothermParameters instead of Dict # pylint: disable=fixme
        "forcefield": "UFF",  # str, Forcefield of the structure
        "ff_tailcorr": True,  # bool, Apply tail corrections
        "ff_shift": False,  # bool, Shift or truncate at cutoff
        "ff_cutoff": 12.0,  # float, CutOff truncation for the VdW interactions (Angstrom)
        "temperature": 300,  # float, Temperature of the simulation
        "temperature_list": None,  # list, to be used by IsothermMultiTempWorkChain
        "zeopp_volpo_samples": int(1e5),  # int, Number of samples for VOLPO calculation (per UC volume)
        "zeopp_block_samples": int(100),  # int, Number of samples for BLOCK calculation (per A^3)
        "raspa_minKh": 1e-10,  # float, If Henry coefiicient < raspa_minKh do not run the isotherm (mol/kg/Pa)
        "raspa_verbosity": 10,  # int, Print stats every: number of cycles / raspa_verbosity
        "raspa_widom_cycles": int(1e5),  # int, Number of widom cycles
        "raspa_gcmc_init_cycles": int(1e3),  # int, Number of GCMC initialization cycles
        "raspa_gcmc_prod_cycles": int(1e4),  # int, Number of GCMC production cycles
        "lcd_max" : 15.0, # mandatory
        "pld_scale" : 1.0, # mandatory
        "pressure_list": [0.1e5,1.0e5],  # list, Pressure list for the isotherm (bar): if given it will skip  guess
        "ideal_selectivity_threshold": 10.0, #mandatory if protocol is relative.
        "ideal_selectivity_protocol": 'ignore', # ignore, loose, and tight!
    })

class SmartHTSWorkChain(WorkChain):
    """
    The ZeoppRaspaWorkChain is designed to perform zeo++ and
    RASPA calculations.
    """

    @classmethod
    def define(cls, spec):
        super(SmartHTSWorkChain, cls).define(spec)
        """
        Define workflow specification.
        This is the most important method of a Workchain, which defines the
        inputs it takes, the logic of the execution and the outputs
        that are generated in the process.
        """
        # SmartHTSWorkChain inputs!
        spec.input("structure", valid_type=CifData, required=True, help="Input structure in cif format")
        spec.input("parameters", valid_type=Dict, default=HTSPARAMETERS_DEFAULT, required=True, help='It provides the parameters which control the decision making behavior of workchain.')
        spec.input("mixture", valid_type=Str, required=True, help='It correspnds to the name of yaml file which provides componenets information!')


        # Exposing the remaining inputs!
        spec.expose_inputs(ZeoppCalculation, namespace='zeopp', include=['code', 'metadata'])
        spec.expose_inputs(RaspaBaseWorkChain, namespace='raspa_base', exclude=['raspa.structure', 'raspa.parameters'])

        # Redefine after finishing.
        spec.outline(
            cls.setup,
            cls.run_zeopp,
            cls.inspect_zeopp_calc,
            if_(cls.should_run_widom)(
                cls.run_raspa_widom,
                cls.inspect_widom_calc,
                if_(cls.should_run_gcmc)(
                    cls.init_raspa_gcmc,
                    while_(cls.should_run_another_gcmc)(
                        cls.run_raspa_gcmc,
                        cls.parse_raspa_gcmc,
                    ),
                )
            ),
            cls.return_results,
        )

        # to be returned
        spec.outputs.dynamic = True

    def setup(self):
        """
        Initialize variables and setup screening protocol!
        """
        # Getting the components dict.
        if isinstance(self.inputs.mixture, Str):
            self.ctx.components = get_components_dict(self.inputs.mixture)
        else:
            raise Exception('Mixture is not provided properly!')

        # Getting and updating the HTS parameters!
        self.ctx.parameters = aiida_dict_merge(HTSPARAMETERS_DEFAULT, self.inputs.parameters)

        # Get integer temperature in context for easy reports
        self.ctx.temperature = int(round(self.ctx.parameters['temperature']))

    def run_zeopp(self):
        """
        It performs the zeopp pore diameter calculation.
        """
        # Required inputs
        zeopp_inputs = self.exposed_inputs(ZeoppCalculation, 'zeopp')
        zeopp_inputs['structure'] = self.inputs.structure
        zeopp_inputs['metadata']['label'] = 'ZeoppVolpoBlock'
        zeopp_inputs['metadata']['call_link_label'] = 'run_zeopp'
        zeopp_inputs['metadata']['description'] = 'Called by SmartHTSWorkChain'
        zeopp_inputs['atomic_radii'] = get_atomic_radii(self.ctx.parameters)

        for key, value in self.ctx.components.get_dict().items():
            comp = value['name']
            proberad = Float(value['proberad'])
            # zeopp_inputs['parameters'] = get_zeopp_parameters(self.ctx.components, self.ctx.parameters, Str(comp))
            zeopp_inputs['parameters'] = get_zeopp_parameters(self.ctx.components, self.ctx.parameters, proberad)
            running = self.submit(ZeoppCalculation, **zeopp_inputs)
            zeopp_label = "zeopp_{}".format(comp)
            self.report("Running zeo++ block and volpo Calculation<{}>".format(running.id))
            self.to_context(**{zeopp_label:running})

    def inspect_zeopp_calc(self):
        """Checks if all zeopp calculations are finished ok."""
        for key, value in self.ctx.components.get_dict().items():
            zeopp_label = "zeopp_{}".format(value['name'])
            assert self.ctx[zeopp_label].is_finished_ok

    def should_run_widom(self):
        """Decided whether to run widom or not"""
        self.ctx.should_run_widom = []
        self.ctx.geom = {}
        lcd_lim = self.ctx.parameters["lcd_max"]
        for key, value in self.ctx.components.get_dict().items():
            comp = value['name']
            zeopp_label = "zeopp_{}".format(comp)
            pld_lim = value["proberad"] * self.ctx.parameters["pld_scale"]
            self.ctx.geom[comp] = get_geometric_output(self.ctx[zeopp_label].outputs.output_parameters)
            pld_component = self.ctx.geom[comp]["Largest_free_sphere"]
            lcd_component = self.ctx.geom[comp]["Largest_included_sphere"]
            if (lcd_component <= lcd_lim) and (pld_component >= pld_lim) and (self.ctx.geom[comp]['is_porous']):
                self.report("ALL pre-selection conditions are satisfied: Calculate Henry coefficients")
                self.report("Found {} blocking spheres".format(self.ctx.geom[comp]['Number_of_blocking_spheres']))
                if self.ctx.geom[comp]['Number_of_blocking_spheres'] > 0:
                    self.out_many(self.exposed_outputs(self.ctx[zeopp_label], ZeoppCalculation))
                self.ctx.should_run_widom.append(True)
            else:
                self.report("All/Some of pre-selection criteria are NOT met: terminate!")
                self.ctx.should_run_widom.append(False)

        self.out("geometric_output", self.ctx.geom)
        return all(self.ctx.should_run_widom)

    def _get_widom_param(self):
        """Write Raspa input parameters from scratch, for a Widom calculation"""

        param = {
            "GeneralSettings": {
                "SimulationType": "MonteCarlo",
                "NumberOfInitializationCycles": 0,
                "NumberOfCycles": self.ctx.parameters['raspa_widom_cycles'],
                "PrintPropertiesEvery": self.ctx.parameters['raspa_widom_cycles'] / self.ctx.parameters['raspa_verbosity'],
                "PrintEvery": int(1e10),
                "RemoveAtomNumberCodeFromLabel": True,  # be careful!
                "Forcefield": self.ctx.parameters['forcefield'],
                "CutOff": self.ctx.parameters['ff_cutoff'],
            },
            "System": {
                "framework_1": {
                    "type": "Framework",
                    "ExternalTemperature": self.ctx.parameters['temperature'],
                }
            },
            "Component": {},
        }

        # Check particular conditions and settings
        mult = check_resize_unit_cell(self.inputs.structure, 2 * self.ctx.parameters['ff_cutoff'])
        param["System"]["framework_1"]["UnitCells"] = "{} {} {}".format(mult[0], mult[1], mult[2])

        return param

    def run_raspa_widom(self):
        """
        It generates the ParameterData for RASPA and use the blocking spheres if there is any.
        With current plugin we can have more than one component to be used for Widom insertion.
        """
        self.ctx.raspa_inputs = self.exposed_inputs(RaspaBaseWorkChain, 'raspa_base')
        self.ctx.raspa_inputs['metadata']['label'] = "RaspaWidom"
        self.ctx.raspa_inputs['metadata']['description'] = "Called by SmartHTSWorkChain"
        self.ctx.raspa_inputs['metadata']['call_link_label'] = "run_raspa_widom"

        self.ctx.raspa_inputs['raspa']['framework'] = {"framework_1": self.inputs.structure}

        for key, value in self.ctx.components.get_dict().items():
            comp = value['name']
            zeopp_label = "zeopp_{}".format(comp)
            self.ctx.raspa_param = self._get_widom_param()
            self.ctx.raspa_inputs["raspa"]["block_pocket"] = {}
            self.ctx.raspa_param["System"]["framework_1"]["HeliumVoidFraction"] = self.ctx.geom[comp]["POAV_Volume_fraction"]
            self.ctx.raspa_param["Component"][comp] = {}
            self.ctx.raspa_param["Component"][comp]["MoleculeDefinition"] = value['forcefield']
            self.ctx.raspa_param["Component"][comp]["WidomProbability"] = 1.0
            if self.ctx.geom[comp]['Number_of_blocking_spheres'] > 0:
                self.ctx.raspa_inputs["raspa"]["block_pocket"] = {comp + "_block_file": self.ctx[zeopp_label].outputs.block}
                self.ctx.raspa_param["Component"][comp]["BlockPocketsFileName"] = comp + "_block_file"
            if value['charged']:
                self.ctx.raspa_param["GeneralSettings"].update({"UseChargesFromCIFFile":"yes", "ChargeMethod": "Ewald", "EwaldPrecision": 1e-6})

            self.ctx.raspa_inputs['raspa']['parameters'] = Dict(dict=self.ctx.raspa_param)
            running = self.submit(RaspaBaseWorkChain, **self.ctx.raspa_inputs)
            widom_label = "widom_{}".format(comp)
            self.report("Running Raspa Widom @ {}K for the Henry coefficient".format(self.ctx.temperature))
            self.to_context(**{widom_label:running})

    def inspect_widom_calc(self):
        """
        Checks if all zeopp calculations are finished ok.
        """
        for key, value in self.ctx.components.get_dict().items():
            widom_label = "widom_{}".format(value['name'])
            assert self.ctx[widom_label].is_finished_ok

    def should_run_gcmc(self):
        """
        Should be updated!
        """
        self.ctx.should_run_gcmc = []
        if self.ctx.parameters['ideal_selectivity_protocol'] == 'ignore':
            self.ctx.should_run_gcmc.append(True)

        if self.ctx.parameters['ideal_selectivity_protocol'] == 'loose':
            widom_label_comp1 = "widom_{}".format(self.ctx.components.get_dict()['comp1']['name'])
            widom_label_comp2 = "widom_{}".format(self.ctx.components.get_dict()['comp2']['name'])
            output1 = self.ctx[widom_label_comp1].outputs.output_parameters.get_dict()
            output2 = self.ctx[widom_label_comp2].outputs.output_parameters.get_dict()
            self.ctx.kh_comp1 = output1["framework_1"]["components"][self.ctx.components.get_dict()['comp1']['name']]["henry_coefficient_average"]
            self.ctx.kh_comp2 = output2["framework_1"]["components"][self.ctx.components.get_dict()['comp2']['name']]["henry_coefficient_average"]
            self.ctx.ideal_selectivity = self.ctx.kh_comp1 / self.ctx.kh_comp2
            self.ctx.ideal_selectivity_threshold = self.ctx.parameters["ideal_selectivity_threshold"]
            if self.ctx.ideal_selectivity >= self.ctx.ideal_selectivity_threshold:
                self.report("Ideal selectivity is greater than threshold: compute the GCMC")
                self.ctx.should_run_gcmc.append(True)
            else:
                self.report("Ideal selectivity is less than threshold: DO NOT compute the GCMC")
                self.ctx.should_run_gcmc.append(False)

        if self.ctx.parameters['ideal_selectivity_protocol'] == 'tight':
            widom_label_comp1 = "widom_{}".format(self.ctx.components.get_dict()['comp1']['name'])
            output1 = self.ctx[widom_label_comp1].outputs.output_parameters.get_dict()
            self.ctx.kh_comp1 = output1["framework_1"]["components"][self.ctx.components.get_dict()['comp1']['name']]["henry_coefficient_average"]
            for key, value in self.ctx.components.get_dict().items():
                comp = value['name']
                widom_label = "widom_{}".format(comp)
                output = self.ctx[widom_label].outputs.output_parameters.get_dict()
                self.ctx.kh_comp = output["framework_1"]["components"][comp]["henry_coefficient_average"]
                self.ctx.ideal_selectivity = self.ctx.kh_comp1 / self.ctx.kh_comp
                if self.ctx.ideal_selectivity > self.ctx.ideal_selectivity_threshold:
                    self.ctx.should_run_gcmc.append(True)
                else:
                    self.ctx.should_run_gcmc.append(False)

        return all(self.ctx.should_run_gcmc)

    def init_raspa_gcmc(self):
        """Initialize RASPA gcmc"""
        # resetting the component section for having multi-comp simulation.
        self.ctx.raspa_param["Component"] = {}
        self.ctx.raspa_param["Component"] = {item:{} for index, item in enumerate(list(self.ctx.components.get_dict()))}

        # Initializate counter and set restart to None
        self.ctx.current_p_index = 0
        self.ctx.restart_raspa_calc = None

        self.ctx.pressures = choose_pressure_points(self.ctx.parameters)

        self.report("<{}> points are chosen for GCMC calculation".format(len(self.ctx.pressures)))

        # Iterating over components and update the input.
        for key, value in self.ctx.components.get_dict().items():
            comp = value['name']
            mol_frac = value['molfraction']
            singlebead = value['singlebead']
            self.ctx.raspa_param["Component"][comp] = self.ctx.raspa_param["Component"].pop(key)
            self.ctx.raspa_param["Component"][comp]["MolFraction"] = float(mol_frac)
            self.ctx.raspa_param["Component"][comp]["TranslationProbability"] = 0.5
            # Only adds RotationProbability move if it is not singlebead model.
            if not singlebead:
                self.ctx.raspa_param["Component"][comp]["RotationProbability"] = 0.5
            self.ctx.raspa_param["Component"][comp]["ReinsertionProbability"] = 0.5
            self.ctx.raspa_param["Component"][comp]["SwapProbability"] = 1.0
            self.ctx.raspa_param["Component"][comp]["IdentityChangeProbability"] = 1.0
            self.ctx.raspa_param["Component"][comp]["NumberOfIdentityChanges"] = len(list(self.ctx.components.get_dict()))
            self.ctx.raspa_param["Component"][comp]["IdentityChangesList"] = [i for i in range(len(list(self.ctx.components.get_dict())))]
        return

    def should_run_another_gcmc(self):
        """
        We run another raspa calculation only if the current iteration is
        smaller than the total number of pressures we want to compute.
        """
        return self.ctx.current_p_index < len(self.ctx.pressures)

    def run_raspa_gcmc(self):
        """
        It submits Raspa calculation to RaspaBaseWorkchain.
        """
        self.ctx.raspa_inputs['metadata']['label'] = "RaspaGCMC_{}".format(self.ctx.current_p_index + 1)
        self.ctx.raspa_inputs['metadata']['description'] = 'Called by SmartHTSWorkChain'
        self.ctx.raspa_inputs['metadata']['call_link_label'] = "run_raspa_gcmc_{}".format(self.ctx.current_p_index + 1)

        self.ctx.raspa_param["System"]["framework_1"]["ExternalPressure"] = self.ctx.pressures[self.ctx.current_p_index]

        # Create the input dictionary
        self.ctx.raspa_inputs["raspa"]["block_pocket"] ={}

        for key, value in self.ctx.components.get_dict().items():
            comp = value['name']
            zeopp_label = "zeopp_{}".format(comp)
            bp_label = comp + "_block_file"

            if self.ctx.geom[comp]['Number_of_blocking_spheres'] > 0:
                self.ctx.raspa_param["Component"][comp]["BlockPocketsFileName"] = {}
                self.ctx.raspa_param["Component"][comp]["BlockPocketsFileName"]["framework_1"] = bp_label
                self.ctx.raspa_inputs["raspa"]["block_pocket"][bp_label] = self.ctx[zeopp_label].outputs.block

        if self.ctx.restart_raspa_calc is not None:
            self.ctx.raspa_inputs["raspa"]['retrieved_parent_folder'] = self.ctx.restart_raspa_calc

        self.ctx.raspa_inputs['raspa']['parameters'] = Dict(dict=self.ctx.raspa_param)

        # Create the calculation process and launch it
        running = self.submit(RaspaBaseWorkChain, **self.ctx.raspa_inputs)
        self.report("pk: <{}> | Submitted RaspaBaseWorkchain at p(bar)={:.3f} <{} of {}>".format(running.pk, self.ctx.pressures[self.ctx.current_p_index]/1e5, self.ctx.current_p_index+1, len(self.ctx.pressures)))
        # self.ctx.nruns += 1
        return ToContext(raspa_gcmc=running)

    def parse_raspa_gcmc(self):
        """
        Extract the pressure and loading average of the last completed and converged
        RASPA calculation.
        """
        pressure = self.ctx.raspa_param["System"]["framework_1"]["ExternalPressure"]/1e5
        output_gcmc = self.ctx.raspa_gcmc.outputs.output_parameters.get_dict()
        self.ctx.restart_raspa_calc = self.ctx.raspa_gcmc.outputs['retrieved']
        # Creating the loading empty dictionary only at first run.
        if self.ctx.current_p_index == 0:
            self.ctx.loading = {}
            self.ctx.enthalpy = {}
            for key, value in self.ctx.components.get_dict().items():
                comp = value['name']
                self.ctx.loading[comp] = []
                self.ctx.enthalpy[comp] = []

        # Iterating over components and append the loadings and
        # error bars to the dictionary.
        conv1 =  1.0 / 120.273
        for key, value in self.ctx.components.get_dict().items():
            comp = value['name']
            conv2 = output_gcmc["framework_1"]["components"][comp]["conversion_factor_molec_uc_to_mol_kg"]
            loading_average_comp = conv2 * output_gcmc["framework_1"]["components"][comp]["loading_absolute_average"]
            loading_dev_comp = conv2 * output_gcmc["framework_1"]["components"][comp]["loading_absolute_dev"]
            enthalpy_average_comp = conv1 * output_gcmc["framework_1"]["components"][comp]["enthalpy_of_adsorption_average"]
            enthalpy_dev_comp = conv1 * output_gcmc["framework_1"]["components"][comp]["enthalpy_of_adsorption_dev"]

            self.ctx.loading[comp].append([pressure, loading_average_comp, loading_dev_comp])
            self.ctx.enthalpy[comp].append([pressure, enthalpy_average_comp, enthalpy_dev_comp])

        # Simulation is converged that we are parsing. Increase the pressure index.
        self.ctx.current_p_index += 1
        return

    def return_results(self):
        """
        Attach the results to the output.
        """
        # TODO: Benefit from calcfuntion!
        # Create empty results dictionary.
        result_dict = {}

        # Zeopp section
        try:
            # Output of Di and Df

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
            for key, value in self.ctx.components.get_dict().items():
                comp = value['name']
                zeopp_label = "zeopp_{}".format(comp)
                output_zeo = self.ctx[zeopp_label].outputs.output_parameters.get_dict()

                # result_dict["Largest_included_sphere"] = self.ctx[zeopp_label].outputs.output_parameters.get_dict()["Largest_included_sphere"]
                # result_dict["Largest_free_sphere"] = self.ctx[zeopp_label].outputs.output_parameters.get_dict()["Largest_free_sphere"]

                result_dict["POAV_Volume_fraction"][comp] = output_zeo["POAV_Volume_fraction"]
                result_dict["PONAV_Volume_fraction"][comp] = output_zeo["PONAV_Volume_fraction"]
                result_dict["POAV"][comp] = output_zeo["POAV_cm^3/g"]
                result_dict["GASA"][comp] = output_zeo["ASA_m^2/g"]
                # result_dict["VASA"][comp] = output_zeo["ASA_m^2/cm^3"]
                result_dict["NGASA"][comp] = output_zeo["NASA_m^2/g"]
                result_dict["NVASA"][comp] = output_zeo["NASA_m^2/cm^3"]
                result_dict["Channel_surface_area"][comp] = output_zeo["Channel_surface_area_A^2"]
                result_dict["Pocket_surface_area"][comp] = output_zeo["Pocket_surface_area_A^2"]
                # result_dict["number_blocking_spheres"][comp] = self.ctx.number_blocking_spheres[comp]
                result_dict["Density"] = output_zeo["Density"]
        except AttributeError:
            pass

        # RASPA widom Section
        try:
            # Getting the output parameters of widom calculation.
            # Creating the needed keys within the result dictionary for RASPA Widom results.
            result_dict["henry_coefficient_average"] = {}
            result_dict["henry_coefficient_dev"] = {}
            result_dict["adsorption_energy_average"] = {}
            result_dict["adsorption_energy_dev"] = {}
            result_dict["temperature"] = self.ctx.raspa_param["System"]["framework_1"]["ExternalTemperature"]
            result_dict["temperature_unit"] = "K"

            # Iterating over components and extract results.
            for key, value in self.ctx.components.get_dict().items():
                comp = value['name']
                widom_label = "widom_{}".format(comp)
                output_widom = self.ctx[widom_label].outputs.output_parameters.get_dict()

                compout = output_widom["framework_1"]["components"][comp]
                result_dict["henry_coefficient_unit"] = compout["henry_coefficient_unit"]
                result_dict["adsorption_energy_unit"] = compout["adsorption_energy_widom_unit"]
                result_dict["henry_coefficient_average"][comp] = compout["henry_coefficient_average"]
                result_dict["henry_coefficient_dev"][comp] = compout["henry_coefficient_dev"]
                result_dict["adsorption_energy_average"][comp] = compout["adsorption_energy_widom_average"]
                result_dict["adsorption_energy_dev"][comp] = compout["adsorption_energy_widom_dev"]
        except AttributeError:
            pass

        # RASPA Section
        try:
            # Getting the output parameters of converged GCMC calculation.
            output_gcmc = self.ctx.raspa_gcmc.outputs.output_parameters.get_dict()
            # Creating the needed keys within the result dictionary for RASPA Widom results.
            result_dict["conversion_factor_molec_uc_to_cm3stp_cm3"] = {}
            result_dict["conversion_factor_molec_uc_to_gr_gr"] = {}
            result_dict["conversion_factor_molec_uc_to_mol_kg"] = {}
            result_dict["isotherm_loading_header"] = ["Pressure(bar)", "Loading_average(mol/kg)", "Loading_deviation(mol/kg)"]
            result_dict["isotherm_loading"] = {}
            result_dict["mol_fraction"] = {}
            result_dict["isotherm_enthalpy"] = {}
            # Iterating over components and extract results.
            for key, value in self.ctx.components.get_dict().items():
                comp = value['name']
                #mol_frac = value['molfraction']
                result_dict["conversion_factor_molec_uc_to_cm3stp_cm3"][comp] = output_gcmc["framework_1"]["components"][comp]["conversion_factor_molec_uc_to_cm3stp_cm3"]
                result_dict["conversion_factor_molec_uc_to_gr_gr"][comp] = output_gcmc["framework_1"]["components"][comp]["conversion_factor_molec_uc_to_gr_gr"]
                result_dict["conversion_factor_molec_uc_to_mol_kg"][comp] = output_gcmc["framework_1"]["components"][comp]["conversion_factor_molec_uc_to_mol_kg"]
                result_dict["mol_fraction"][comp] = output_gcmc["framework_1"]["components"][comp]["mol_fraction"]
                result_dict["isotherm_loading"][comp] = self.ctx.loading[comp]
                result_dict["isotherm_enthalpy"][comp] = self.ctx.enthalpy[comp]
        except AttributeError:
            pass

        # Blocking Spheres Section
        # for key, value in self.ctx.raspa_comp.items():
        #     comp = value['name']
        #     zeopp_label = "zeopp_{}".format(comp)
        #     bp_label = comp + "_block_file"
        #     # Only keeping non-empty block pocket files.
        #     if self.ctx.number_blocking_spheres[comp] > 0:
        #         self.out("blocking_spheres_{}".format(comp), self.ctx[zeopp_label].outputs.block)

        # Finalizing the results and report!
        self.out("results", Dict(dict=result_dict))
        self.report("Workchain completed successfully! | Result Dict is <{}>".format(self.outputs["results"].pk))
        return
# EOF
