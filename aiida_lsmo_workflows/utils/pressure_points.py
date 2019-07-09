def choose_pressure_points(Kh, qsat, dpa, dpmax, pmax):
    """Model the isotherm as single-site langmuir and return the most important
    pressure points to evaluate for an isotherm. Returns a list of pressures.

    :param Kh: Henry coefficient (mol/kg/Pa)
    :param qsat: saturations loading (mol/kg)
    :param dpa: precision of the sampling at low pressure (0.1 is a good one)
    :param dpmax: maximum distance between two pressure points (Pa)
    :param pmax: max pressure to sample (Pa)
    """
    R = 8.314/1000 #(kJ/mol/K)
    b = Kh/qsat #(1/Pa)
    pmin = 0.001e5 #(Pa)
    p = [pmin]
    while True:
        pold = p[-1]
        dp = min(dpmax,dpa*(b*pold**2+2*pold+1/b))
        pnew = pold+dp
        if pnew <= pmax:
            p.append(pnew)
        else:
            p.append(pmax)
            return p

@workfunction
def update_raspa_parameters(parameters, pressure):
    """Store input parameters of Raspa for given pressure.

    Note: In order to keep the provenance of both Raspa calculations,
    changing the pressure force us to create a new ParameterData node.
    "workfunctions" will take care of linking the user-provided ParameterData
    node to the new one containing the pressure.
    """
    param_dict = parameters.get_dict()
    param_dict['GeneralSettings']['ExternalPressure'] = pressure.value
    return ParameterData(dict=param_dict)
