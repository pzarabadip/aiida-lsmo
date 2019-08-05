def choose_pressure_points(**kwargs):
    """Model the isotherm as single-site langmuir and return the most important
    pressure points to evaluate for an isotherm. Returns a list of pressures.

    :param Kh: Henry coefficient (mol/kg/Pa)
    :param qsat: saturations loading (mol/kg)
    :param dpa: precision of the sampling at low pressure (0.1 is a good one)
    :param dpmax: maximum distance between two pressure points (Pa)
    :param pmax: max pressure to sample (Pa)
    """
    dynamic = kwargs.pop('dynamic')
    full = kwargs.pop('full')


    if dynamic and full:
        Kh = kwargs.pop('Kh')
        qsat = kwargs.pop('qsat')
        dpmax = kwargs.pop('dpmax')
        pmin = kwargs.pop('pmin')
        dpa = kwargs.pop('dpa')
        pmax = kwargs.pop('pmax')
        R = 8.314/1000 #(kJ/mol/K)
        b = Kh/qsat #(1/Pa)
        # pmin = 0.001e5 #(Pa)
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
    elif (not dynamic) and full:
        pmin = kwargs.pop('pmin')
        dpa = kwargs.pop('dpa')
        pmax = kwargs.pop('pmax')
        p = [pmin]
        while True:
            pold = p[-1]
            pnew = pold + dpa
            if pnew < pmax:
                p.append(pnew)
            else:
                p.append(pmax)
                return p
    else:
        p = kwargs.pop('selected_pressures')
        return p
