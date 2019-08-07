def choose_temp_points(T_min,T_max, dT):
    """Model the isotherm as single-site langmuir and return the most important
    pressure points to evaluate for an isotherm. Returns a list of pressures.

    :param Kh: Henry coefficient (mol/kg/Pa)
    :param qsat: saturations loading (mol/kg)
    :param dpa: precision of the sampling at low pressure (0.1 is a good one)
    :param dpmax: maximum distance between two pressure points (Pa)
    :param pmax: max pressure to sample (Pa)
    """
    import numpy as np

    T_mean = np.mean([T_min,T_max])
    T = set()
    T = [T_min]
    while True:
        Told = T[-1]
        Tnew = Told + dT
        if Tnew < T_max:
            T.append(Tnew)
        else:
            T.append(T_max)
            if T_mean not in T:
                T.append(T_mean)
            T.sort()
            return T, T_mean
