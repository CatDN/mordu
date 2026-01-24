"""Algorithms for the calculation of common thermodynamic values

Includes:

- saturation pressure of a pure fluid at a specific temperature
"""

from .eos import EOS
from .symbols import *
from .utilities import multi_root_x, cubic_root_density

import numpy as np

def calc_saturation_pressure(eos: EOS, temperatures:list, max_iter: int = 100, tol: float = 1e-6, debug: bool = False) -> list:
    """For the calculation of the saturation pressure of a specific fluid
    from a specific EOS at specific temperatures.

    It uses an iterative method to calculate the saturation pressure
    at specific temperatures at a given tolerance.

    Parameters
    ----------
    eos: EOS
        The specific EOS to be used to calculate the saturation pressure
    temperatures: list
        List of temperatures in [K] at which to calculate saturation pressure.
        Values between the triple and critical points only.
    max_iter: int
        Maximum number of iterations for the iterative algorithm.
    tol: float = 1e-6
        Tolerance for iterative method
    debug: bool = False
        Debug flag. If 'True' then useful debugging print statements will appear.

    Returns
    -------
    saturation_pressure: list
        List of saturation pressure in [Pa] at each corresponding temperature.
    """
    
    # guesses and stores
    saturation_pressure = []
    pressure_guess = eos.fluid.P_t
   
    # equations and expressions from EOS
    P_equation = sp.lambdify((rho, P, T), P-eos.pressure)
    fugacity_coefficient = sp.lambdify((rho, T), eos.fugacity_coefficient)
    
    # root finding method
    if eos.name in ["vdW", "PR", "RK", "MSRK"]:
        root_finding = cubic_root_density
        args = [eos]
    else:
        root_finding = multi_root_x
        densities = np.logspace(-1, 5, int(1e4))
        args = [P_equation, densities]

    density_roots = root_finding(*args, args=(pressure_guess, temperatures[0]))

    for temperature in temperatures:

        # iteration counter
        counter = 0

        # multiplication factor is inversely proportional to distance to critical temperature
        factor = 1 + max(5e-3, (1 - (temperature-eos.fluid.T_t)/(eos.fluid.T_c - eos.fluid.T_t)))

        if debug:
            print(f"factor = {factor}")

        while len(density_roots)<3 and pressure_guess < 5e7:
            # calculate the density roots from pressure and temperature
            pressure_guess = pressure_guess*factor
            density_roots = root_finding(*args, args=(pressure_guess, temperature))

        if debug:
            print(f"guess = {pressure_guess}")

        # calculate the fugacity coefficients
        phi_v = fugacity_coefficient(density_roots[0], temperature)
        phi_l = fugacity_coefficient(density_roots[-1], temperature)

        # check their ratio
        while abs(phi_v/phi_l - 1) > tol and pressure_guess < 5e7 and counter < max_iter:
            if phi_v > phi_l:   # the vapour wants to escape, so decrease pressure
                pressure_guess = pressure_guess - pressure_guess*min(0.2, abs(phi_v/phi_l-1)**(1)) # min(1e2, pressure_guess*0.01)
            if phi_l > phi_v:   # the liquid want to escape so increase the pressure
                pressure_guess = pressure_guess + pressure_guess*min(0.2, abs(phi_v/phi_l-1)**(1)) # min(1e2, pressure_guess*0.01)

            # recalculate density roots
            density_roots = root_finding(*args, args=(pressure_guess, temperature))

            # make the factor dependent on the error metric
            factor = 1 + abs(phi_v/phi_l - 1)**1.5  # power of 1.5 slows down the factor the smaller the error metric
            # recheck for number of roots
            while len(density_roots)<3 and pressure_guess < 5e7:
                # calculate the density roots from pressure and temperature
                pressure_guess = pressure_guess*factor
                density_roots = root_finding(*args, args=(pressure_guess, temperature))

            # calculate the fugacity coefficients
            phi_v = fugacity_coefficient(density_roots[0],temperature)
            phi_l = fugacity_coefficient(density_roots[-1],temperature)

            counter += 1

            if debug:
                print(f"phi_l = {phi_l}, phi_v = {phi_v}")
                print(f"roots = {density_roots}, pressure = {pressure_guess}, phi_v/phi_l = {phi_v/phi_l}")

        if pressure_guess>5e7 or counter>= max_iter:
            print(f"Pressure point couldn't be found, temperature = {temperature} K")
            saturation_pressure += [0]
            pressure_guess = max(max(saturation_pressure), eos.fluid.P_t)
            continue

        saturation_pressure += [pressure_guess]

        if debug:
            print(f"pressure({temperature} K) = {pressure_guess/1e5} bar")

    return saturation_pressure