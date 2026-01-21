"""Algorithms for the calculation of common thermodynamic values

"""

from .eos import EOS
from .symbols import *
from .utilities import multi_root_x, cubic_root_density
import numpy as np

def calc_saturation_pressure(eos: EOS, temperatures:list, tol: float = 1e-6, debug: bool = False) -> list:
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

    for temperature in temperatures:
        density_roots = []
        while len(density_roots)<3 and pressure_guess < 5e7:
            # calculate the density roots from pressure and temperature
            density_roots = root_finding(*args, pressure_guess, temperature)
            pressure_guess = pressure_guess*1.0005

        if debug:
            print(f"guess = {pressure_guess}")

        # calculate the fugacity coefficients
        phi_v = fugacity_coefficient(density_roots[0], temperature)
        phi_l = fugacity_coefficient(density_roots[-1], temperature)

        # check their ratio
        while abs(phi_v/phi_l - 1) > tol and pressure_guess < 5e7:
            if phi_v > phi_l:   # the vapour wants to escape, so increase pressure
                pressure_guess = pressure_guess - pressure_guess*min(0.8, abs(phi_v/phi_l-1)) # min(1e2, pressure_guess*0.01)
            if phi_l > phi_v:   # the liquid want to escape so decrease the pressure
                pressure_guess = pressure_guess + pressure_guess*min(0.8, abs(phi_v/phi_l-1)) # min(1e2, pressure_guess*0.01)

            # recalculate density roots
            density_roots = cubic_root_density(eos, pressure_guess, temperature)

            # calculate the fugacity coefficients
            phi_v = fugacity_coefficient(density_roots[0],temperature)
            phi_l = fugacity_coefficient(density_roots[-1],temperature)

            if debug:
                print(f"phi_l = {phi_l}, phi_v = {phi_v}")
                print(f"roots = {density_roots}, pressure = {pressure_guess}, phi_v/phi_l = {phi_v/phi_l}")

        if pressure_guess>5e7:
            print("Pressure point couldn't be found")
            saturation_pressure += [0]
            pressure_guess = max(max(saturation_pressure), eos.fluid.T_t)
            continue

        saturation_pressure += [pressure_guess]

        if debug:
            print(f"pressure({temperature} K) = {pressure_guess/1e5} bar")

    return saturation_pressure