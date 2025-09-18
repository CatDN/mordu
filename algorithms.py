# 11/08/2025
# algorithms for the calculation of specific quantities

#libraries
from scipy import optimize
import numpy as np
from symbols import *
from other_functions import multi_root


# aid functions for algorithms

############################################################################### Cubic polynomials for density root finding of cubic EOS
def vdW_coefficients(a_function, b_function, pressure, temperature, z1, z2):
    #see log page 192
    a = a_function(temperature, z1, z2)
    b = b_function(z1, z2)

    coeff_0 = pressure
    coeff_1 = -pressure * b - R * temperature
    coeff_2 = a
    coeff_3 = - a * b

    return np.array([coeff_3, coeff_2, coeff_1, coeff_0])

def PR_coefficients(a_function, b_function, pressure, temperature, z1, z2):
    #see log page 192
    a = a_function(temperature, z1, z2)
    b = b_function(z1, z2)

    coeff_0 = pressure
    coeff_1 = pressure * b - R * temperature
    coeff_2 = -3*pressure*b**2 - 2*R*temperature*b + a
    coeff_3 = b**3*pressure + R*temperature*b**2 - a * b

    return np.array([coeff_3, coeff_2, coeff_1, coeff_0])

def RK_coefficients(a_function, b_function, pressure, temperature, z1, z2):
    #see log page 192
    a = a_function(temperature, z1, z2)
    b = b_function(z1, z2)

    coeff_0 = pressure
    coeff_1 = - R * temperature
    coeff_2 = -pressure*b**2 - R*temperature*b + a
    coeff_3 = - a * b

    return np.array([coeff_3, coeff_2, coeff_1, coeff_0])

def SRK_coefficients(a_function, b_function, pressure, temperature, z1, z2):
    #see log page 192
    a = a_function(temperature, z1, z2)
    b = b_function(z1, z2)

    coeff_0 = pressure
    coeff_1 = -R * temperature
    coeff_2 = -pressure*b**2 - b*R*temperature + a
    coeff_3 = - a * b

    return np.array([coeff_3, coeff_2, coeff_1, coeff_0])

def MSRK_coefficients(a_function, b_function, pressure, temperature, z1, z2):
    #see log page 192
    a = a_function(temperature, z1, z2)
    b = b_function(z1, z2)

    coeff_0 = pressure
    coeff_1 = -R * temperature
    coeff_2 = -pressure*b**2 - b*R*temperature + a
    coeff_3 = - a * b

    return np.array([coeff_3, coeff_2, coeff_1, coeff_0])

############################################################################### Wilson approximations
# Wilson approximation root finding problem for dew and bubble point pressure
def Wilson_P(variables, temperature, z1, z2, fluid_1, fluid_2, is_dew=True):
    """
    Calculate the equilibrium equations for the Wilson model regarding pressure.

    Parameters:
    variables (tuple): A tuple containing K1, K2, and pressure.
    temperature (float): The temperature at which to calculate the equilibrium.
    z1 (float): The mole fraction of component 1 in either the vapour or the liquid phase.
    z2 (float): The mole fraction of component 2 in either the vapour or the liquid phase.
    fluid_1 (object): An object representing the first fluid with properties P_c, omega, and T_c.
    fluid_2 (object): An object representing the second fluid with properties P_c, omega, and T_c.
    is_dew (bool): A flag indicating whether to calculate dew point (True) or bubble point (False).

    Returns:
    tuple: A tuple containing the calculated equilibrium equations (eq1, eq2, eq3), where:
        eq1 is the equilibrium equation for component 1,
        eq2 is the equilibrium equation for component 2,
        eq3 is the overall mass balance equation.
    """
    K1, K2, pressure = variables
    eq1 = K1 - fluid_1.P_c * 10**(7/3 * (1 + fluid_1.omega) * (1 - fluid_1.T_c / temperature)) / pressure
    eq2 = K2 - fluid_2.P_c * 10**(7/3 * (1 + fluid_2.omega) * (1 - fluid_2.T_c / temperature)) / pressure
    if is_dew:
        eq3 = z1 / K1 + z2 / K2 - 1
    else:
        eq3 = z1 * K1 + z2 * K2 - 1
    return (eq1, eq2, eq3)
    
# Solving for Wilson approximation values
def Wilson_P_solve(temperature, fluid_1, fluid_2, wilson_pressure=1e5,**kwargs):

    print("Solving Wilson approximation ...")

    K1_guess, K2_guess, pressure_guess = 1, 1, 1e5

    while pressure_guess<1e8:
        try:
            #dew
            wilson_dew = optimize.root(Wilson_P, args = (temperature, 0, 1, fluid_1, fluid_2, True), method="broyden1", x0=(K1_guess, K2_guess, pressure_guess)).x

            #bubble
            wilson_bubble = optimize.root(Wilson_P, args = (temperature, 0, 1, fluid_1, fluid_2, False), method="broyden1", x0=(K1_guess, K2_guess, pressure_guess)).x
            
            # show both results to user
            print(f"Wilson dew: {wilson_dew}, Wilson bubble: {wilson_bubble}")

            if abs(wilson_bubble[2]-wilson_dew[2])/wilson_bubble[2]>1 or wilson_dew[0]<0 or wilson_dew[1]<0 or wilson_dew[2]<0 or wilson_bubble[0]<0 or wilson_bubble[1]<0 or wilson_bubble[2]<0:
                raise ValueError("Wilson approximation has failed due to pressure_guess magnitude")
            
        except ValueError:
            pressure_guess = pressure_guess*1e1
            print(f"Increasing magnitude of pressure guess, pressure_guess = {pressure_guess/1e5} [bar]")

        else:
            K1_wilson, K2_wilson, pressure_wilson = wilson_bubble
            print(f"1st guess [K1, K2, pressure] = {K1_wilson}, {K2_wilson}, {pressure_wilson}")
            break
    
    return K1_wilson, K2_wilson, pressure_wilson

# Wilson approximation root finding problem for dew and bubble point temperature
def Wilson_T(variables, pressure, z1, z2, fluid_1, fluid_2, is_dew=True):
    """
    Calculate the equilibrium equations for the Wilson model regarding temperature.

    Parameters:
    variables (tuple): A tuple containing K1, K2, and temperature.
    pressure (float): The pressure at which the equilibrium is calculated.
    z1 (float): The mole fraction of component 1 in either the vapour or the liquid phase.
    z2 (float): The mole fraction of component 2 in either the vapour or the liquid phase.
    fluid_1 (object): An object representing the first fluid with properties P_c, omega, and T_c.
    fluid_2 (object): An object representing the second fluid with properties P_c, omega, and T_c.
    is_dew (bool): A flag indicating whether to calculate dew point (True) or bubble point (False).

    Returns:
    tuple: A tuple containing the calculated equilibrium equations (eq1, eq2, eq3), where:
           eq1 is the equilibrium equation for component 1,
           eq2 is the equilibrium equation for component 2,
           eq3 is the overall mass balance equation.
    """
    K1, K2, temperature = variables
    eq1 = K1 - fluid_1.P_c * 10**(7/3 * (1 + fluid_1.omega) * (1 - fluid_1.T_c / temperature)) / pressure
    eq2 = K2 - fluid_2.P_c * 10**(7/3 * (1 + fluid_2.omega) * (1 - fluid_2.T_c / temperature)) / pressure
    if is_dew:
        eq3 = z1 / K1 + z2 / K2 - 1
    else:
        eq3 = z1 * K1 + z2 * K2 - 1
    return (eq1, eq2, eq3)

############################################################################### Density root finding
# density root finding, multi-root
def calc_density_multi_root(EOS, args: tuple = (), bracket: list = [2, 5],  n: int = 50):
    f = EOS.pressure_equation
    roots = multi_root(f=f, bracket=bracket, args=args, n=n)
    # print(roots)
    if len(roots) == 0:
        raise ValueError(f"No roots found in the specified bracket, bracket = {bracket} and n = {n}, with args = {args}")
    return roots

# density root finding, cubic polynomial coefficients
def calc_density_coefficients(EOS, pressure: float = None, temperature: float = None, z1: float = None, z2: float = None, coefficient_function: callable = None):
    a_function = EOS.alpha_r.a_function
    b_function = EOS.alpha_r.b_function
    
    coeffs = coefficient_function(a_function, b_function, pressure, temperature, z1, z2)
    roots = np.roots(coeffs)
    # print(roots)

    roots = np.sort(roots[~np.iscomplex(roots)].astype(float))  # Remove complex roots
    roots = roots[roots > 0]  # Keep only positive roots
    if len(roots) == 0:
        raise ValueError("No positive real roots found.")
    return roots

#calculate vapour density
def calc_density_vap(method, EOS, pressure, temperature, y1, y2, **kwargs):
    # multi root method
    if method == "multi-root":
        density = calc_density_multi_root(EOS, (pressure, temperature, y1, y2), **kwargs)[0]

    # cubic polynomial method
    if method == "coefficients":
        density = calc_density_coefficients(EOS, pressure, temperature, y1, y2, **kwargs)[0]

    return density

#calculate liquid density
def calc_density_liq(method, EOS, pressure, temperature, x1, x2, **kwargs):
    # multi root method
    if method == "multi-root":
        density = calc_density_multi_root(EOS, (pressure, temperature, x1, x2), **kwargs)[-1]

    # cubic polynomial method
    if method == "coefficients":
        density = calc_density_coefficients(EOS, pressure, temperature, x1, x2, **kwargs)[-1]

    return density
    
############################################################################### VLE algorithms
# VLE: dew pressure algorithm (watchout for function namespaces)
def VLE_P_dew(EOS, temperature, method= "multi-root", tolerance=1e-6, step=1e-3, pressure_factor=0.2, pressure_break=1e10, **kwargs):
    """
    Calculates the dew point pressure and compositions of a mixture using an equation of state (EOS).

    This function iteratively solves for the vapor-liquid equilibrium (VLE) of a mixture at a given temperature.
    It adjusts the compositions and pressures until the system reaches equilibrium, defined by the specified tolerance.

    Parameters:
    mix_EOS (object): An object representing the mixture's equation of state, which provides methods for calculating fugacity coefficients and pressure.
    temperature (float): The temperature at which the VLE calculations are performed (in Kelvin).
    tolerance (float, optional): The convergence tolerance for the sum of mole fractions (default is 1e-6).
    step (float, optional): The step size for adjusting the vapor composition (default is 1e-3).
    density_brackets (list of tuples, optional): Brackets for density calculations (default is [(1e1, 1e5), (1e2, 1e5)]).
    n (int, optional): Number of points for density calculations (default is 1000).
    pressure_jump (float, optional): The pressure adjustment step for the iterative solver (default is 5e5).
    pressure_break (float, optional): The maximum pressure limit for the calculations (default is 1e10).

    Returns:
    tuple: A tuple containing two lists:
        - compositions (list): The calculated mole fractions of the liquid phase at equilibrium.
        - pressures (list): The corresponding pressures at which the compositions were calculated.
    """

    
    def loop_summation_change(x_t, x, pressure, phi_vap, density_vap, density_liq, **kwargs):
        """
        Adjusts the composition vector `x` until the sum of its elements converges 
        to a specified tolerance level. The function iteratively normalizes `x`, 
        recalculates liquid density and fugacity, and updates `x` based on the 
        calculated equilibrium constants.

        Parameters:
        x_t (float): The initial sum of the composition vector `x`.
        x (list of float): The composition vector representing the mole fractions.
        pressure (float): The pressure at which the calculations are performed.
        phi_vap (list of float): The vapor fugacity coefficients.

        Returns:
        tuple: A tuple containing the final sum of `x` and the updated composition vector `x`.
        """
        # print("tiny loop 1")
        #has the sum of x changed
        while abs(x_t - sum(x)) > tolerance:
            #recalculate the sum
            x_t = sum(x)

            #normalise x
            x = x / sum(x)
            #calculate new liquid density
            density_liq = calc_density_liq(method, EOS, pressure, temperature, *x,  **kwargs)
            # print("densities = ", density_vap, density_liq)

            #calculate new liquid fugacity
            phi_liq = [phi_1(density_liq, temperature, *x), phi_2(density_liq, temperature, *x)]
            # print("fugacities = ", phi_vap, phi_liq)

            #recalculate K
            K = [phi_liq[i] / phi_vap[i] for i in range(2)]
            # print("loop_summation_change namespace:", phi_vap)
            #recalculate x from y and K
            x = [y[i] / K[i] for i in range(2)]   
            # print("x = ", x)
            # print("sum of x = ", sum(x))

        return x_t, x, density_vap, density_liq

    def loop_mole_fraction_1(x_t, x, pressure, density_vap, density_liq, **kwargs):
        """
        Adjusts the mole fractions of a mixture until their sum is approximately equal to 1.

        This function iteratively modifies the pressure based on the sum of the mole fractions
        and recalculates the densities and fugacity coefficients until the sum of the mole fractions
        is within a specified tolerance.

        Parameters:
        x_t (float): The initial sum of mole fractions.
        x (list): A list of mole fractions for the liquid phase.
        pressure (float): The pressure at which equilibrium is being calculated.

        Returns:
        tuple: A tuple containing the updated mole fractions (list), 
               the adjusted pressure (float), 
               the vapor density (float), 
               and the liquid density (float).
        """
        # print("tiny loop 2")
        # is the sum of x 1
        while abs(sum(x) - 1) > tolerance:
            # recalculate sum of x
            x_t = sum(x)

            # if higher than 1 then decrease pressure (too much liquid)
            if sum(x) > 1:
                pressure = pressure - pressure*pressure_factor * min((1, abs(sum(x) - 1)))

            # else, if less than 1, then increase pressure (too little liquid)
            elif sum(x) < 1:
                pressure = pressure + pressure*pressure_factor * min((1, abs(sum(x) - 1)))
            # print("pressure = ", pressure)

            # normalise x
            x = x/sum(x)

            # recalculate densities
            # calculate densities
            density_vap = calc_density_vap(method, EOS, pressure, temperature, *y, **kwargs)
            density_liq = calc_density_liq(method, EOS, pressure, temperature, *x, **kwargs)
            # print("densities = ", density_vap, density_liq)

            # calculate fugacity coefficients
            phi_vap = [phi_1(density_vap, temperature, *y), phi_2(density_vap, temperature, *y)]
            phi_liq = [phi_1(density_liq, temperature, *x), phi_2(density_liq, temperature, *x)]
            # print("fugacities = ", phi_vap, phi_liq)

            # calculate K from fugacity
            K = [phi_liq[i] / phi_vap[i] for i in range(2)]
            # print("K = ", K)

            # recalculate x from y and K
            x = [y[i] / K[i] for i in range(2)]
            # print("x = ", x)
            # print("loop_mole_fraction_1 namespace:", phi_vap)

            x_t, x, density_vap, density_liq = loop_summation_change(x_t, x, pressure, phi_vap, density_vap, density_liq, **kwargs)

        # print(density_vap, density_liq)

        return x, pressure, density_vap, density_liq


    # mixture
    mix = EOS.mixture

    #fugacity coefficient functions
    fugacity_coefficients = EOS.fugacity_coefficients
    phi_1 = sp.lambdify((rho, T, z1, z2), fugacity_coefficients[0])
    phi_2 = sp.lambdify((rho, T, z1, z2), fugacity_coefficients[1])


    print("")
    print("*"*50, " Dew solver ", "*"*50)
    print("")

    # stores for compositions and pressures
    compositions = []
    pressures = []

    #initialise pressure and y
    pressure = 0
    y = [0, 1]

    # 1st guess for solver from relevant Wilson approximation
    K1_guess, K2_guess, pressure_guess = Wilson_P_solve(temperature, mix.fluid_1, mix.fluid_2, **kwargs)
    

    print("")
    print("y, Pressure, Liquid Density, Vapour Density")
    # composition loop, break if y[0] reaches 1 or pressure hits specified limit
    while pressure < pressure_break and y[0] < 1:
        # print("main loop")
        #define vapour composition

        #calculate a first guess for K1, K2 and pressure from the Wilson approximation
        K1, K2, pressure = K1_guess, K2_guess, pressure_guess

        #calculate x from y and K
        x = [y[0]/K1, y[1]/K2]

        #sum of x
        x_t= sum(x)
        # print("x_t = ", x_t)

        # normalise x
        # x = x/sum(x)

        # print("x on main loop before fugacity = ", x)

        #calculate densities
        density_vap = calc_density_vap(method, EOS, pressure, temperature, *y, **kwargs)
        density_liq = calc_density_liq(method, EOS, pressure, temperature, *x, **kwargs)
        # print("densities on main loop = ", density_vap, density_liq)

        #calculate gugacity coefficients
        phi_vap = [phi_1(density_vap, temperature, *y), phi_2(density_vap, temperature, *y)]
        phi_liq = [phi_1(density_liq, temperature, *x), phi_2(density_liq, temperature, *x)]
        # print("fugacities on main loop = ", phi_vap, phi_liq)

        #calculate K from fugacity
        K = [phi_liq[i]/phi_vap[i] for i in range(2)]
        # print("K on main loop= ", K)

        #recalculate x from y and K
        x = [y[i]/K[i] for i in range(2)]
        # print("x on main loop = ",x)

        # ====================================== LOOPS
        # print("main namespace 1:", density_liq, density_vap)

        # 1st summation loop
        x_t, x, density_vap, density_liq = loop_summation_change(x_t, x, pressure, phi_vap, density_vap, density_liq, **kwargs)
        # print("main namespace 2:", density_liq, density_vap)

        x, pressure, density_vap, density_liq = loop_mole_fraction_1(x_t, x, pressure, density_vap, density_liq, **kwargs)
        # print("main namespace 3:", density_liq, density_vap)


        # ====================================== Results
        print(f"[{y[0]:.4f}, {y[1]:.4f}]: ", f"{pressure:.2f}", f"{density_liq:.2f}", f"{density_vap:.2f}")
        pressures += [pressure]
        compositions += [y[0]]

        #update the first_guess
        if x[0]==0 or np.isnan(x[0]):
            K1_guess=1
        else:
            K1_guess = y[0]/x[0]
        
        if x[1]==0 or np.isnan(x[1]):
            K2_guess=1
        else:
            K2_guess = y[1]/x[1]

        # # pressure_guess dependent on composition step
        if len(pressures)<2:
            y = [y[0]+step, y[1]-step]
            pressure_guess = pressure
        else:
            pressure_guess = pressure + (pressures[-2]- pressures[-1])/(compositions[-2] - compositions[-1])*step
            y = [y[0]+step, y[1]-step]

        # disable if needed
        # y = [y[0]+step, y[1]-step]
        # pressure_guess = pressure

    return compositions, pressures

# VLE: dew pressure algorithm (watchout for function namespaces)
def VLE_P_bubble(EOS, temperature, method= "multi-root", tolerance=1e-6, step=1e-3, pressure_factor=0.2, pressure_break=1e10, **kwargs):
    """
    Calculates the bubble point pressure and compositions of a mixture using an equation of state (EOS).

    This function iteratively solves for the vapor-liquid equilibrium (VLE) of a mixture at a given temperature.
    It adjusts the compositions and pressures until the system reaches equilibrium, defined by the specified tolerance.

    Parameters:
    mix_EOS (object): An object representing the mixture's equation of state, which provides methods for calculating fugacity coefficients and pressure.
    temperature (float): The temperature at which the VLE calculations are performed (in [K]).
    tolerance (float, optional): The convergence tolerance for the sum of mole fractions (default is 1e-6).
    step (float, optional): The step size for adjusting the vapor composition (default is 1e-3).
    density_brackets (list of tuples, optional): Brackets for density calculations (default is [(1e1, 1e5), (1e2, 1e5)]).
    n (int, optional): Number of points for density calculations (default is 1000).
    pressure_jump (float, optional): The pressure adjustment step for the iterative solver (default is 5e5 [Pa]).
    pressure_break (float, optional): The maximum pressure limit for the calculations (default is 1e10 [Pa]).

    Returns:
    tuple: A tuple containing two lists:
        - compositions (list): The calculated mole fractions of the liquid phase at equilibrium.
        - pressures (list): The corresponding pressures at which the compositions were calculated.
    """

    
    def loop_summation_change(y_t, y, pressure, density_vap, density_liq, **kwargs):
        """
        Adjusts the composition vector `y` until the sum of its elements converges 
        to a specified tolerance level. The function iteratively normalizes `y`, 
        recalculates liquid density and fugacity, and updates `y` based on the 
        calculated equilibrium constants.

        Parameters:
        y_t (float): The initial sum of the composition vector `y`.
        y (list of float): The composition vector representing the vapour mole fractions.
        pressure (float): The pressure at which the calculations are performed.
        phi_vap (list of float): The vapor fugacity coefficients.
        density_vap (float): vapour density
        density_liq (float): liquid density

        Returns:
        tuple: A tuple containing the final sum of `y` and the updated composition vector `y`.
        """
        # print("tiny loop 1")
        #has the sum of y changed
        while abs(y_t - sum(y)) > tolerance:
            #recalculte the sum
            y_t = sum(y)

            #normalise y
            y = y / sum(y)
            #calculate new liquid density
            density_liq = calc_density_liq(method, EOS, pressure, temperature, *x, **kwargs)

            #calculate new liquid fugacity
            phi_vap = [phi_1(density_vap, temperature, *y), phi_2(density_vap, temperature, *y)]
            phi_liq = [phi_1(density_liq, temperature, *x), phi_2(density_liq, temperature, *x)]

            #recalculate K
            K = [phi_liq[i] / phi_vap[i] for i in range(2)]
            # print("K = ", K)
            #recalculate y from x and K
            y = [x[i] * K[i] for i in range(2)]   
            # print("y = ", y)
            # print("sum of y = ", sum(xy))


        return y_t, y, density_vap, density_liq

    def loop_mole_fraction_1(y_t, y, pressure, density_vap, density_liq, **kwargs):
        """
        Adjusts the vapour mole fractions of a mixture until their sum is approximately equal to 1.

        This function iteratively modifies the pressure based on the sum of the mole fractions
        and recalculates the densities and fugacity coefficients until the sum of the mole fractions
        is within a specified tolerance.

        Parameters:
        y_t (float): The initial sum of vapour mole fractions.
        y (list): A list of mole fractions for the liquid phase.
        pressure (float): The pressure at which equilibrium is being calculated.
        density_vap (float): The vapour density
        density_liq (float): The liquid density

        Returns:
        tuple: A tuple containing the updated mole fractions (list), 
               the adjusted pressure (float), 
               the vapor density (float), 
               and the liquid density (float).
        """
        # print("tiny loop 2")
        # is the sum of y 1
        while abs(sum(y) - 1) > tolerance:
            # recalculate sum of y
            y_t = sum(y)

            # if higher than 1 then increase pressure (too much vapour)
            if sum(y) > 1:
                pressure = pressure + pressure*pressure_factor * min((1, abs(sum(y) - 1)))

            # else, if less than 1, then decerease pressure (too little vapour)
            elif sum(y) < 1:
                pressure = pressure - pressure*pressure_factor * min((1, abs(sum(y) - 1)))

            
            # print("pressure = ", pressure)
            # print("y = ", y)

            # normalise y
            y = y / sum(y)

            # recalculate densities
            density_vap = calc_density_vap(method, EOS, pressure, temperature, *y, **kwargs)
            density_liq = calc_density_liq(method, EOS, pressure, temperature, *x, **kwargs)
            # print("densities = ",density_vap, density_liq)

            # calculate fugacity coefficients
            phi_vap = [phi_1(density_vap, temperature, *y), phi_2(density_vap, temperature, *y)]
            phi_liq = [phi_1(density_liq, temperature, *x), phi_2(density_liq, temperature, *x)]
            # print("fugacities = ", phi_vap, phi_liq)

            # calculate K from fugacity
            K = [phi_liq[i] / phi_vap[i] for i in range(2)]
            # print("K = ", K)

            # recalculate y from x and K
            y = [x[i] * K[i] for i in range(2)]
            # print("y before summation change loop: ", y)
            # print("loop_mole_fraction_1 namespace:", phi_vap)

            y_t, y, density_vap, density_liq = loop_summation_change(y_t, y, pressure, density_vap, density_liq, **kwargs)

            # print("y = ", y)
        # print(density_vap, density_liq)

        return y, pressure, density_vap, density_liq

    # density root finding method
    # mixture
    mix = EOS.mixture

    #fugacity coefficient functions
    fugacity_coefficients = EOS.fugacity_coefficients
    phi_1 = sp.lambdify((rho, T, z1, z2), fugacity_coefficients[0])
    phi_2 = sp.lambdify((rho, T, z1, z2), fugacity_coefficients[1])



    print("")
    print("*"*50, " Bubble solver ", "*"*50)
    print("")

    # stores for compositions and pressures
    compositions = []
    pressures = []

    #initialise pressure and y
    pressure = 0
    x = [0, 1]

    # 1st guess for solver from relevant Wilson approximation
    K1_guess, K2_guess, pressure_guess = Wilson_P_solve(temperature, mix.fluid_1, mix.fluid_2, **kwargs)

    print("")
    print("y, Pressure, Liquid Density, Vapour Density")
    # composition loop, break if x[0] reaches 1 or pressure hits specified limit
    while pressure < pressure_break and x[0] < 1:
        # print("main loop")
        #define vapour composition

        #calculate a first guess for K1, K2 and pressure from the Wilson approximation
        K1, K2, pressure = K1_guess, K2_guess, pressure_guess

        #calculate y from x and K
        y = [x[0]*K1, x[1]*K2]

        #sum of y
        y_t = sum(y)
        # print("y_t = ", y_t)

        #calculate densities
        density_vap = calc_density_vap(method, EOS, pressure, temperature, *y, **kwargs)
        density_liq = calc_density_liq(method, EOS, pressure, temperature, *x, **kwargs)
        # print(density_vap, density_liq)

        #calculate gugacity coefficients
        phi_vap = [phi_1(density_vap, temperature, *y), phi_2(density_vap, temperature, *y)]
        phi_liq = [phi_1(density_liq, temperature, *x), phi_2(density_liq, temperature, *x)]
        # print("fugacities = ", phi_vap, phi_liq)

        #calculate K from fugacity
        K = [phi_liq[i]/phi_vap[i] for i in range(2)]
        # print("K = ", K)

        #recalculate x from x and K
        y = [x[i]*K[i] for i in range(2)]
        # print("y = ",y)

        # ====================================== LOOPS
        # print("main namespace 1:", density_liq, density_vap)

        # 1st summation loop
        y_t, y, density_vap, density_liq = loop_summation_change(y_t, y, pressure, density_vap, density_liq,**kwargs)
        # print("main namespace 2:", density_liq, density_vap)

        y, pressure, density_vap, density_liq = loop_mole_fraction_1(y_t, y, pressure, density_vap, density_liq, **kwargs)
        # print("main namespace 3:", density_liq, density_vap)


        # ====================================== Results
        print(f"[{x[0]:.4f}, {x[1]:.4f}]: ", f"{pressure:.2f}", f"{density_liq:.2f}", f"{density_vap:.2f}")
        pressures += [pressure]
        compositions += [x[0]]

        #update the first_guess
        if y[0]==0 or np.isnan(y[0]):
            K1_guess=1
        else:
            K1_guess = y[0]/x[0]
        
        if y[1]==0 or np.isnan(y[1]):
            K2_guess = 1
        else:
            K2_guess = y[1]/x[1]
        # print(f"x = {x}, and y = {y}")


        # # pressure_guess dependent on composition step
        # if len(pressures)<2:
        #     y = [x[0]+step, x[1]-step]
        #     pressure_guess = pressure
        # else:
        #     pressure_guess = pressure + (pressures[-2]- pressures[-1])/(compositions[-2] - compositions[-1])*step
        #     x = [x[0]+step, x[1]-step]

        # disable if needed
        x = [x[0]+step, x[1]-step]
        pressure_guess = pressure

        # print("guess = ", K1_guess, K2_guess, pressure_guess)

    return compositions, pressures


