#13/01/2025

#other useful functions to use with the EOS and EOS data
import numpy as np
import pandas as pd
from scipy.spatial import distance
from scipy import optimize

from .symbols import *
from .purefluid import PureFluid
from .eos import EOS

# find multiple roots of a function given logarithmically spaced x values (univariate function only)
# recommended over multi_root if used multiple times with the same x
def multi_root_x(f: callable = None, x: np.ndarray = None, args: tuple = (), tol: float = 1) -> np.ndarray:
    """ Find all roots of f in given an x array.

    Fine-grained root finding is performed with `scipy.optimize.root_scalar`.

    Parameters
    ----------
    f: Callable
        Function to be evaluated
    x: np.ndarray
        Array with values at which callable f is evaluated.
    args: Iterable, optional
        Iterable passed to `f` for evaluation
    tol: float = 1
        Absolute tolerance for the discontinuity check of the root

    Returns
    -------
    roots: np.ndarray
        Array containing all unique roots that were found in `bracket`.
    """
    
    # Evaluate function in given x
    y = f(x, *args)

    # Find where adjacent signs are not equal
    sign = np.sign(y)
    sign_changes = np.where(sign[:-1] != sign[1:])[0]

    # find roots around sign changes usign the brentq method
    roots = [optimize.root_scalar(f=f,method ="brentq", args = args, bracket=(x[s], x[s+1]), xtol=1e-12, rtol=1e-12).root for s in sign_changes]
    roots = np.array(roots)     # array conversion needed for discontinuity check

    # check for discontinuities and erase if present
    roots = roots[f(roots, *args)<tol]

    # sort the roots in ascending order
    roots.sort()

    return roots   

# find multiple roots of a function between a given bracket (univariate functions only)
def multi_root(f: callable = None, bracket: list = [], args: tuple = (), n: int = 0) -> np.ndarray:
    """Find all roots within a given bracket

    Uses logarithmic spacing to create an x array where f is evaluated
    The roots are then found between all the points at which f changes sign
    using the brentq method.

    Parameters
    ----------
    f: callable
        The function for which root finding is necessary
    bracket: list
        A list of two elements which specifies the x limits of the root finding
    args: tuple
        Any additional arguments the callable f requires
    n: int
        Amount of x points to create for root finding, the more x points the more
        it is assured all roots will be found, but the less efficient the code is

    Returns
    -------
    roots: np.ndarray
        All the roots found in the specified bracket
    """
    x = np.logspace(*bracket, n)

    roots = multi_root_x(f, x, args)

    return roots

#
def cubic_root_density(eos: EOS, pressure: float, temperature: float) -> np.ndarray:
    if eos.name not in ["vdW", "PR", "RK", "MSRK"]:
        raise AttributeError(f"{eos.name} is not one of the cubic EOS available.")
    
    def vdW_coefficients(a_function, b_function, pressure, temperature):
        #see log page 192
        a = a_function(temperature)
        b = b_function()

        coeff_0 = pressure
        coeff_1 = -pressure * b - R * temperature
        coeff_2 = a
        coeff_3 = - a * b

        return np.array([coeff_3, coeff_2, coeff_1, coeff_0])
    
    def PR_coefficients(a_function, b_function, pressure, temperature):
        #see log page 192
        a = a_function(temperature)
        b = b_function()

        coeff_0 = pressure
        coeff_1 = pressure * b - R * temperature
        coeff_2 = -3*pressure*b**2 - 2*R*temperature*b + a
        coeff_3 = b**3*pressure + R*temperature*b**2 - a * b

        return np.array([coeff_3, coeff_2, coeff_1, coeff_0])

    def RK_coefficients(a_function, b_function, pressure, temperature):
        #see log page 192
        a = a_function(temperature)
        b = b_function()

        coeff_0 = pressure
        coeff_1 = - R * temperature
        coeff_2 = -pressure*b**2 - R*temperature*b + a
        coeff_3 = - a * b

        return np.array([coeff_3, coeff_2, coeff_1, coeff_0])

    def SRK_coefficients(a_function, b_function, pressure, temperature):
        #see log page 192
        a = a_function(temperature)
        b = b_function()

        coeff_0 = pressure
        coeff_1 = -R * temperature
        coeff_2 = -pressure*b**2 - b*R*temperature + a
        coeff_3 = - a * b

        return np.array([coeff_3, coeff_2, coeff_1, coeff_0])

    def MSRK_coefficients(a_function, b_function, pressure, temperature):
        #see log page 192
        a = a_function(temperature)
        b = b_function()

        coeff_0 = pressure
        coeff_1 = -R * temperature
        coeff_2 = -pressure*b**2 - b*R*temperature + a
        coeff_3 = - a * b

        return np.array([coeff_3, coeff_2, coeff_1, coeff_0])

    coefficient_function = eval(f"{eos.name}_coefficients")

    a_function = eos.alpha_r.a_function
    b_function = eos.alpha_r.b_function
    
    coeffs = coefficient_function(a_function, b_function, pressure, temperature)
    roots = np.roots(coeffs)
    # print(roots)

    roots = np.sort(roots[~np.iscomplex(roots)].astype(float))  # Remove complex roots
    roots = roots[roots > 0]  # Keep only positive roots
    if len(roots) == 0:
        raise ValueError("No positive real roots found.")
    return roots

#choose the best root from a list by knowing the experimental value
def choose_root(roots: list, experimental_value=0):
    roots = np.array(roots)

    difference = np.abs(experimental_value-roots)

    try:
        best_index = np.nanargmin(difference)
        best_root = roots[best_index]

    except ValueError:
        best_root = np.nan

    return best_root

#see zandbox for trial and error
def calc_rho_inverse_distance(df_missing_rho, df_experimental, pascals=1e4, neighbours=3):
    #restrict experimental dataframe to Paper, P, T, rho and deep copy it
    df_experimental = df_experimental.copy(deep=True)[["Paper", "P", "T", "rho"]]

    #restrict the other dataframes to Paper, P,T and deep copy it
    df_missing_rho = df_missing_rho.copy(deep=True)[["Paper", "P", "T"]]

    #create og index columns
    df_missing_rho["index"] = df_missing_rho.index
    df_experimental["index"] = df_experimental.index
    
    #merge the two dataframes on pressure and temperature
    df = pd.merge(df_experimental, df_missing_rho, on=["Paper","P", "T"], how="outer")

    #sort by temperature and then by pressure
    df = df.sort_values(by=["T","P"])
    #delete pressure and temperature duplicates
    # df = df.drop_duplicates(subset=["P", "T"])

    #reset the indices
    df = df.reset_index()

    #create existance column
    df["exist"] = df["rho"].notna()

    #reduce the pressure by the pascals
    df["P"] = df["P"]/pascals

    #calculate distances between points
    distances = distance.cdist(df[["T","P"]], df[["T","P"]], "euclidean")

    #multiply the distances by the truth vector
    distances_true = (distances.T * np.array(df["exist"])).T

    #make the distances a dataframe
    df_distances = pd.DataFrame(data=distances_true)

    #get the point numbers for which you would like to linearly interpolate
    point_numbers = df.index[df["exist"]==False]

    #list to keep track of interpolation values
    rho_interp = []

    #for unknown point
    for i in point_numbers:
        distance_column = df_distances[i]

        #get the smallest neighbours nonzero values
        smallest_n = distance_column.loc[distance_column!=0].nsmallest(neighbours)

        #calculate the inverse distance density and add it to the list
        rho_interp +=[sum(df.loc[smallest_n.index, "rho"]*1/smallest_n.values)/sum(1/smallest_n.values)]

    #add the interpolated values to the dataframe
    df.loc[point_numbers,["rho"]] = rho_interp

    #return the pressure back to normal
    df["P"] = df["P"]*pascals

    #drop all the rows where the index_y is Nan
    df = df[df["index_y"].notna()]

    # set index to the index of the dataframe which is missing density
    df = df.set_index("index_y").sort_index().reset_index(drop=True)

    return df[["Paper", "P", "T", "rho"]]   #return the  dataframe but only Paper, P, T, rho

def calc_Psat(fluid: PureFluid, eos: EOS, temperature_list: list, bracket:tuple =(), n: int = int(1e3), debug: bool = False):
    """Calculate the saturation pressure of a fluid at given temperatures.

    Parameters
    ----------
    fluid: PureFluid
        The fluid object containing properties like critical and triple point temperatures.
    eos: EOS
        The equation of state object used for pressure and fugacity calculations.
    temperature_list: list
        A list of temperatures at which to calculate saturation pressures.

    Returns
    -------
    P_sat: list
        A list of calculated saturation pressures corresponding to the input temperatures.
    """
    
    pressure_equation = sp.utilities.lambdify((rho, T, P), P - eos.pressure)
    fugacity_coefficient_function = sp.utilities.lambdify((rho, T), eos.fugacity_coefficient)
    densities = np.logspace(*bracket, n)

    def get_pressure_guess(temperature: float) -> float:
        """Get an initial pressure guess for a specified temperature.

        Parameters
        ----------
        temperature: float
            The temperature for which to estimate the pressure.

        Returns
        -------
        pressure_guess: float
            An initial pressure guess.
        """
        pressure_guess = fluid.P_t# + (fluid.P_c - fluid.P_t) * (temperature - fluid.T_t) / (fluid.T_c - fluid.T_t)
        n_roots = 0
        while n_roots < 2:
            pressure_guess += 1e4
            density_roots = multi_root_x(pressure_equation, densities, args=(temperature, pressure_guess))
            n_roots = len(density_roots)
        return pressure_guess

    def calc_saturation_pressure(pressure: float, temperature: float) -> float:
        """ This is the function for which finding the root yields the saturation pressure

        Parameters
        ----------
        pressure: float
            The initial pressure guess.
        temperature:  float
            The temperature for which to calculate saturation pressure.

        Returns
        -------
        eq: float
            The result of the fugacity ratio equation.
        """
        density_roots = multi_root_x(pressure_equation, densities, args=(temperature, pressure))
        fugacities = fugacity_coefficient_function(density_roots, temperature)
        fugacities = fugacities[fugacities != max(fugacities)]
        eq = fugacities[0] / fugacities[1] - 1
        return eq

    P_sat = []
    pressure_guess = get_pressure_guess(temperature_list[0])

    for temperature in temperature_list:
        if temperature > fluid.T_c or temperature < fluid.T_t:
            print("Temperature is outside the critical and triple point bounds")
            break
        try:
            sol = optimize.root(calc_saturation_pressure, x0=[pressure_guess], args=(temperature,), method='broyden1', options={"xtol": 1e2})
            solution = sol.x
            pressure_guess = sol.x[0]
        except IndexError:
            if debug:
                print("IndexError: no solution has been found due to the root finding method reaching a point with less than two density roots")
                print("Assigning 0 to return value")
            solution = [0]
        P_sat += list(solution)

        if debug:
            print(f"P_sat({temperature} K) = {solution} Pa")

    return P_sat


